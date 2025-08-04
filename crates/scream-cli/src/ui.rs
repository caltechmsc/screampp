use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressState, ProgressStyle};
use screampp::engine::progress::Progress;
use std::sync::Arc;
use std::time::Duration;
use tokio::sync::{mpsc, watch};
use tracing::warn;

#[derive(Debug)]
pub enum UiEvent {
    Progress(Progress),
    Log(String),
}

pub struct UiManager {
    mp: Arc<MultiProgress>,
    state: BarState,
    event_receiver: mpsc::Receiver<UiEvent>,
    shutdown_receiver: watch::Receiver<bool>,
    _sentinel_bar: ProgressBar,
}

#[derive(Default)]
struct BarState {
    active_bar: Option<ProgressBar>,
    base_message: String,
}

impl UiManager {
    pub fn new() -> (Self, mpsc::Sender<UiEvent>, watch::Sender<bool>) {
        let (event_sender, event_receiver) = mpsc::channel(1024);
        let (shutdown_sender, shutdown_receiver) = watch::channel(false);
        let mp = Arc::new(MultiProgress::new());
        mp.set_draw_target(ProgressDrawTarget::stderr_with_hz(12));
        let _sentinel_bar = mp.add(ProgressBar::hidden());
        let manager = Self {
            mp,
            state: BarState::default(),
            event_receiver,
            shutdown_receiver,
            _sentinel_bar,
        };

        (manager, event_sender, shutdown_sender)
    }

    pub async fn run(mut self) {
        loop {
            tokio::select! {
                Some(event) = self.event_receiver.recv() => {
                    self.handle_event(event);
                }
                result = self.shutdown_receiver.changed() => {
                    if result.is_err() || *self.shutdown_receiver.borrow() {
                        break;
                    }
                }
            }
        }
        self._sentinel_bar.finish_and_clear();
    }

    fn handle_event(&mut self, event: UiEvent) {
        match event {
            UiEvent::Log(msg) => {
                self.mp.println(msg).ok();
            }
            UiEvent::Progress(progress) => self.handle_progress(progress),
        }
    }

    fn handle_progress(&mut self, progress: Progress) {
        match progress {
            Progress::PhaseStart { name } => {
                if let Some(bar) = self.state.active_bar.take() {
                    bar.finish_and_clear();
                }

                let pb = self.mp.add(ProgressBar::new_spinner());
                pb.enable_steady_tick(Duration::from_millis(80));
                pb.set_style(Self::spinner_style());
                pb.set_message(name.to_string());

                self.state.base_message = name.to_string();
                self.state.active_bar = Some(pb);
            }
            Progress::PhaseFinish => {
                if let Some(bar) = self.state.active_bar.take() {
                    bar.finish_and_clear();
                }

                let final_message = format!("✓ {}", self.state.base_message);
                self.mp.println(final_message).ok();

                self.state.base_message.clear();
            }
            Progress::TaskStart { total } => {
                if let Some(bar) = self.state.active_bar.as_ref() {
                    bar.reset();
                    bar.set_style(Self::bar_style());
                    bar.set_length(total);
                    bar.set_message(self.state.base_message.clone());
                    bar.disable_steady_tick();
                }
            }
            Progress::TaskIncrement { amount } => {
                if let Some(bar) = self.state.active_bar.as_ref() {
                    bar.inc(amount);
                }
            }
            Progress::TaskFinish => {
                if let Some(bar) = self.state.active_bar.as_ref() {
                    bar.finish();
                }
            }
            Progress::StatusUpdate { text } => {
                if let Some(bar) = self.state.active_bar.as_ref() {
                    bar.set_message(format!("{} ({})", self.state.base_message, text));
                }
            }
            Progress::Message(msg) => {
                self.mp.println(format!("  {}", msg)).ok();
            }
        }
    }

    fn spinner_style() -> ProgressStyle {
        ProgressStyle::with_template("{spinner:.green} {msg}")
            .expect("Invalid template")
            .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
    }

    fn bar_style() -> ProgressStyle {
        ProgressStyle::with_template("{msg:<45} [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .expect("Invalid template")
            .with_key(
                "eta",
                |state: &ProgressState, w: &mut dyn std::fmt::Write| {
                    write!(w, "{:.1}s", state.eta().as_secs_f64()).unwrap();
                },
            )
            .progress_chars("━╸ ")
    }
}

#[derive(Clone)]
pub struct CliProgressHandler {
    sender: mpsc::Sender<UiEvent>,
}

impl CliProgressHandler {
    pub fn new(sender: mpsc::Sender<UiEvent>) -> Self {
        Self { sender }
    }

    pub fn get_callback(&self) -> screampp::engine::progress::ProgressCallback<'static> {
        let sender = self.sender.clone();
        Box::new(move |progress: Progress| {
            if let Err(e) = sender.try_send(UiEvent::Progress(progress)) {
                warn!("Failed to send progress update to UI channel: {}", e);
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tokio::sync::mpsc;

    #[tokio::test]
    async fn cli_progress_handler_sends_progress_updates() {
        let (sender, mut receiver) = mpsc::channel(10);
        let handler = CliProgressHandler::new(sender);
        let callback = handler.get_callback();
        let progress = Progress::PhaseStart { name: "Test" };

        callback(progress.clone());

        let received = receiver.recv().await.unwrap();
        match received {
            UiEvent::Progress(p) => assert_eq!(p, progress),
            _ => panic!("Incorrect event type received"),
        }
    }

    #[tokio::test]
    async fn ui_manager_handles_log_event() {
        let (ui_manager, sender, _shutdown_sender) = UiManager::new();
        let mp = ui_manager.mp.clone();
        mp.set_draw_target(ProgressDrawTarget::hidden());

        let log_message = "Test log message".to_string();
        sender
            .send(UiEvent::Log(log_message.clone()))
            .await
            .unwrap();
    }

    #[tokio::test]
    async fn ui_manager_handles_phase_start_and_finish() {
        let (mut ui_manager, sender, _shutdown_sender) = UiManager::new();
        ui_manager.mp.set_draw_target(ProgressDrawTarget::hidden());

        sender
            .send(UiEvent::Progress(Progress::PhaseStart { name: "Testing" }))
            .await
            .unwrap();
        let event = ui_manager.event_receiver.recv().await.unwrap();
        ui_manager.handle_event(event);
        assert!(ui_manager.state.active_bar.is_some());
        assert_eq!(ui_manager.state.base_message, "Testing");

        sender
            .send(UiEvent::Progress(Progress::PhaseFinish))
            .await
            .unwrap();
        let event = ui_manager.event_receiver.recv().await.unwrap();
        ui_manager.handle_event(event);
        assert!(ui_manager.state.active_bar.is_none());
        assert!(ui_manager.state.base_message.is_empty());
    }

    #[tokio::test]
    async fn ui_manager_handles_task_lifecycle() {
        let (mut ui_manager, sender, _shutdown_sender) = UiManager::new();
        ui_manager.mp.set_draw_target(ProgressDrawTarget::hidden());

        sender
            .send(UiEvent::Progress(Progress::PhaseStart {
                name: "Task Phase",
            }))
            .await
            .unwrap();
        let event = ui_manager.event_receiver.recv().await.unwrap();
        ui_manager.handle_event(event);

        sender
            .send(UiEvent::Progress(Progress::TaskStart { total: 100 }))
            .await
            .unwrap();
        let event = ui_manager.event_receiver.recv().await.unwrap();
        ui_manager.handle_event(event);
        {
            let bar = ui_manager.state.active_bar.as_ref().unwrap();
            assert_eq!(bar.length(), Some(100));
            assert_eq!(bar.position(), 0);
        }

        sender
            .send(UiEvent::Progress(Progress::TaskIncrement { amount: 10 }))
            .await
            .unwrap();
        let event = ui_manager.event_receiver.recv().await.unwrap();
        ui_manager.handle_event(event);
        {
            let bar = ui_manager.state.active_bar.as_ref().unwrap();
            assert_eq!(bar.position(), 10);
        }

        sender
            .send(UiEvent::Progress(Progress::TaskFinish))
            .await
            .unwrap();
        let event = ui_manager.event_receiver.recv().await.unwrap();
        ui_manager.handle_event(event);
        let bar = ui_manager.state.active_bar.as_ref().unwrap();
        assert!(bar.is_finished());
    }

    #[tokio::test]
    async fn ui_manager_handles_status_update() {
        let (mut ui_manager, sender, _shutdown_sender) = UiManager::new();
        ui_manager.mp.set_draw_target(ProgressDrawTarget::hidden());

        sender
            .send(UiEvent::Progress(Progress::PhaseStart {
                name: "Status Phase",
            }))
            .await
            .unwrap();
        let event = ui_manager.event_receiver.recv().await.unwrap();
        ui_manager.handle_event(event);

        let update_text = "doing something important".to_string();
        sender
            .send(UiEvent::Progress(Progress::StatusUpdate {
                text: update_text.clone(),
            }))
            .await
            .unwrap();
        let event = ui_manager.event_receiver.recv().await.unwrap();
        ui_manager.handle_event(event);

        let bar = ui_manager.state.active_bar.as_ref().unwrap();
        assert_eq!(bar.message(), "Status Phase (doing something important)");
    }

    #[tokio::test]
    async fn ui_manager_shuts_down_on_signal() {
        let (ui_manager, _sender, shutdown_sender) = UiManager::new();
        ui_manager.mp.set_draw_target(ProgressDrawTarget::hidden());

        let handle = tokio::spawn(ui_manager.run());

        shutdown_sender.send(true).unwrap();

        let result = tokio::time::timeout(Duration::from_secs(1), handle).await;
        assert!(result.is_ok(), "UI manager did not shut down in time");
    }
}

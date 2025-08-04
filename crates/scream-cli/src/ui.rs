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

                self.state.active_bar = Some(pb);
                self.state.base_message = name.to_string();
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
                    bar.set_style(Self::bar_style());
                    bar.set_length(total);
                    bar.set_position(0);
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
    use screampp::engine::progress::Progress;
    use tokio::sync::mpsc;

    fn setup_manager() -> (UiManager, mpsc::Sender<UiEvent>) {
        let (manager, sender) = UiManager::new();
        manager.mp.set_draw_target(ProgressDrawTarget::hidden());
        (manager, sender)
    }

    #[test]
    fn handle_phase_start_creates_new_spinner() {
        let (mut manager, _) = setup_manager();
        assert!(manager.state.active_bar.is_none());

        manager.handle_event(UiEvent::Progress(Progress::PhaseStart {
            name: "Test Phase".into(),
        }));

        assert!(manager.state.active_bar.is_some());
        let bar = manager.state.active_bar.as_ref().unwrap();
        assert_eq!(bar.message(), "Test Phase");
        assert_eq!(manager.state.base_message, "Test Phase");
    }

    #[test]
    fn handle_phase_start_replaces_existing_bar() {
        let (mut manager, _) = setup_manager();
        manager.handle_event(UiEvent::Progress(Progress::PhaseStart {
            name: "First Phase".into(),
        }));
        assert_eq!(
            manager.state.active_bar.as_ref().unwrap().message(),
            "First Phase"
        );

        manager.handle_event(UiEvent::Progress(Progress::PhaseStart {
            name: "Second Phase".into(),
        }));

        assert!(manager.state.active_bar.is_some());
        let second_bar = manager.state.active_bar.as_ref().unwrap();
        assert_eq!(second_bar.message(), "Second Phase");
        assert_eq!(manager.state.base_message, "Second Phase");
    }

    #[test]
    fn handle_phase_finish_clears_active_bar() {
        let (mut manager, _) = setup_manager();
        manager.handle_event(UiEvent::Progress(Progress::PhaseStart {
            name: "Test Phase".into(),
        }));
        assert!(manager.state.active_bar.is_some());

        manager.handle_event(UiEvent::Progress(Progress::PhaseFinish));

        assert!(manager.state.active_bar.is_none());
        assert!(manager.state.base_message.is_empty());
    }

    #[test]
    fn handle_task_start_configures_bar_for_task() {
        let (mut manager, _) = setup_manager();
        manager.handle_event(UiEvent::Progress(Progress::PhaseStart {
            name: "Test Phase".into(),
        }));

        manager.handle_event(UiEvent::Progress(Progress::TaskStart { total: 100 }));

        let bar = manager.state.active_bar.as_ref().unwrap();
        assert_eq!(bar.length(), Some(100));
        assert_eq!(bar.position(), 0);
    }

    #[test]
    fn handle_task_increment_updates_bar_position() {
        let (mut manager, _) = setup_manager();
        manager.handle_event(UiEvent::Progress(Progress::PhaseStart {
            name: "Test Phase".into(),
        }));
        manager.handle_event(UiEvent::Progress(Progress::TaskStart { total: 100 }));

        manager.handle_event(UiEvent::Progress(Progress::TaskIncrement { amount: 10 }));

        let bar = manager.state.active_bar.as_ref().unwrap();
        assert_eq!(bar.position(), 10);
    }

    #[test]
    fn handle_task_finish_completes_bar() {
        let (mut manager, _) = setup_manager();
        manager.handle_event(UiEvent::Progress(Progress::PhaseStart {
            name: "Test Phase".into(),
        }));
        manager.handle_event(UiEvent::Progress(Progress::TaskStart { total: 100 }));

        manager.handle_event(UiEvent::Progress(Progress::TaskFinish));

        let bar = manager.state.active_bar.as_ref().unwrap();
        assert!(bar.is_finished());
    }

    #[test]
    fn handle_status_update_changes_bar_message() {
        let (mut manager, _) = setup_manager();
        manager.handle_event(UiEvent::Progress(Progress::PhaseStart {
            name: "Test Phase".into(),
        }));

        manager.handle_event(UiEvent::Progress(Progress::StatusUpdate {
            text: "doing something".into(),
        }));

        let bar = manager.state.active_bar.as_ref().unwrap();
        assert_eq!(bar.message(), "Test Phase (doing something)");
    }

    #[tokio::test]
    async fn cli_progress_handler_sends_progress_event() {
        let (sender, mut receiver) = mpsc::channel(1);
        let handler = CliProgressHandler::new(sender);
        let callback = handler.get_callback();
        let progress = Progress::PhaseStart {
            name: "Testing".into(),
        };

        callback(progress.clone());

        let event = receiver.recv().await.unwrap();
        if let UiEvent::Progress(p) = event {
            if let Progress::PhaseStart { name } = p {
                assert_eq!(name, "Testing");
            } else {
                panic!("Incorrect progress variant received");
            }
        } else {
            panic!("Incorrect event type received");
        }
    }

    #[test]
    fn handle_log_event_prints_message() {
        let (mut manager, _) = setup_manager();
        manager.handle_event(UiEvent::Log("Test log message".to_string()));
    }

    #[test]
    fn handle_progress_message_prints_indented_message() {
        let (mut manager, _) = setup_manager();
        manager.handle_event(UiEvent::Progress(Progress::Message(
            "Test progress message".to_string(),
        )));
    }
}

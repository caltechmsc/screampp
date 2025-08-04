use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressState, ProgressStyle};
use screampp::engine::progress::Progress;
use std::sync::Arc;
use std::time::Duration;
use tokio::sync::mpsc;
use tracing::warn;

#[derive(Debug)]
pub enum UiEvent {
    Progress(Progress),
    Log(String),
}

pub struct UiManager {
    mp: Arc<MultiProgress>,
    state: BarState,
    receiver: mpsc::Receiver<UiEvent>,
}

#[derive(Default)]
struct BarState {
    active_bar: Option<ProgressBar>,
    base_message: String,
}

impl UiManager {
    pub fn new() -> (Self, mpsc::Sender<UiEvent>) {
        let (sender, receiver) = mpsc::channel(1024);
        let mp = Arc::new(MultiProgress::new());
        mp.set_draw_target(ProgressDrawTarget::stderr_with_hz(12));

        let manager = Self {
            mp,
            state: BarState::default(),
            receiver,
        };

        (manager, sender)
    }

    pub async fn run(mut self) {
        while let Some(event) = self.receiver.recv().await {
            self.handle_event(event);
        }
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

use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressState, ProgressStyle};
use screampp::engine::progress::{Progress, ProgressCallback};
use std::sync::{Arc, Mutex};
use std::time::Duration;
use tracing::warn;

#[derive(Default)]
struct BarState {
    active_bar: Option<ProgressBar>,
    base_message: String,
}

#[derive(Clone)]
pub struct CliProgressHandler {
    mp: Arc<MultiProgress>,
    state: Arc<Mutex<BarState>>,
}

impl CliProgressHandler {
    pub fn new() -> Self {
        let mp = MultiProgress::new();
        mp.set_draw_target(ProgressDrawTarget::stderr_with_hz(12));
        Self {
            mp: Arc::new(mp),
            state: Arc::new(Mutex::new(BarState::default())),
        }
    }

    pub fn log(&self, msg: &str) {
        self.mp.println(msg).ok();
    }

    pub fn get_callback(&self) -> ProgressCallback<'static> {
        let mp_clone = self.mp.clone();
        let state_clone = self.state.clone();

        Box::new(move |progress: Progress| {
            let Ok(mut state) = state_clone.lock() else {
                warn!("Progress bar mutex was poisoned; cannot update UI.");
                return;
            };

            match progress {
                Progress::PhaseStart { name } => {
                    if let Some(bar) = state.active_bar.take() {
                        bar.finish_and_clear();
                    }

                    let pb = mp_clone.add(ProgressBar::new_spinner());
                    pb.enable_steady_tick(Duration::from_millis(80));
                    pb.set_style(Self::spinner_style());
                    pb.set_message(name.to_string());

                    state.active_bar = Some(pb);
                    state.base_message = name.to_string();
                }
                Progress::PhaseFinish => {
                    if let Some(bar) = state.active_bar.take() {
                        bar.finish_and_clear();
                    }

                    let final_message = format!("✓ {}", state.base_message);
                    mp_clone.println(final_message).ok();

                    state.base_message.clear();
                }
                Progress::TaskStart { total } => {
                    if let Some(bar) = state.active_bar.as_ref() {
                        bar.set_style(Self::bar_style());
                        bar.set_length(total);
                        bar.set_position(0);
                        bar.disable_steady_tick();
                    }
                }
                Progress::TaskIncrement { amount } => {
                    if let Some(bar) = state.active_bar.as_ref() {
                        bar.inc(amount);
                    }
                }
                Progress::TaskFinish => {
                    if let Some(bar) = state.active_bar.as_ref() {
                        bar.finish();
                        bar.set_style(Self::spinner_style());
                        bar.set_message(state.base_message.clone());
                        bar.enable_steady_tick(Duration::from_millis(80));
                    }
                }
                Progress::StatusUpdate { text } => {
                    if let Some(bar) = state.active_bar.as_ref() {
                        bar.set_message(format!("{} ({})", state.base_message, text));
                    }
                }
                Progress::Message(msg) => {
                    mp_clone.println(format!("  {}", msg)).ok();
                }
            }
        })
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

impl Default for CliProgressHandler {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use screampp::engine::progress::Progress;
    use std::thread;

    #[test]
    fn handler_initializes_in_a_clean_state() {
        let handler = CliProgressHandler::new();
        let pb = handler.pb.lock().unwrap();
        assert_eq!(pb.length(), Some(0));
        assert!(pb.is_finished());
    }

    #[test]
    fn callback_updates_progress_bar_state() {
        let handler = CliProgressHandler::new();
        let callback = handler.get_callback();

        callback(Progress::PhaseStart { name: "Test Phase" });
        {
            let pb = handler.pb.lock().unwrap();
            assert_eq!(pb.message(), "Test Phase");
            assert!(!pb.is_finished());
            assert_eq!(pb.length(), Some(0));
        }

        callback(Progress::TaskStart { total_steps: 100 });
        {
            let pb = handler.pb.lock().unwrap();
            assert_eq!(pb.length(), Some(100));
            assert_eq!(pb.position(), 0);
        }

        callback(Progress::TaskIncrement);
        {
            let pb = handler.pb.lock().unwrap();
            assert_eq!(pb.position(), 1);
        }

        callback(Progress::TaskFinish);
        {
            let pb = handler.pb.lock().unwrap();
            assert!(pb.is_finished());
            assert_eq!(pb.position(), 100);
        }

        callback(Progress::PhaseFinish);
        {
            let pb = handler.pb.lock().unwrap();
            assert_eq!(pb.message(), "✓ Done");
        }
    }

    #[test]
    fn callback_is_thread_safe() {
        let handler = CliProgressHandler::new();
        let callback = handler.get_callback();

        thread::spawn(move || {
            callback(Progress::PhaseStart {
                name: "Thread Test",
            });
            callback(Progress::TaskIncrement);
            callback(Progress::PhaseFinish);
        })
        .join()
        .unwrap();

        let pb = handler.pb.lock().unwrap();
        assert!(pb.is_finished());
        assert_eq!(pb.message(), "✓ Done");
    }
}

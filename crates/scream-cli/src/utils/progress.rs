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

    fn get_active_bar(handler: &CliProgressHandler) -> Option<ProgressBar> {
        handler.state.lock().unwrap().active_bar.clone()
    }

    fn get_base_message(handler: &CliProgressHandler) -> String {
        handler.state.lock().unwrap().base_message.clone()
    }

    #[test]
    fn new_handler_initializes_in_a_clean_state() {
        let handler = CliProgressHandler::new();

        let state = handler.state.lock().unwrap();

        assert!(state.active_bar.is_none());
        assert!(state.base_message.is_empty());
    }

    #[test]
    fn phase_start_creates_a_new_spinner_and_sets_base_message() {
        let handler = CliProgressHandler::new();
        let callback = handler.get_callback();

        callback(Progress::PhaseStart { name: "Test Phase" });

        let bar = get_active_bar(&handler).expect("Bar should be active");
        assert_eq!(bar.message(), "Test Phase");
        assert_eq!(get_base_message(&handler), "Test Phase");
        assert!(!bar.is_finished());
        assert!(!bar.length().is_some_and(|l| l > 0));
    }

    #[test]
    fn task_start_transforms_spinner_into_progress_bar() {
        let handler = CliProgressHandler::new();
        let callback = handler.get_callback();
        callback(Progress::PhaseStart { name: "Processing" });

        callback(Progress::TaskStart { total: 100 });

        let bar = get_active_bar(&handler).expect("Bar should still be active");
        assert_eq!(bar.length(), Some(100));
        assert_eq!(bar.position(), 0);
    }

    #[test]
    fn task_finish_reverts_progress_bar_to_spinner_state() {
        let handler = CliProgressHandler::new();
        let callback = handler.get_callback();
        callback(Progress::PhaseStart { name: "Heavy Task" });
        callback(Progress::TaskStart { total: 100 });
        callback(Progress::TaskIncrement { amount: 100 });

        callback(Progress::TaskFinish);

        let bar = get_active_bar(&handler).expect("Bar should remain active");
        assert!(bar.is_finished());
        assert_eq!(bar.message(), "Heavy Task");
    }

    #[test]
    fn phase_finish_clears_bar_and_prints_permanent_message() {
        let handler = CliProgressHandler::new();
        let callback = handler.get_callback();
        callback(Progress::PhaseStart { name: "Finalizing" });

        let _bar_before_finish = get_active_bar(&handler).unwrap();

        callback(Progress::PhaseFinish);

        assert!(get_active_bar(&handler).is_none());
        assert_eq!(get_base_message(&handler), "");
    }

    #[test]
    fn status_update_combines_base_message_with_new_text() {
        let handler = CliProgressHandler::new();
        let callback = handler.get_callback();
        callback(Progress::PhaseStart {
            name: "Clash Resolution",
        });

        callback(Progress::StatusUpdate {
            text: "Pass 1/10, Clashes: 5".to_string(),
        });

        let bar = get_active_bar(&handler).expect("Bar should be active");
        assert_eq!(bar.message(), "Clash Resolution (Pass 1/10, Clashes: 5)");
        assert_eq!(get_base_message(&handler), "Clash Resolution");
    }

    #[test]
    fn handler_is_thread_safe_and_maintains_state() {
        let handler = CliProgressHandler::new();
        let callback1 = handler.get_callback();
        let callback2 = handler.get_callback();

        let thread1 = thread::spawn(move || {
            callback1(Progress::PhaseStart {
                name: "Parallel Work",
            });
            callback1(Progress::TaskStart { total: 10 });
            for _ in 0..5 {
                callback1(Progress::TaskIncrement { amount: 1 });
                thread::sleep(Duration::from_millis(10));
            }
        });

        let thread2 = thread::spawn(move || {
            for i in 0..3 {
                callback2(Progress::StatusUpdate {
                    text: format!("Update {}", i),
                });
                thread::sleep(Duration::from_millis(15));
            }
        });

        thread1.join().unwrap();
        thread2.join().unwrap();

        let bar = get_active_bar(&handler).expect("Bar should be active");
        assert_eq!(bar.position(), 5);
        assert_eq!(bar.message(), "Parallel Work (Update 2)");
    }
}

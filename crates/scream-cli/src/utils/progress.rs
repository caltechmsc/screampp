use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use screampp::engine::progress::{Progress, ProgressCallback};
use std::sync::{Arc, Mutex};
use std::time::Duration;
use tracing::warn;

const SPINNER_TICK_MS: u64 = 80;

#[derive(Clone)]
pub struct CliProgressHandler {
    pb: Arc<Mutex<ProgressBar>>,
}

impl CliProgressHandler {
    pub fn new() -> Self {
        let pb = ProgressBar::new(0)
            .with_style(Self::spinner_style())
            .with_message("Initializing...");
        pb.set_draw_target(indicatif::ProgressDrawTarget::stderr());
        pb.disable_steady_tick();
        pb.finish_and_clear();

        Self {
            pb: Arc::new(Mutex::new(pb)),
        }
    }

    pub fn get_callback(&self) -> ProgressCallback<'static> {
        let pb_clone = self.pb.clone();

        Box::new(move |progress: Progress| {
            let Ok(mut pb_guard) = pb_clone.lock() else {
                warn!("Progress bar mutex was poisoned. Cannot update progress.");
                return;
            };

            match progress {
                Progress::PhaseStart { name } => {
                    pb_guard.reset();
                    pb_guard.set_length(0);
                    pb_guard.set_style(Self::spinner_style());
                    pb_guard.enable_steady_tick(Duration::from_millis(SPINNER_TICK_MS));
                    pb_guard.set_message(name.to_string());
                }
                Progress::PhaseFinish => {
                    pb_guard.disable_steady_tick();
                    pb_guard.finish_with_message("✓ Done");
                }
                Progress::TaskStart { total_steps } => {
                    pb_guard.disable_steady_tick();
                    pb_guard.reset();
                    pb_guard.set_length(total_steps);
                    pb_guard.set_position(0);
                    pb_guard.set_style(Self::bar_style());
                }
                Progress::TaskIncrement => {
                    pb_guard.inc(1);
                }
                Progress::TaskFinish => {
                    if pb_guard.position() < pb_guard.length().unwrap_or(0) {
                        pb_guard.set_position(pb_guard.length().unwrap_or(0));
                    }
                    pb_guard.finish();
                }
                Progress::Message(msg) => {
                    if !pb_guard.is_finished() {
                        pb_guard.println(format!("  {}", msg));
                    } else {
                        pb_guard.set_message(msg);
                    }
                }
            }
        })
    }

    fn spinner_style() -> ProgressStyle {
        ProgressStyle::with_template("{spinner:.green} {msg}")
            .expect("Failed to create spinner style template")
    }

    fn bar_style() -> ProgressStyle {
        ProgressStyle::with_template("{msg:<20} [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .expect("Failed to create bar style template")
            .with_key(
                "eta",
                |state: &ProgressState, w: &mut dyn std::fmt::Write| {
                    write!(w, "{:.1}s", state.eta().as_secs_f64()).unwrap()
                },
            )
            .progress_chars("##-")
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

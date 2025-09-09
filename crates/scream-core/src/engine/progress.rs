/// Represents different types of progress events during SCREAM++ operations.
///
/// This enum defines the various progress reporting events that can be emitted
/// during molecular optimization and analysis operations. It allows external
/// code to track the progress of long-running computations and provide user feedback.
#[derive(Debug, Clone, PartialEq)]
pub enum Progress {
    /// Signals the start of a new computational phase.
    ///
    /// # Arguments
    ///
    /// * `name` - A static string identifier for the phase being started.
    PhaseStart { name: &'static str },
    /// Signals the completion of the current computational phase.
    PhaseFinish,

    /// Signals the start of a task with a known total amount of work.
    ///
    /// # Arguments
    ///
    /// * `total` - The total number of work units expected for this task.
    TaskStart { total: u64 },
    /// Signals incremental progress on the current task.
    ///
    /// # Arguments
    ///
    /// * `amount` - The number of work units completed in this increment.
    TaskIncrement { amount: u64 },
    /// Signals the completion of the current task.
    TaskFinish,

    /// Provides a status update message for the current operation.
    ///
    /// # Arguments
    ///
    /// * `text` - A descriptive message about the current status.
    StatusUpdate { text: String },

    /// Provides a general informational message.
    ///
    /// The message text to display.
    Message(String),
}

/// A callback function type for handling progress events.
///
/// This type alias defines the signature for functions that can receive and
/// process progress events. The callback must be thread-safe (Send + Sync) and
/// can capture variables from its environment.
pub type ProgressCallback<'a> = Box<dyn Fn(Progress) + Send + Sync + 'a>;

/// Provides a mechanism for reporting progress events during computations.
///
/// This struct manages progress reporting for SCREAM++ operations. It can optionally
/// hold a callback function that will be invoked whenever progress events are reported.
/// If no callback is provided, progress events are silently ignored, allowing the same
/// code to work with or without progress monitoring.
#[derive(Default)]
pub struct ProgressReporter<'a> {
    /// Optional callback function to handle progress events.
    callback: Option<ProgressCallback<'a>>,
}

impl<'a> ProgressReporter<'a> {
    /// Creates a new progress reporter with no callback.
    ///
    /// This constructor creates a reporter that will silently ignore all
    /// progress events. This is useful when progress monitoring is not needed
    /// or when the callback will be set later.
    ///
    /// # Return
    ///
    /// Returns a new `ProgressReporter` instance with no callback configured.
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates a new progress reporter with the specified callback.
    ///
    /// # Arguments
    ///
    /// * `callback` - The callback function to invoke for progress events.
    ///
    /// # Return
    ///
    /// Returns a new `ProgressReporter` instance configured with the provided callback.
    pub fn with_callback(callback: ProgressCallback<'a>) -> Self {
        Self {
            callback: Some(callback),
        }
    }

    /// Reports a progress event by invoking the configured callback.
    ///
    /// If no callback is configured, this method does nothing. The method is
    /// marked as inline to minimize overhead when progress reporting is disabled.
    ///
    /// # Arguments
    ///
    /// * `event` - The progress event to report.
    #[inline]
    pub fn report(&self, event: Progress) {
        if let Some(cb) = &self.callback {
            cb(event);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{Arc, Mutex};

    #[test]
    fn new_reporter_with_no_callback_does_not_panic_on_report() {
        let reporter = ProgressReporter::new();

        reporter.report(Progress::PhaseStart { name: "Test Phase" });
    }

    #[test]
    fn with_callback_reporter_invokes_the_provided_callback() {
        let received_event = Arc::new(Mutex::new(None));
        let received_event_clone = Arc::clone(&received_event);

        let callback = Box::new(move |event: Progress| {
            *received_event_clone.lock().unwrap() = Some(event);
        });

        let reporter = ProgressReporter::with_callback(callback);
        let event_to_send = Progress::TaskStart { total: 100 };

        reporter.report(event_to_send.clone());

        let received = received_event.lock().unwrap();
        assert!(received.is_some(), "Callback was not invoked");

        if let Some(Progress::TaskStart { total }) = received.as_ref() {
            assert_eq!(*total, 100);
        } else {
            panic!("Received event of the wrong type: {:?}", received);
        }
    }

    #[test]
    fn reporter_is_thread_safe_and_handles_concurrent_reports() {
        let received_count = Arc::new(Mutex::new(0));

        let callback = Box::new({
            let received_count_clone = Arc::clone(&received_count);
            move |_event: Progress| {
                *received_count_clone.lock().unwrap() += 1;
            }
        });

        let reporter = ProgressReporter::with_callback(callback);

        std::thread::scope(|s| {
            for _ in 0..10 {
                s.spawn(|| {
                    reporter.report(Progress::TaskIncrement { amount: 1 });
                });
            }
        });

        assert_eq!(*received_count.lock().unwrap(), 10);
    }

    #[test]
    fn progress_enum_variants_are_constructible() {
        let _ = Progress::PhaseStart { name: "phase" };
        let _ = Progress::PhaseFinish;
        let _ = Progress::TaskStart { total: 42 };
        let _ = Progress::TaskIncrement { amount: 1 };
        let _ = Progress::TaskFinish;
        let _ = Progress::StatusUpdate {
            text: "status".to_string(),
        };
        let _ = Progress::Message("msg".to_string());
    }
}

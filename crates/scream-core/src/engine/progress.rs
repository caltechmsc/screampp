#[derive(Debug, Clone, PartialEq)]
pub enum Progress {
    PhaseStart { name: &'static str },
    PhaseFinish,

    TaskStart { total: u64 },
    TaskIncrement { amount: u64 },
    TaskFinish,

    StatusUpdate { text: String },

    Message(String),
}

pub type ProgressCallback<'a> = Box<dyn Fn(Progress) + Send + Sync + 'a>;

#[derive(Default)]
pub struct ProgressReporter<'a> {
    callback: Option<ProgressCallback<'a>>,
}

impl<'a> ProgressReporter<'a> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_callback(callback: ProgressCallback<'a>) -> Self {
        Self {
            callback: Some(callback),
        }
    }

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

#[derive(Debug, Clone)]
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
    fn progress_reporter_new_has_no_callback() {
        let reporter = ProgressReporter::new();
        reporter.report(Progress::PhaseStart { name: "init" });
    }

    #[test]
    fn progress_reporter_with_callback_invokes_callback() {
        let called = Arc::new(Mutex::new(false));
        let called_clone = Arc::clone(&called);
        let reporter = ProgressReporter::with_callback(Box::new(move |event| {
            if let Progress::PhaseFinish = event {
                *called_clone.lock().unwrap() = true;
            }
        }));
        reporter.report(Progress::PhaseFinish);
        assert!(*called.lock().unwrap());
    }

    #[test]
    fn report_with_no_callback_does_not_panic() {
        let reporter = ProgressReporter::new();
        reporter.report(Progress::TaskFinish);
    }

    #[test]
    fn progress_callback_receives_correct_event() {
        let received = Arc::new(Mutex::new(None));
        let received_clone = Arc::clone(&received);
        let reporter = ProgressReporter::with_callback(Box::new(move |event| {
            if let Progress::Message(msg) = event {
                *received_clone.lock().unwrap() = Some(msg);
            }
        }));
        reporter.report(Progress::Message("hello".to_string()));
        assert_eq!(*received.lock().unwrap(), Some("hello".to_string()));
    }

    #[test]
    fn progress_enum_variants_are_constructible() {
        let _ = Progress::PhaseStart { name: "phase" };
        let _ = Progress::PhaseFinish;
        let _ = Progress::TaskStart { total_steps: 42 };
        let _ = Progress::TaskIncrement;
        let _ = Progress::TaskFinish;
        let _ = Progress::Message("msg".to_string());
    }
}

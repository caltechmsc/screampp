#[derive(Debug, Clone)]
pub enum Progress {
    PhaseStart { name: &'static str },
    PhaseFinish,

    TaskStart { total_steps: u64 },
    TaskIncrement,
    TaskFinish,

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

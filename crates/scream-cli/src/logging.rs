use crate::error::{CliError, Result};
use crate::utils::progress::CliProgressHandler;
use std::fs::File;
use std::path::PathBuf;
use tracing::Subscriber;
use tracing_subscriber::{
    filter::LevelFilter,
    fmt::{self},
    prelude::*,
    registry::LookupSpan,
};

struct IndicatifLayer {
    progress_handler: CliProgressHandler,
}

impl<S> tracing_subscriber::Layer<S> for IndicatifLayer
where
    S: Subscriber + for<'a> LookupSpan<'a>,
{
    fn on_event(
        &self,
        event: &tracing::Event<'_>,
        _ctx: tracing_subscriber::layer::Context<'_, S>,
    ) {
        let level = *event.metadata().level();
        let timestamp = chrono::Local::now().format("%Y-%m-%dT%H:%M:%S%.3fZ");

        let mut message = String::new();
        let mut visitor = StringVisitor(&mut message);
        event.record(&mut visitor);

        let log_line = format!("[{}] ({}): {}", timestamp, level, message);

        self.progress_handler.log(&log_line);
    }
}

struct StringVisitor<'a>(&'a mut String);
impl<'a> tracing::field::Visit for StringVisitor<'a> {
    fn record_debug(&mut self, field: &tracing::field::Field, value: &dyn std::fmt::Debug) {
        if field.name() == "message" {
            self.0.push_str(&format!("{:?}", value));
        }
    }
}

pub fn setup_logging(
    verbosity: u8,
    quiet: bool,
    log_file: &Option<PathBuf>,
    progress_handler: &CliProgressHandler,
) -> Result<()> {
    let console_level_filter = if quiet {
        LevelFilter::OFF
    } else {
        match verbosity {
            0 => LevelFilter::WARN,
            1 => LevelFilter::INFO,
            2 => LevelFilter::DEBUG,
            _ => LevelFilter::TRACE,
        }
    };

    let console_layer = IndicatifLayer {
        progress_handler: progress_handler.clone(),
    }
    .with_filter(console_level_filter);

    let file_layer = if let Some(path) = log_file {
        let file = File::create(path).map_err(CliError::Io)?;
        let layer = fmt::layer()
            .with_writer(file)
            .with_ansi(false)
            .with_thread_ids(true)
            .with_target(true);
        Some(layer.with_filter(LevelFilter::TRACE))
    } else {
        None
    };

    tracing_subscriber::registry()
        .with(console_layer)
        .with(file_layer)
        .init();

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use serial_test::serial;
    use std::sync::Once;
    use std::thread;
    use std::time::Duration;
    use tracing::{debug, error, info, trace, warn};

    static INIT: Once = Once::new();

    fn ensure_global_logger_is_set() {
        INIT.call_once(|| {
            setup_logging(3, false, &None).expect("Failed to set up global logger for tests");
        });
    }

    #[test]
    #[serial]
    fn initialization_and_macros_work() {
        ensure_global_logger_is_set();

        error!("This is an error");
        warn!("This is a warning");
        info!("This is info");
        debug!("This is debug");
        trace!("This is trace");
    }

    #[test]
    #[serial]
    fn file_logging_can_be_added_to_global_logger() {
        let temp_dir = tempfile::tempdir().unwrap();
        let log_path = temp_dir.path().join("test.log");

        let file = File::create(log_path.clone()).unwrap();
        let file_layer = fmt::layer()
            .with_writer(file)
            .with_ansi(false)
            .with_thread_ids(true);
        let subscriber = tracing_subscriber::registry().with(file_layer);

        tracing::subscriber::with_default(subscriber, || {
            let test_message = "Message for file-only test.";
            debug!("{}", test_message);
        });

        thread::sleep(Duration::from_millis(100));

        let content = std::fs::read_to_string(log_path).unwrap();
        assert!(content.contains("Message for file-only test."));
        assert!(content.contains("DEBUG"));
        assert!(content.contains("ThreadId"));
    }

    #[test]
    #[serial]
    fn invalid_log_file_path_propagates_error() {
        let invalid_path = PathBuf::from("/");

        if cfg!(unix) && invalid_path.is_dir() {
            let result = setup_logging(0, false, &Some(invalid_path));
            assert!(matches!(result, Err(CliError::Io(_))));
        }
    }
}

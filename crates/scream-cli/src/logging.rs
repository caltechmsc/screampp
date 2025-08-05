use crate::error::{CliError, Result};
use crate::ui::UiEvent;
use std::fs::File;
use std::path::PathBuf;
use tokio::sync::mpsc;
use tracing::Subscriber;
use tracing_subscriber::{
    filter::LevelFilter,
    fmt::{self},
    prelude::*,
    registry::LookupSpan,
};

struct ChannelLayer {
    sender: mpsc::Sender<UiEvent>,
}

impl<S> tracing_subscriber::Layer<S> for ChannelLayer
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

        let _ = self.sender.try_send(UiEvent::Log(log_line));
    }
}

struct StringVisitor<'a>(&'a mut String);
impl<'a> tracing::field::Visit for StringVisitor<'a> {
    fn record_debug(&mut self, field: &tracing::field::Field, value: &dyn std::fmt::Debug) {
        if field.name() == "message" {
            let mut formatted = format!("{:?}", value);
            if formatted.starts_with('"') && formatted.ends_with('"') {
                formatted = formatted[1..formatted.len() - 1].to_string();
            }
            self.0.push_str(&formatted);
        }
    }
}

pub fn setup_logging(
    verbosity: u8,
    quiet: bool,
    log_file: &Option<PathBuf>,
    ui_sender: mpsc::Sender<UiEvent>,
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

    let console_layer = ChannelLayer { sender: ui_sender }.with_filter(console_level_filter);

    let file_layer = if let Some(path) = log_file {
        let file = File::create(path).map_err(CliError::Io)?;
        let layer = fmt::layer()
            .with_writer(file)
            .with_ansi(false)
            .with_thread_ids(true)
            .with_target(true)
            .compact();
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

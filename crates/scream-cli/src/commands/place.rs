use crate::cli::PlaceArgs;
use crate::config::PartialPlacementConfig;
use crate::data::DataManager;
use crate::error::{CliError, Result};
use crate::utils::progress::CliProgressHandler;
use screampp::{
    core::io::{bgf::BgfFile, traits::MolecularFile},
    engine::progress::ProgressReporter,
    workflows,
};
use std::path::{Path, PathBuf};
use tracing::{info, warn};

pub async fn run(args: PlaceArgs) -> Result<()> {
    info!("Initializing data manager...");
    let data_manager = DataManager::new()?;

    let partial_config = PartialPlacementConfig::from_file(&args.config)?;
    info!("Merging configuration from file and CLI arguments...");
    let final_config = partial_config.merge_with_cli(&args, &data_manager)?;

    info!("Loading input structure from {:?}", &args.input);
    let (system, metadata) =
        BgfFile::read_from_path(&args.input).map_err(|e| CliError::FileParsing {
            path: args.input.clone(),
            source: e.into(),
        })?;

    let progress_handler = CliProgressHandler::new();
    let reporter = ProgressReporter::with_callback(progress_handler.get_callback());

    println!("Starting side-chain placement...");
    info!("Invoking the core placement workflow...");

    let solutions =
        tokio::task::block_in_place(|| workflows::place::run(&system, &final_config, &reporter))?;

    info!(
        "Workflow finished, received {} solution(s).",
        solutions.len()
    );

    if solutions.is_empty() {
        warn!("Workflow completed but found no valid solutions.");
        println!("Warning: SCREAM finished but found no valid solutions.");
    } else {
        println!(
            "Workflow complete. Writing {} solution(s)...",
            solutions.len()
        );

        for (i, solution) in solutions.iter().enumerate() {
            let output_path = generate_output_path(&args.output, i + 1, solutions.len());
            info!(
                "Writing solution {} (Energy: {:.4}) to {:?}",
                i + 1,
                solution.energy,
                &output_path
            );

            BgfFile::write_to(
                &solution.state.system,
                &metadata,
                &mut std::fs::File::create(&output_path)?,
            )
            .map_err(|e| CliError::FileParsing {
                path: output_path.clone(),
                source: e.into(),
            })?;

            if i == 0 {
                println!(
                    "✓ Best solution (Energy: {:.4} kcal/mol) written to: {}",
                    solution.energy,
                    output_path.display()
                );
            } else {
                println!(
                    "  Solution {} (Energy: {:.4} kcal/mol) written to: {}",
                    i + 1,
                    solution.energy,
                    output_path.display()
                );
            }
        }
    }

    Ok(())
}

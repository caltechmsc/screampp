use crate::cli::PlaceArgs;
use crate::config::PartialPlacementConfig;
use crate::data::DataManager;
use crate::error::{CliError, Result};
use crate::utils::progress::CliProgressHandler;
use screampp::{
    core::io::{bgf::BgfFile, traits::MolecularFile},
    engine::{progress::ProgressReporter, state::Solution},
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
    let (mut system, metadata) =
        BgfFile::read_from_path(&args.input).map_err(|e| CliError::FileParsing {
            path: args.input.clone(),
            source: e.into(),
        })?;

    let progress_handler = CliProgressHandler::new();
    let reporter = ProgressReporter::with_callback(progress_handler.get_callback());

    println!("Starting side-chain placement...");
    info!("Invoking the core placement workflow...");

    let solutions = tokio::task::block_in_place(|| {
        workflows::place::run(&mut system, &final_config, &reporter)
    })?;

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
            let output_path =
                generate_output_path(&args.output, solution, i + 1, solutions.len());
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

fn generate_output_path(
    base_path_template: &Path,
    solution: &Solution,
    index: usize,
    total: usize,
) -> PathBuf {
    if total <= 1 {
        return base_path_template.to_path_buf();
    }

    let template_str = match base_path_template.to_str() {
        Some(s) => s,
        None => {
            return generate_indexed_path(base_path_template, index);
        }
    };

    let contains_placeholder = ["{n}", "{i}", "{N}", "{total}", "{energy}"]
        .iter()
        .any(|p| template_str.contains(p));

    if contains_placeholder {
        let path_str = template_str
            .replace("{n}", &index.to_string())
            .replace("{i}", &index.to_string())
            .replace("{N}", &total.to_string())
            .replace("{total}", &total.to_string())
            .replace("{energy}", &format!("{:.2}", solution.energy));
        PathBuf::from(path_str)
    } else {
        generate_indexed_path(base_path_template, index)
    }
}

fn generate_indexed_path(base_path: &Path, index: usize) -> PathBuf {
    let stem = base_path
        .file_stem()
        .unwrap_or_default()
        .to_str()
        .unwrap_or("");
    let extension = base_path
        .extension()
        .unwrap_or_default()
        .to_str()
        .unwrap_or("");

    let new_filename = if extension.is_empty() {
        format!("{}-best-{}", stem, index)
    } else {
        format!("{}-best-{}.{}", stem, index, extension)
    };

    base_path.with_file_name(new_filename)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_generate_output_path_single_solution() {
        let path = PathBuf::from("output.bgf");
        let new_path = generate_output_path(&path, 1, 1);
        assert_eq!(new_path, path);
    }

    #[test]
    fn test_generate_output_path_multiple_solutions() {
        let path = PathBuf::from("/tmp/results/output.bgf");
        let new_path_1 = generate_output_path(&path, 1, 3);
        let new_path_2 = generate_output_path(&path, 2, 3);

        assert_eq!(new_path_1, PathBuf::from("/tmp/results/output-best-1.bgf"));
        assert_eq!(new_path_2, PathBuf::from("/tmp/results/output-best-2.bgf"));
    }

    #[test]
    fn test_generate_output_path_no_extension() {
        let path = PathBuf::from("my_output");
        let new_path = generate_output_path(&path, 5, 10);
        assert_eq!(new_path, PathBuf::from("my_output-best-5"));
    }

    #[test]
    fn test_generate_output_path_with_dots_in_stem() {
        let path = PathBuf::from("protein.v1.bgf");
        let new_path = generate_output_path(&path, 2, 2);
        assert_eq!(new_path, PathBuf::from("protein.v1-best-2.bgf"));
    }
}

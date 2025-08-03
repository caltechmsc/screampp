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
            let output_path = generate_output_path(&args.output, solution, i + 1, solutions.len());
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
    use screampp::core::models::system::MolecularSystem;
    use screampp::engine::state::{Solution, SolutionState};
    use std::collections::HashMap;
    use std::path::PathBuf;

    fn mock_solution(energy: f64) -> Solution {
        Solution {
            energy,
            state: SolutionState {
                system: MolecularSystem::new(),
                rotamers: HashMap::new(),
            },
        }
    }

    const PLACEHOLDERS: &[&str] = &["{n}", "{i}", "{N}", "{total}", "{energy}"];

    #[test]
    fn generate_indexed_path_should_add_suffix_to_filename_with_extension() {
        let base_path = PathBuf::from("/path/to/result.bgf");
        let result = generate_indexed_path(&base_path, 3);
        assert_eq!(result, PathBuf::from("/path/to/result-best-3.bgf"));
    }

    #[test]
    fn generate_indexed_path_should_add_suffix_to_filename_without_extension() {
        let base_path = PathBuf::from("output");
        let result = generate_indexed_path(&base_path, 5);
        assert_eq!(result, PathBuf::from("output-best-5"));
    }

    #[test]
    fn generate_indexed_path_should_handle_multiple_dots_in_filename() {
        let base_path = PathBuf::from("my.protein.v1.pdb");
        let result = generate_indexed_path(&base_path, 1);
        assert_eq!(result, PathBuf::from("my.protein.v1-best-1.pdb"));
    }

    #[test]
    fn generate_output_path_should_return_original_path_when_total_is_one() {
        let path_template = PathBuf::from("single_result.bgf");
        let solution = mock_solution(-10.0);

        let result = generate_output_path(&path_template, &solution, 1, 1);

        assert_eq!(result, path_template);
    }

    #[test]
    fn generate_output_path_should_use_indexed_fallback_when_multiple_solutions_and_no_placeholder()
    {
        let path_template = PathBuf::from("fallback.bgf");
        let solution = mock_solution(-20.0);

        let result = generate_output_path(&path_template, &solution, 2, 5);

        assert_eq!(result, PathBuf::from("fallback-best-2.bgf"));
    }

    #[test]
    fn generate_output_path_should_replace_index_placeholders_n_and_i() {
        let solution = mock_solution(-30.0);
        let path_with_n = PathBuf::from("solution-{n}.bgf");
        let path_with_i = PathBuf::from("solution-{i}.bgf");

        let result_n = generate_output_path(&path_with_n, &solution, 4, 10);
        let result_i = generate_output_path(&path_with_i, &solution, 4, 10);

        assert_eq!(result_n, PathBuf::from("solution-4.bgf"));
        assert_eq!(result_i, PathBuf::from("solution-4.bgf"));
    }

    #[test]
    fn generate_output_path_should_replace_total_placeholders_n_and_total() {
        let solution = mock_solution(-40.0);
        let path_with_cap_n = PathBuf::from("solution_1_of_{N}.bgf");
        let path_with_total = PathBuf::from("solution_1_of_{total}.bgf");

        let result_cap_n = generate_output_path(&path_with_cap_n, &solution, 1, 7);
        let result_total = generate_output_path(&path_with_total, &solution, 1, 7);

        assert_eq!(result_cap_n, PathBuf::from("solution_1_of_7.bgf"));
        assert_eq!(result_total, PathBuf::from("solution_1_of_7.bgf"));
    }

    #[test]
    fn generate_output_path_should_replace_energy_placeholder_with_two_decimals() {
        let solution = mock_solution(-55.5555);
        let path_template = PathBuf::from("energy_is_{energy}.bgf");

        let result = generate_output_path(&path_template, &solution, 1, 2);

        assert_eq!(result, PathBuf::from("energy_is_-55.56.bgf"));
    }

    #[test]
    fn generate_output_path_should_replace_all_placeholders_at_once() {
        let solution = mock_solution(-99.9);
        let path_template = PathBuf::from("run_{i}_of_{N}_E={energy}.pdb");

        let result = generate_output_path(&path_template, &solution, 9, 10);

        assert_eq!(result, PathBuf::from("run_9_of_10_E=-99.90.pdb"));
    }
}

use super::{config::InputFiles, error::ScreamError};
use crate::core::forcefield::params::Forcefield;
use crate::core::io::bgf::BgfFile;
use crate::core::io::traits::MolecularFile;
use crate::core::models::system::MolecularSystem;
use crate::core::rotamers::placement::{PlacementRegistry, load_placement_registry};
use tracing::info;

#[derive(Debug)]
pub struct ScreamContext {
    pub system: MolecularSystem,
    pub forcefield: Forcefield,
    pub placement_registry: PlacementRegistry,
}

impl ScreamContext {
    pub fn new(paths: &InputFiles) -> Result<Self, ScreamError> {
        // --- 1. Load Molecular Structure ---
        info!(path = ?paths.structure, "Loading molecular system...");
        let (system, _) =
            BgfFile::read_from_path(&paths.structure).map_err(ScreamError::BfgParse)?;

        // --- 2. Load All Forcefield-Related Parameters ---
        info!("Loading forcefield parameters...");
        let forcefield = Forcefield::load(
            &paths.non_bonded_params,
            &paths.delta_params,
            &paths.charge_params,
            &paths.topology_params,
        )
        .map_err(|e| ScreamError::Load(e.to_string()))?;

        // --- 3. Load Placement Information ---
        info!(path = ?paths.placement_params, "Loading placement registry...");
        let placement_registry = load_placement_registry(&paths.placement_params)
            .map_err(|e| ScreamError::Load(e.to_string()))?;

        info!("ScreamContext created successfully.");
        Ok(Self {
            system,
            forcefield,
            placement_registry,
        })
    }
}

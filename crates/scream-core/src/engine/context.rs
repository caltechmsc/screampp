use super::config::InputPaths;
use super::error::ScreamError;
use super::progress::ProgressUpdate;
use crate::core::forcefield::parameterization::Parameterizer;
use crate::core::forcefield::params::Forcefield;
use crate::core::forcefield::placement::{PlacementRegistry, load_placement_registry};
use crate::core::io::bgf::BgfFile;
use crate::core::io::traits::MolecularFile;
use crate::core::models::system::MolecularSystem;
use crate::core::rotamers::library::RotamerLibrary;
use tracing::{debug, info, instrument};

#[derive(Debug)]
pub struct ScreamContext {
    pub system: MolecularSystem,
    pub forcefield: Forcefield,
    pub rotamer_library: RotamerLibrary,
    pub placement_registry: PlacementRegistry,
}

impl ScreamContext {
    #[instrument(name = "context_creation", skip_all)]
    pub fn new<F>(
        paths: &InputPaths,
        s_factor: f64,
        mut on_progress: F,
    ) -> Result<Self, ScreamError>
    where
        F: FnMut(ProgressUpdate),
    {
        info!("Initializing SCREAM context...");
        on_progress(ProgressUpdate::New(
            "Preparing context...".to_string(),
            Some(6),
        ));

        // --- Step 1: Load Forcefield Parameters ---
        on_progress(ProgressUpdate::Status("Loading forcefield parameters..."));
        debug!(
            "Forcefield path: {:?}, Delta: {:?}, Charge: {:?}, Topology: {:?}",
            paths.forcefield_file, paths.delta_file, paths.charge_file, paths.topology_file
        );
        let forcefield = Forcefield::load(
            &paths.forcefield_file,
            &paths.delta_file,
            &paths.charge_file,
            &paths.topology_file,
        )
        .map_err(|e| ScreamError::Forcefield(e.to_string()))?;
        on_progress(ProgressUpdate::Increment(1));

        // --- Step 2: Load Placement Registry ---
        on_progress(ProgressUpdate::Status("Loading placement registry..."));
        debug!("Placement registry path: {:?}", paths.placement_file);
        let placement_registry = load_placement_registry(&paths.placement_file)?;
        on_progress(ProgressUpdate::Increment(1));

        // --- Step 3: Load Molecular System ---
        on_progress(ProgressUpdate::Status("Loading molecular system..."));
        debug!("Structure file path: {:?}", paths.structure_file);
        let (mut system, _) = BgfFile::read_from_path(&paths.structure_file)?;
        on_progress(ProgressUpdate::Increment(1));

        // --- Step 4: Parameterize the System ---
        on_progress(ProgressUpdate::Status("Parameterizing system..."));
        let parameterizer = Parameterizer::new(forcefield.clone(), s_factor);
        parameterizer.parameterize_system(&mut system)?;
        on_progress(ProgressUpdate::Increment(1));

        // --- Step 5: Load and Parameterize Rotamer Library ---
        on_progress(ProgressUpdate::Status("Loading rotamer library..."));
        debug!("Rotamer library path: {:?}", paths.rotamer_library_file);
        let rotamer_library = RotamerLibrary::load(
            &paths.rotamer_library_file,
            &placement_registry,
            &forcefield,
            s_factor,
        )?;
        on_progress(ProgressUpdate::Increment(1));

        // --- Step 6: Finalization ---
        on_progress(ProgressUpdate::Status("Context created successfully."));
        info!("SCREAM context successfully initialized.");
        on_progress(ProgressUpdate::Increment(1));
        on_progress(ProgressUpdate::Done);

        Ok(Self {
            system,
            forcefield,
            rotamer_library,
            placement_registry,
        })
    }
}

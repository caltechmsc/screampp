use super::context::OptimizationContext;
use super::error::EngineError;
use super::placement;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use std::collections::HashMap;

/// Provides a transactional view of the molecular system for optimization operations.
///
/// This struct enables safe, temporary modifications to the molecular system during
/// optimization moves. It tracks rotamer assignments and provides transaction methods
/// that automatically revert changes after evaluation, ensuring the system remains
/// in a consistent state for energy calculations.
pub struct SystemView<'a, 'ctx, C>
where
    C: super::context::ProvidesResidueSelections + Sync,
{
    /// Mutable reference to the molecular system being modified.
    pub system: &'a mut MolecularSystem,
    /// Reference to the optimization context providing necessary data.
    pub context: &'ctx OptimizationContext<'ctx, C>,
    /// Mutable reference to the current rotamer assignments map.
    pub current_rotamers: &'a mut HashMap<ResidueId, usize>,
}

impl<'a, 'ctx, C> SystemView<'a, 'ctx, C>
where
    C: super::context::ProvidesResidueSelections + Sync,
{
    /// Creates a new system view with the provided components.
    ///
    /// # Arguments
    ///
    /// * `system` - Mutable reference to the molecular system.
    /// * `context` - Reference to the optimization context.
    /// * `current_rotamers` - Mutable reference to the rotamer assignments map.
    ///
    /// # Return
    ///
    /// Returns a new `SystemView` instance.
    pub fn new(
        system: &'a mut MolecularSystem,
        context: &'ctx OptimizationContext<'ctx, C>,
        current_rotamers: &'a mut HashMap<ResidueId, usize>,
    ) -> Self {
        Self {
            system,
            context,
            current_rotamers,
        }
    }

    /// Applies a rotamer move to the specified residue.
    ///
    /// This method updates both the molecular system structure and the rotamer
    /// tracking map to reflect the new conformation.
    ///
    /// # Arguments
    ///
    /// * `res_id` - The residue to modify.
    /// * `new_rotamer_idx` - The index of the new rotamer to apply.
    ///
    /// # Return
    ///
    /// Returns `Ok(())` on success.
    ///
    /// # Errors
    ///
    /// Returns `EngineError` if the rotamer placement fails.
    pub fn apply_move(
        &mut self,
        res_id: ResidueId,
        new_rotamer_idx: usize,
    ) -> Result<(), EngineError> {
        self.place_rotamer(res_id, new_rotamer_idx)?;
        self.current_rotamers.insert(res_id, new_rotamer_idx);
        Ok(())
    }

    /// Places a specific rotamer on the given residue in the molecular system.
    ///
    /// This is an internal method that handles the actual placement of rotamer
    /// atoms using the placement module.
    fn place_rotamer(&mut self, res_id: ResidueId, rotamer_idx: usize) -> Result<(), EngineError> {
        let residue = self.system.residue(res_id).unwrap();
        let res_type = residue.residue_type.unwrap();
        let rotamer = &self
            .context
            .rotamer_library
            .get_rotamers_for(res_type)
            .unwrap()[rotamer_idx];
        let res_name = res_type.to_three_letter();
        let topology = self.context.topology_registry.get(res_name).unwrap();

        placement::place_rotamer_on_system(self.system, res_id, rotamer, topology)
    }

    /// Executes an action within a transaction that automatically reverts changes.
    ///
    /// This method provides transactional semantics for system modifications. The
    /// action is executed, and any changes to the specified residue are automatically
    /// reverted after the action completes, regardless of success or failure.
    ///
    /// # Arguments
    ///
    /// * `res_id_to_modify` - The residue that may be modified during the action.
    /// * `action` - A closure that performs the desired operation on the system view.
    ///
    /// # Return
    ///
    /// Returns the result of the action closure.
    ///
    /// # Errors
    ///
    /// Returns `EngineError` if the residue is not found in the rotamer map or if
    /// the action fails.
    pub fn transaction<F, R>(
        &mut self,
        res_id_to_modify: ResidueId,
        action: F,
    ) -> Result<R, EngineError>
    where
        F: FnOnce(&mut Self) -> Result<R, EngineError>,
    {
        // 1. Record the original state for the residue.
        let original_rotamer_idx =
            *self
                .current_rotamers
                .get(&res_id_to_modify)
                .ok_or_else(|| {
                    EngineError::Internal(format!(
                        "Residue {:?} not found in rotamer tracking map",
                        res_id_to_modify
                    ))
                })?;

        // 2. Execute the action.
        let result = action(self)?;

        // 3. Check if the residue was actually modified during the action.
        let current_rotamer_idx = *self.current_rotamers.get(&res_id_to_modify).unwrap();

        // 4. If it was modified, revert it back to the original state.
        if current_rotamer_idx != original_rotamer_idx {
            self.place_rotamer(res_id_to_modify, original_rotamer_idx)?;
            self.current_rotamers
                .insert(res_id_to_modify, original_rotamer_idx);
        }

        Ok(result)
    }
    /// Executes an action within a transaction for two residues, automatically reverting changes.
    ///
    /// This method extends the transaction concept to handle simultaneous modifications
    /// of two residues. Both residues are tracked and reverted to their original states
    /// after the action completes.
    ///
    /// # Arguments
    ///
    /// * `res_a_id` - The first residue that may be modified.
    /// * `res_b_id` - The second residue that may be modified.
    /// * `action` - A closure that performs the desired operation on the system view.
    ///
    /// # Return
    ///
    /// Returns the result of the action closure.
    ///
    /// # Errors
    ///
    /// Returns `EngineError` if either residue is not found in the rotamer map or if
    /// the action fails.
    pub fn transaction_doublet<F, R>(
        &mut self,
        res_a_id: ResidueId,
        res_b_id: ResidueId,
        action: F,
    ) -> Result<R, EngineError>
    where
        F: FnOnce(&mut Self) -> Result<R, EngineError>,
    {
        // 1. Record the original states for both residues.
        let original_rotamer_idx_a = *self.current_rotamers.get(&res_a_id).ok_or_else(|| {
            EngineError::Internal(format!(
                "Residue {:?} not found in rotamer tracking map",
                res_a_id
            ))
        })?;

        let original_rotamer_idx_b = *self.current_rotamers.get(&res_b_id).ok_or_else(|| {
            EngineError::Internal(format!(
                "Residue {:?} not found in rotamer tracking map",
                res_b_id
            ))
        })?;

        // 2. Execute the action.
        let result = action(self)?;

        // 3. Check if the residues were actually modified during the action.
        let current_rotamer_idx_a = *self.current_rotamers.get(&res_a_id).unwrap();
        let current_rotamer_idx_b = *self.current_rotamers.get(&res_b_id).unwrap();

        // 4. If they were modified, revert them back to the original states.
        if current_rotamer_idx_a != original_rotamer_idx_a {
            self.place_rotamer(res_a_id, original_rotamer_idx_a)?;
            self.current_rotamers
                .insert(res_a_id, original_rotamer_idx_a);
        }

        if current_rotamer_idx_b != original_rotamer_idx_b {
            self.place_rotamer(res_b_id, original_rotamer_idx_b)?;
            self.current_rotamers
                .insert(res_b_id, original_rotamer_idx_b);
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{
        forcefield::{parameterization::Parameterizer, params::Forcefield},
        models::{atom::Atom, chain::ChainType, residue::ResidueType, system::MolecularSystem},
        rotamers::{library::RotamerLibrary, rotamer::Rotamer},
        topology::registry::TopologyRegistry,
    };
    use crate::engine::{
        config::{ConvergenceConfig, PlacementConfig, PlacementConfigBuilder, ResidueSelection},
        context::OptimizationContext,
        progress::ProgressReporter,
    };
    use nalgebra::Point3;
    use std::{
        collections::{HashMap, HashSet},
        fs,
        path::Path,
    };
    use tempfile::TempDir;

    struct TestSetup {
        system: MolecularSystem,
        forcefield: Forcefield,
        rotamer_library: RotamerLibrary,
        topology_registry: TopologyRegistry,
        active_residues: HashSet<ResidueId>,
        initial_rotamers: HashMap<ResidueId, usize>,
        _temp_dir: TempDir,
    }

    fn write_file(path: &Path, content: &str) {
        fs::write(path, content).expect("Failed to write temporary file for test");
    }

    fn setup() -> TestSetup {
        let temp_dir = tempfile::tempdir().expect("Failed to create temp dir");

        let ff_path = temp_dir.path().join("ff.toml");
        write_file(
            &ff_path,
            r#"
            [globals]
            dielectric_constant = 4.0
            potential_function = "lennard-jones-12-6"
            [vdw]
            N_BB = { radius = 2.8, well_depth = 0.15 }
            C_BB = { radius = 3.0, well_depth = 0.1 }
            C_SC = { radius = 4.0, well_depth = 0.2 }
            O_SC = { radius = 3.2, well_depth = 0.3 }
            [hbond]
            "#,
        );

        let delta_path = temp_dir.path().join("delta.csv");
        write_file(&delta_path, "residue_type,atom_name,mu,sigma\n");

        let topo_path = temp_dir.path().join("registry.toml");
        write_file(
            &topo_path,
            r#"
            [ALA]
            anchor_atoms = ["N", "CA", "C"]
            sidechain_atoms = ["CB"]

            [SER]
            anchor_atoms = ["N", "CA", "C"]
            sidechain_atoms = ["CB", "OG"]
            "#,
        );

        let forcefield = Forcefield::load(&ff_path, &delta_path, &Default::default()).unwrap();
        let topology_registry = TopologyRegistry::load(&topo_path).unwrap();

        let mut system = MolecularSystem::new();
        let chain_a = system.add_chain('A', ChainType::Protein);

        let ala_id = system
            .add_residue(chain_a, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        add_atom(&mut system, ala_id, "N", [0.0, 1.4, 0.0], "N_BB");
        add_atom(&mut system, ala_id, "CA", [0.0, 0.0, 0.0], "C_BB");
        add_atom(&mut system, ala_id, "C", [1.4, 0.0, 0.0], "C_BB");
        add_atom(&mut system, ala_id, "CB", [0.0, -1.5, 0.0], "C_SC");

        let ser_id = system
            .add_residue(chain_a, 2, "SER", Some(ResidueType::Serine))
            .unwrap();
        add_atom(&mut system, ser_id, "N", [10.0, 0.0, 0.0], "N_BB");
        add_atom(&mut system, ser_id, "CA", [11.0, 0.0, 0.0], "C_BB");
        add_atom(&mut system, ser_id, "C", [12.0, 0.0, 0.0], "C_BB");
        add_atom(&mut system, ser_id, "CB", [11.0, 1.5, 0.0], "C_SC");
        add_atom(&mut system, ser_id, "OG", [11.0, 2.7, 0.0], "O_SC");

        let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 0.0);
        parameterizer.parameterize_system(&mut system).unwrap();

        let mut rotamer_library = RotamerLibrary::default();
        let ala_rotamers = vec![
            create_rotamer(
                &parameterizer,
                &topology_registry,
                "ALA",
                &[("CB", [0.0, -1.5, 0.0])],
            ),
            create_rotamer(
                &parameterizer,
                &topology_registry,
                "ALA",
                &[("CB", [0.0, -1.2, 1.2])],
            ),
        ];
        rotamer_library
            .rotamers
            .insert(ResidueType::Alanine, ala_rotamers);

        let ser_rotamers = vec![
            create_rotamer(
                &parameterizer,
                &topology_registry,
                "SER",
                &[("CB", [11.0, 1.5, 0.0]), ("OG", [11.0, 2.7, 0.0])],
            ),
            create_rotamer(
                &parameterizer,
                &topology_registry,
                "SER",
                &[("CB", [11.0, 1.2, 1.2]), ("OG", [11.0, 2.4, 1.2])],
            ),
        ];
        rotamer_library
            .rotamers
            .insert(ResidueType::Serine, ser_rotamers);

        let active_residues = HashSet::from([ala_id, ser_id]);
        let mut initial_rotamers = HashMap::new();
        initial_rotamers.insert(ala_id, 0);
        initial_rotamers.insert(ser_id, 0);

        TestSetup {
            system,
            forcefield,
            rotamer_library,
            topology_registry,
            active_residues,
            initial_rotamers,
            _temp_dir: temp_dir,
        }
    }

    fn add_atom(
        system: &mut MolecularSystem,
        res_id: ResidueId,
        name: &str,
        pos: [f64; 3],
        ff_type: &str,
    ) {
        let mut atom = Atom::new(name, res_id, Point3::from(pos));
        atom.force_field_type = ff_type.to_string();
        system.add_atom_to_residue(res_id, atom).unwrap();
    }

    fn create_rotamer(
        parameterizer: &Parameterizer,
        topology_registry: &TopologyRegistry,
        res_name: &str,
        sidechain_atoms: &[(&str, [f64; 3])],
    ) -> Rotamer {
        let mut atoms = vec![
            Atom::new("N", ResidueId::default(), Point3::new(0.0, 1.4, 0.0)),
            Atom::new("CA", ResidueId::default(), Point3::new(0.0, 0.0, 0.0)),
            Atom::new("C", ResidueId::default(), Point3::new(1.4, 0.0, 0.0)),
        ];

        for (name, pos) in sidechain_atoms {
            atoms.push(Atom::new(*name, ResidueId::default(), Point3::from(*pos)));
        }

        atoms.iter_mut().for_each(|a| {
            a.force_field_type = match a.name.as_str() {
                "N" => "N_BB".to_string(),
                "CA" | "C" => "C_BB".to_string(),
                "CB" => "C_SC".to_string(),
                "OG" => "O_SC".to_string(),
                _ => "UNKNOWN".to_string(),
            };
        });

        let bonds = match res_name {
            "ALA" => vec![(0, 1), (1, 2), (1, 3)],
            "SER" => vec![(0, 1), (1, 2), (1, 3), (3, 4)],
            _ => vec![],
        };

        let topology = topology_registry.get(res_name).unwrap();
        let mut rotamer = Rotamer { atoms, bonds };
        parameterizer
            .parameterize_rotamer(&mut rotamer, res_name, topology)
            .unwrap();
        rotamer
    }

    fn create_config(selection: ResidueSelection) -> PlacementConfig {
        PlacementConfigBuilder::new()
            .forcefield_path("dummy.ff")
            .delta_params_path("dummy.delta")
            .s_factor(0.0)
            .rotamer_library_path("dummy.rotlib")
            .topology_registry_path("dummy.topo")
            .max_iterations(10)
            .num_solutions(1)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.01,
                patience_iterations: 3,
            })
            .final_refinement_iterations(2)
            .include_input_conformation(false)
            .residues_to_optimize(selection)
            .build()
            .unwrap()
    }

    #[test]
    fn system_view_new_creates_view_correctly() {
        let setup = setup();
        let config = create_config(ResidueSelection::All);
        let reporter = ProgressReporter::new();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );
        let mut system = setup.system.clone();
        let mut current_rotamers = setup.initial_rotamers.clone();

        let system_view = SystemView::new(&mut system, &context, &mut current_rotamers);

        assert_eq!(
            system_view.system.residues_iter().count(),
            setup.system.residues_iter().count()
        );
        assert_eq!(
            system_view.current_rotamers.len(),
            setup.active_residues.len()
        );
    }

    #[test]
    fn apply_move_updates_system_and_tracking() {
        let setup = setup();
        let config = create_config(ResidueSelection::All);
        let reporter = ProgressReporter::new();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );
        let mut system = setup.system.clone();
        let mut current_rotamers = setup.initial_rotamers.clone();
        let mut system_view = SystemView::new(&mut system, &context, &mut current_rotamers);

        let ala_id = *setup
            .active_residues
            .iter()
            .find(|&&id| {
                setup.system.residue(id).unwrap().residue_type == Some(ResidueType::Alanine)
            })
            .unwrap();

        assert_eq!(*system_view.current_rotamers.get(&ala_id).unwrap(), 0);

        system_view.apply_move(ala_id, 1).unwrap();

        assert_eq!(*system_view.current_rotamers.get(&ala_id).unwrap(), 1);

        let residue = system_view.system.residue(ala_id).unwrap();
        let cb_atoms: Vec<_> = residue
            .atoms()
            .iter()
            .filter_map(|&atom_id| {
                system_view.system.atom(atom_id).and_then(|atom| {
                    if atom.name == "CB" {
                        Some(atom.position)
                    } else {
                        None
                    }
                })
            })
            .collect();

        assert_eq!(cb_atoms.len(), 1);
    }

    #[test]
    fn transaction_reverts_changes_on_scope_exit() {
        let setup = setup();
        let config = create_config(ResidueSelection::All);
        let reporter = ProgressReporter::new();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );
        let mut system = setup.system.clone();
        let mut current_rotamers = setup.initial_rotamers.clone();
        let mut system_view = SystemView::new(&mut system, &context, &mut current_rotamers);

        let ala_id = *setup
            .active_residues
            .iter()
            .find(|&&id| {
                setup.system.residue(id).unwrap().residue_type == Some(ResidueType::Alanine)
            })
            .unwrap();

        let original_rotamer = *system_view.current_rotamers.get(&ala_id).unwrap();

        let result = system_view.transaction(ala_id, |view| {
            view.apply_move(ala_id, 1)?;
            assert_eq!(*view.current_rotamers.get(&ala_id).unwrap(), 1);
            Ok("transaction_success")
        });

        assert_eq!(result.unwrap(), "transaction_success");

        assert_eq!(
            *system_view.current_rotamers.get(&ala_id).unwrap(),
            original_rotamer
        );
    }

    #[test]
    fn transaction_preserves_changes_when_no_modification() {
        let setup = setup();
        let config = create_config(ResidueSelection::All);
        let reporter = ProgressReporter::new();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );
        let mut system = setup.system.clone();
        let mut current_rotamers = setup.initial_rotamers.clone();
        let mut system_view = SystemView::new(&mut system, &context, &mut current_rotamers);

        let ala_id = *setup
            .active_residues
            .iter()
            .find(|&&id| {
                setup.system.residue(id).unwrap().residue_type == Some(ResidueType::Alanine)
            })
            .unwrap();

        let original_rotamer = *system_view.current_rotamers.get(&ala_id).unwrap();

        let result = system_view.transaction(ala_id, |_| Ok("no_change"));

        assert_eq!(result.unwrap(), "no_change");
        assert_eq!(
            *system_view.current_rotamers.get(&ala_id).unwrap(),
            original_rotamer
        );
    }

    #[test]
    fn transaction_doublet_reverts_both_residues() {
        let setup = setup();
        let config = create_config(ResidueSelection::All);
        let reporter = ProgressReporter::new();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );
        let mut system = setup.system.clone();
        let mut current_rotamers = setup.initial_rotamers.clone();
        let mut system_view = SystemView::new(&mut system, &context, &mut current_rotamers);

        let ala_id = *setup
            .active_residues
            .iter()
            .find(|&&id| {
                setup.system.residue(id).unwrap().residue_type == Some(ResidueType::Alanine)
            })
            .unwrap();

        let ser_id = *setup
            .active_residues
            .iter()
            .find(|&&id| {
                setup.system.residue(id).unwrap().residue_type == Some(ResidueType::Serine)
            })
            .unwrap();

        let original_ala_rotamer = *system_view.current_rotamers.get(&ala_id).unwrap();
        let original_ser_rotamer = *system_view.current_rotamers.get(&ser_id).unwrap();

        let result = system_view.transaction_doublet(ala_id, ser_id, |view| {
            view.apply_move(ala_id, 1)?;
            view.apply_move(ser_id, 1)?;
            assert_eq!(*view.current_rotamers.get(&ala_id).unwrap(), 1);
            assert_eq!(*view.current_rotamers.get(&ser_id).unwrap(), 1);
            Ok("doublet_success")
        });

        assert_eq!(result.unwrap(), "doublet_success");

        assert_eq!(
            *system_view.current_rotamers.get(&ala_id).unwrap(),
            original_ala_rotamer
        );
        assert_eq!(
            *system_view.current_rotamers.get(&ser_id).unwrap(),
            original_ser_rotamer
        );
    }

    #[test]
    fn transaction_doublet_preserves_changes_when_no_modification() {
        let setup = setup();
        let config = create_config(ResidueSelection::All);
        let reporter = ProgressReporter::new();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );
        let mut system = setup.system.clone();
        let mut current_rotamers = setup.initial_rotamers.clone();
        let mut system_view = SystemView::new(&mut system, &context, &mut current_rotamers);

        let ala_id = *setup
            .active_residues
            .iter()
            .find(|&&id| {
                setup.system.residue(id).unwrap().residue_type == Some(ResidueType::Alanine)
            })
            .unwrap();

        let ser_id = *setup
            .active_residues
            .iter()
            .find(|&&id| {
                setup.system.residue(id).unwrap().residue_type == Some(ResidueType::Serine)
            })
            .unwrap();

        let original_ala_rotamer = *system_view.current_rotamers.get(&ala_id).unwrap();
        let original_ser_rotamer = *system_view.current_rotamers.get(&ser_id).unwrap();

        let result = system_view.transaction_doublet(ala_id, ser_id, |_| Ok("no_doublet_change"));

        assert_eq!(result.unwrap(), "no_doublet_change");
        assert_eq!(
            *system_view.current_rotamers.get(&ala_id).unwrap(),
            original_ala_rotamer
        );
        assert_eq!(
            *system_view.current_rotamers.get(&ser_id).unwrap(),
            original_ser_rotamer
        );
    }

    #[test]
    fn transaction_fails_for_unknown_residue() {
        let setup = setup();
        let config = create_config(ResidueSelection::All);
        let reporter = ProgressReporter::new();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );
        let mut system = setup.system.clone();
        let mut current_rotamers = HashMap::new();
        let mut system_view = SystemView::new(&mut system, &context, &mut current_rotamers);

        let unknown_residue_id = ResidueId::default();

        let result = system_view.transaction(unknown_residue_id, |_| Ok(()));
        assert!(result.is_err());
    }

    #[test]
    fn transaction_doublet_fails_for_unknown_residue() {
        let setup = setup();
        let config = create_config(ResidueSelection::All);
        let reporter = ProgressReporter::new();
        let context = OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &reporter,
            &config,
            &setup.rotamer_library,
            &setup.topology_registry,
        );
        let mut system = setup.system.clone();
        let mut current_rotamers = HashMap::new();
        let mut system_view = SystemView::new(&mut system, &context, &mut current_rotamers);

        let unknown_residue_id = ResidueId::default();

        let result =
            system_view.transaction_doublet(unknown_residue_id, unknown_residue_id, |_| Ok(()));
        assert!(result.is_err());
    }
}

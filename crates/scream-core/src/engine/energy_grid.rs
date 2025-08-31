use super::cache::ELCache;
use super::error::EngineError;
use super::transaction::SystemView;
use super::utils::query::collect_active_sidechain_atoms;
use crate::core::forcefield::params::Forcefield;
use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::atom::AtomRole;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use itertools::Itertools;
use std::collections::{HashMap, HashSet};
use tracing::{info, trace};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Debug, Clone)]
pub struct EnergyGrid {
    pair_interactions: HashMap<(ResidueId, ResidueId), EnergyTerm>,
    total_residue_interactions: HashMap<ResidueId, EnergyTerm>,
    total_interaction_energy: EnergyTerm,
    current_el_energies: HashMap<ResidueId, EnergyTerm>,
    current_optimization_score: f64,
}

#[derive(Debug, Clone)]
pub struct MoveDelta {
    pub res_id: ResidueId,
    pub new_rotamer_idx: usize,
    pub new_el: EnergyTerm,
    pub new_total_interaction: EnergyTerm,
    pub new_pair_interactions: HashMap<ResidueId, EnergyTerm>,
    pub delta_score: f64,
}

impl EnergyGrid {
    pub fn new(
        system: &MolecularSystem,
        forcefield: &Forcefield,
        active_residues: &HashSet<ResidueId>,
        el_cache: &ELCache,
        initial_rotamers: &HashMap<ResidueId, usize>,
    ) -> Result<Self, EngineError> {
        info!("Initializing EnergyGrid with full energy calculation...");

        let mut pair_interactions = HashMap::new();
        let mut total_residue_interactions: HashMap<ResidueId, EnergyTerm> = active_residues
            .iter()
            .map(|&id| (id, EnergyTerm::default()))
            .collect();

        let scorer = Scorer::new(system, forcefield);

        for pair in active_residues.iter().combinations(2) {
            let res_a_id = *pair[0];
            let res_b_id = *pair[1];

            let atoms_a = collect_active_sidechain_atoms(system, &HashSet::from([res_a_id]));
            let atoms_b = collect_active_sidechain_atoms(system, &HashSet::from([res_b_id]));

            let atoms_a_slice = atoms_a
                .get(&res_a_id)
                .map_or([].as_slice(), |v| v.as_slice());
            let atoms_b_slice = atoms_b
                .get(&res_b_id)
                .map_or([].as_slice(), |v| v.as_slice());

            if atoms_a_slice.is_empty() || atoms_b_slice.is_empty() {
                continue;
            }

            let interaction = scorer.score_interaction(atoms_a_slice, atoms_b_slice)?;

            let key = if res_a_id < res_b_id {
                (res_a_id, res_b_id)
            } else {
                (res_b_id, res_a_id)
            };
            pair_interactions.insert(key, interaction);

            *total_residue_interactions.get_mut(&res_a_id).unwrap() += interaction;
            *total_residue_interactions.get_mut(&res_b_id).unwrap() += interaction;
        }

        let total_interaction_energy = total_residue_interactions
            .values()
            .fold(EnergyTerm::default(), |acc, term| acc + *term)
            * 0.5;

        let mut current_el_energies = HashMap::with_capacity(active_residues.len());
        let mut total_el_energy = EnergyTerm::default();

        for &residue_id in active_residues {
            let residue = system.residue(residue_id).unwrap();
            if let (Some(residue_type), Some(rotamer_idx)) =
                (residue.residue_type, initial_rotamers.get(&residue_id))
            {
                if let Some(el_energy) = el_cache.get(residue_id, residue_type, *rotamer_idx) {
                    current_el_energies.insert(residue_id, *el_energy);
                    total_el_energy += *el_energy;
                } else {
                    current_el_energies.insert(residue_id, EnergyTerm::default());
                }
            }
        }

        let current_optimization_score = total_interaction_energy.total() + total_el_energy.total();

        info!(
            "EnergyGrid initialized. Total optimization score: {:.4}",
            current_optimization_score
        );

        Ok(Self {
            pair_interactions,
            total_residue_interactions,
            total_interaction_energy,
            current_el_energies,
            current_optimization_score,
        })
    }

    pub fn total_score(&self) -> f64 {
        self.current_optimization_score
    }

    pub fn get_el_energy(&self, res_id: ResidueId) -> Option<&EnergyTerm> {
        self.current_el_energies.get(&res_id)
    }

    pub fn calculate_delta_for_move<'a, 'ctx, C>(
        &self,
        res_id: ResidueId,
        new_rotamer_idx: usize,
        system_view: &mut SystemView<'a, 'ctx, C>,
        el_cache: &ELCache,
        active_residues: &HashSet<ResidueId>,
    ) -> Result<MoveDelta, EngineError>
    where
        C: super::context::ProvidesResidueSelections + Sync,
    {
        let (new_total_interaction, new_pair_interactions) =
            system_view.transaction(res_id, |view| {
                view.apply_move(res_id, new_rotamer_idx)?;

                let scorer = Scorer::new(view.system, view.context.forcefield);

                let new_sc_atoms = view
                    .system
                    .residue(res_id)
                    .unwrap()
                    .atoms()
                    .iter()
                    .filter_map(|&id| {
                        view.system.atom(id).and_then(|a| {
                            if a.role == AtomRole::Sidechain {
                                Some(id)
                            } else {
                                None
                            }
                        })
                    })
                    .collect::<Vec<_>>();

                let interaction_sum = EnergyTerm::default();
                let pair_interactions_map = HashMap::new();

                if new_sc_atoms.is_empty() {
                    return Ok((interaction_sum, pair_interactions_map));
                }

                let active_residues_vec: Vec<_> = active_residues
                    .iter()
                    .filter(|&&other_res_id| res_id != other_res_id)
                    .cloned()
                    .collect();

                let compute_interaction = |&other_res_id| {
                    let other_sc_atoms = view
                        .system
                        .residue(other_res_id)
                        .unwrap()
                        .atoms()
                        .iter()
                        .filter_map(|&id| {
                            view.system.atom(id).and_then(|a| {
                                if a.role == AtomRole::Sidechain {
                                    Some(id)
                                } else {
                                    None
                                }
                            })
                        })
                        .collect::<Vec<_>>();

                    if other_sc_atoms.is_empty() {
                        return Ok(None);
                    }

                    scorer
                        .score_interaction(&new_sc_atoms, &other_sc_atoms)
                        .map(|term| Some((other_res_id, term)))
                };

                #[cfg(not(feature = "parallel"))]
                let results: Result<Vec<_>, _> = active_residues_vec
                    .iter()
                    .map(compute_interaction)
                    .collect();

                #[cfg(feature = "parallel")]
                let results: Result<Vec<_>, _> = active_residues_vec
                    .par_iter()
                    .map(compute_interaction)
                    .collect();

                let (interaction_sum, pair_interactions_map) = results?.into_iter().flatten().fold(
                    (EnergyTerm::default(), HashMap::new()),
                    |mut acc, (id, term)| {
                        acc.0 += term;
                        acc.1.insert(id, term);
                        acc
                    },
                );

                Ok((interaction_sum, pair_interactions_map))
            })?;

        let residue_type = system_view
            .system
            .residue(res_id)
            .unwrap()
            .residue_type
            .unwrap();

        let old_el = *self
            .current_el_energies
            .get(&res_id)
            .unwrap_or(&EnergyTerm::default());
        let new_el = *el_cache
            .get(res_id, residue_type, new_rotamer_idx)
            .unwrap_or(&EnergyTerm::default());
        let delta_el = new_el - old_el;

        let old_total_interaction = *self
            .total_residue_interactions
            .get(&res_id)
            .unwrap_or(&EnergyTerm::default());
        let delta_interaction = new_total_interaction - old_total_interaction;
        let delta_score = delta_el.total() + delta_interaction.total();

        trace!(
            "ΔE calculation for res {:?}, rot {}: ΔE_total={:.2} (ΔE_EL={:.2}, ΔE_int={:.2})",
            res_id,
            new_rotamer_idx,
            delta_score,
            delta_el.total(),
            delta_interaction.total()
        );

        Ok(MoveDelta {
            res_id,
            new_rotamer_idx,
            new_el,
            new_total_interaction,
            new_pair_interactions,
            delta_score,
        })
    }

    pub fn apply_move(&mut self, move_delta: MoveDelta) {
        let res_id = move_delta.res_id;

        // 1. Update total score with the pre-calculated delta
        self.current_optimization_score += move_delta.delta_score;

        // 2. Get old interaction terms for diffing
        let _old_el = self
            .current_el_energies
            .get(&res_id)
            .copied()
            .unwrap_or_default();
        let old_total_interaction = self
            .total_residue_interactions
            .get(&res_id)
            .copied()
            .unwrap_or_default();

        // 3. Update EL energy map
        self.current_el_energies.insert(res_id, move_delta.new_el);

        // 4. Update total interaction energy (global sum)
        self.total_interaction_energy = self.total_interaction_energy - old_total_interaction
            + move_delta.new_total_interaction;

        // 5. Update the interaction sum for the moved residue itself
        self.total_residue_interactions
            .insert(res_id, move_delta.new_total_interaction);

        // 6. Update interaction sums for all *other* residues and the pair map
        for (other_res_id, pair_interaction_with_new) in move_delta.new_pair_interactions {
            let key = if res_id < other_res_id {
                (res_id, other_res_id)
            } else {
                (other_res_id, res_id)
            };
            let old_pair_interaction = self
                .pair_interactions
                .get(&key)
                .copied()
                .unwrap_or_default();

            // Update the other residue's total interaction
            if let Some(other_total_int) = self.total_residue_interactions.get_mut(&other_res_id) {
                *other_total_int =
                    *other_total_int - old_pair_interaction + pair_interaction_with_new;
            }

            self.pair_interactions
                .insert(key, pair_interaction_with_new);
        }

        trace!(
            "Applied move for res {:?}. New total score: {:.4}",
            res_id, self.current_optimization_score
        );
    }
}

impl std::ops::Sub for EnergyTerm {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            vdw: self.vdw - rhs.vdw,
            coulomb: self.coulomb - rhs.coulomb,
            hbond: self.hbond - rhs.hbond,
        }
    }
}

impl std::ops::Mul<f64> for EnergyTerm {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            vdw: self.vdw * rhs,
            coulomb: self.coulomb * rhs,
            hbond: self.hbond * rhs,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{
        forcefield::{parameterization::Parameterizer, params::Forcefield, term::EnergyTerm},
        models::{atom::Atom, chain::ChainType, residue::ResidueType, system::MolecularSystem},
        rotamers::{library::RotamerLibrary, rotamer::Rotamer},
        topology::registry::TopologyRegistry,
    };
    use crate::engine::{
        cache::ELCache,
        config::{ConvergenceConfig, PlacementConfigBuilder, ResidueSelection},
        context::OptimizationContext,
        progress::ProgressReporter,
        transaction::SystemView,
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
        el_cache: ELCache,
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

        let mut el_cache = ELCache::new();
        el_cache.insert(
            ala_id,
            ResidueType::Alanine,
            0,
            EnergyTerm::new(1.0, 0.5, 0.0),
        );
        el_cache.insert(
            ala_id,
            ResidueType::Alanine,
            1,
            EnergyTerm::new(1.2, 0.3, 0.0),
        );
        el_cache.insert(
            ser_id,
            ResidueType::Serine,
            0,
            EnergyTerm::new(2.0, 1.0, 0.0),
        );
        el_cache.insert(
            ser_id,
            ResidueType::Serine,
            1,
            EnergyTerm::new(2.1, 0.8, 0.0),
        );

        let active_residues = HashSet::from([ala_id, ser_id]);
        let mut initial_rotamers = HashMap::new();
        initial_rotamers.insert(ala_id, 0);
        initial_rotamers.insert(ser_id, 0);

        TestSetup {
            system,
            forcefield,
            rotamer_library,
            topology_registry,
            el_cache,
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

    fn create_config(selection: ResidueSelection) -> crate::engine::config::PlacementConfig {
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
    fn energy_grid_new_initializes_correctly() {
        let setup = setup();
        let energy_grid = EnergyGrid::new(
            &setup.system,
            &setup.forcefield,
            &setup.active_residues,
            &setup.el_cache,
            &setup.initial_rotamers,
        )
        .unwrap();

        assert_eq!(energy_grid.total_residue_interactions.len(), 2);

        assert!(
            energy_grid
                .current_el_energies
                .contains_key(&setup.active_residues.iter().next().unwrap())
        );

        let total_score = energy_grid.total_score();
        assert!(total_score > 0.0);
    }

    #[test]
    fn energy_grid_total_score_returns_correct_value() {
        let setup = setup();
        let energy_grid = EnergyGrid::new(
            &setup.system,
            &setup.forcefield,
            &setup.active_residues,
            &setup.el_cache,
            &setup.initial_rotamers,
        )
        .unwrap();

        let total_score = energy_grid.total_score();

        let expected_interaction = energy_grid.total_interaction_energy.total();
        let expected_el: f64 = energy_grid
            .current_el_energies
            .values()
            .map(|term| term.total())
            .sum();
        let expected_total = expected_interaction + expected_el;

        assert!((total_score - expected_total).abs() < 1e-9);
    }

    #[test]
    fn calculate_delta_for_move_computes_correct_delta() {
        let setup = setup();
        let energy_grid = EnergyGrid::new(
            &setup.system,
            &setup.forcefield,
            &setup.active_residues,
            &setup.el_cache,
            &setup.initial_rotamers,
        )
        .unwrap();

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

        let move_delta = energy_grid
            .calculate_delta_for_move(
                ala_id,
                1,
                &mut system_view,
                &setup.el_cache,
                &setup.active_residues,
            )
            .unwrap();

        assert_eq!(move_delta.res_id, ala_id);
        assert_eq!(move_delta.new_rotamer_idx, 1);

        assert!(move_delta.delta_score.is_finite());
    }

    #[test]
    fn apply_move_updates_energy_grid_state() {
        let setup = setup();
        let mut energy_grid = EnergyGrid::new(
            &setup.system,
            &setup.forcefield,
            &setup.active_residues,
            &setup.el_cache,
            &setup.initial_rotamers,
        )
        .unwrap();

        let initial_score = energy_grid.total_score();

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

        let move_delta = energy_grid
            .calculate_delta_for_move(
                ala_id,
                1,
                &mut system_view,
                &setup.el_cache,
                &setup.active_residues,
            )
            .unwrap();

        energy_grid.apply_move(move_delta);

        let new_score = energy_grid.total_score();
        assert_ne!(initial_score, new_score);

        let new_el = energy_grid.current_el_energies.get(&ala_id).unwrap();
        let expected_el = setup.el_cache.get(ala_id, ResidueType::Alanine, 1).unwrap();
        assert_eq!(*new_el, *expected_el);
    }

    #[test]
    fn energy_grid_handles_empty_active_residues() {
        let setup = setup();
        let empty_active = HashSet::new();
        let empty_initial = HashMap::new();

        let energy_grid = EnergyGrid::new(
            &setup.system,
            &setup.forcefield,
            &empty_active,
            &setup.el_cache,
            &empty_initial,
        )
        .unwrap();

        assert_eq!(energy_grid.total_residue_interactions.len(), 0);
        assert_eq!(energy_grid.pair_interactions.len(), 0);
        assert_eq!(energy_grid.current_el_energies.len(), 0);
        assert_eq!(energy_grid.total_score(), 0.0);
    }

    #[test]
    fn energy_grid_handles_single_residue() {
        let setup = setup();
        let single_residue = HashSet::from([*setup.active_residues.iter().next().unwrap()]);
        let mut single_initial = HashMap::new();
        single_initial.insert(*single_residue.iter().next().unwrap(), 0);

        let energy_grid = EnergyGrid::new(
            &setup.system,
            &setup.forcefield,
            &single_residue,
            &setup.el_cache,
            &single_initial,
        )
        .unwrap();

        assert_eq!(energy_grid.pair_interactions.len(), 0);
        assert_eq!(energy_grid.total_residue_interactions.len(), 1);

        assert_eq!(energy_grid.current_el_energies.len(), 1);
    }
}

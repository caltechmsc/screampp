use crate::core::forcefield::scoring::Scorer;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::engine::cache::ELCache;
use crate::engine::context::{Context, ProvidesResidueSelections};
use crate::engine::error::EngineError;
use crate::engine::placement;
use tracing::{debug, instrument, trace};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[derive(Debug, Clone, Copy)]
pub struct DoubletResult {
    pub rotamer_idx_a: usize,
    pub rotamer_idx_b: usize,
    pub best_local_energy: f64,
}

#[instrument(skip_all, name = "doublet_optimization_task", fields(res_a = ?res_a_id, res_b = ?res_b_id))]
pub fn run<C: ProvidesResidueSelections + Sync>(
    res_a_id: ResidueId,
    res_b_id: ResidueId,
    system: &crate::core::models::system::MolecularSystem,
    el_cache: &ELCache,
    context: &Context<C>,
) -> Result<DoubletResult, EngineError> {
    let residue_a = system.residue(res_a_id).unwrap();
    let residue_b = system.residue(res_b_id).unwrap();
    let res_type_a = residue_a.res_type.unwrap();
    let res_type_b = residue_b.res_type.unwrap();

    let rotamers_a = context
        .rotamer_library
        .get_rotamers_for(res_type_a)
        .ok_or_else(|| EngineError::RotamerLibrary {
            residue_type: res_type_a.to_string(),
            message: "No rotamers found for doublet optimization.".to_string(),
        })?;

    let rotamers_b = context
        .rotamer_library
        .get_rotamers_for(res_type_b)
        .ok_or_else(|| EngineError::RotamerLibrary {
            residue_type: res_type_b.to_string(),
            message: "No rotamers found for doublet optimization.".to_string(),
        })?;

    let index_pairs: Vec<(usize, usize)> = (0..rotamers_a.len())
        .flat_map(|i| (0..rotamers_b.len()).map(move |j| (i, j)))
        .collect();

    debug!(
        "Optimizing doublet with {}x{} = {} total pairs.",
        rotamers_a.len(),
        rotamers_b.len(),
        index_pairs.len()
    );

    let placement_info_a = context
        .rotamer_library
        .get_placement_info_for(res_type_a)
        .unwrap();
    let placement_info_b = context
        .rotamer_library
        .get_placement_info_for(res_type_b)
        .unwrap();

    #[cfg(not(feature = "parallel"))]
    let iterator = index_pairs.iter();

    #[cfg(feature = "parallel")]
    let iterator = index_pairs.par_iter();

    // Find the pair with the minimum energy
    let best_pair = iterator
        .map(|&(idx_a, idx_b)| {
            let rot_a = &rotamers_a[idx_a];
            let rot_b = &rotamers_b[idx_b];

            let el_a = el_cache
                .get(res_a_id, res_type_a, idx_a)
                .map_or(0.0, |e| e.total());
            let el_b = el_cache
                .get(res_b_id, res_type_b, idx_b)
                .map_or(0.0, |e| e.total());

            let mut temp_system = system.clone();

            placement::place_rotamer_on_system(&mut temp_system, res_a_id, rot_a, placement_info_a)
                .map_err(|e| e)?;
            placement::place_rotamer_on_system(&mut temp_system, res_b_id, rot_b, placement_info_b)
                .map_err(|e| e)?;

            let atoms_a_sidechain_ids: Vec<AtomId> = context
                .rotamer_library
                .get_placement_info_for(res_type_a)
                .unwrap()
                .sidechain_atoms
                .iter()
                .filter_map(|atom_name| {
                    temp_system
                        .residue(res_a_id)
                        .unwrap()
                        .get_atom_id_by_name(atom_name)
                })
                .collect();

            let atoms_b_sidechain_ids: Vec<AtomId> = context
                .rotamer_library
                .get_placement_info_for(res_type_b)
                .unwrap()
                .sidechain_atoms
                .iter()
                .filter_map(|atom_name| {
                    temp_system
                        .residue(res_b_id)
                        .unwrap()
                        .get_atom_id_by_name(atom_name)
                })
                .collect();

            let scorer = Scorer::new(&temp_system, context.forcefield);

            let interaction = scorer
                .score_interaction(&atoms_a_sidechain_ids, &atoms_b_sidechain_ids)?
                .total();

            let local_energy = el_a + el_b + interaction;

            trace!(
                "Pair ({}, {}): E_local = {:.2} (EL_A={:.2}, EL_B={:.2}, Int={:.2})",
                idx_a, idx_b, local_energy, el_a, el_b, interaction
            );

            Ok((local_energy, (idx_a, idx_b)))
        })
        .collect::<Result<Vec<_>, EngineError>>()?
        .into_iter()
        .min_by(|(energy_1, _), (energy_2, _)| {
            energy_1
                .partial_cmp(energy_2)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

    match best_pair {
        Some((best_energy, (best_idx_a, best_idx_b))) => {
            debug!(
                "Found best pair ({}, {}) with local energy {:.2}",
                best_idx_a, best_idx_b, best_energy
            );
            Ok(DoubletResult {
                rotamer_idx_a: best_idx_a,
                rotamer_idx_b: best_idx_b,
                best_local_energy: best_energy,
            })
        }
        None => Err(EngineError::PhaseFailed {
            phase: "Doublet Optimization",
            reason: format!(
                "No rotamer pairs could be evaluated for residues {:?} and {:?}",
                res_a_id, res_b_id
            ),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::params::{Forcefield, GlobalParams, NonBondedParams, VdwParam};
    use crate::core::models::atom::Atom;
    use crate::core::models::chain::ChainType;
    use crate::core::models::residue::ResidueType;
    use crate::core::models::system::MolecularSystem;
    use crate::core::rotamers::library::RotamerLibrary;
    use crate::core::rotamers::placement::PlacementInfo;
    use crate::core::rotamers::rotamer::Rotamer;
    use crate::engine::config::{PlacementConfig, ResidueSelection};
    use crate::engine::progress::ProgressReporter;
    use nalgebra::Point3;
    use std::collections::HashMap;

    fn setup_test_data() -> (
        MolecularSystem,
        ResidueId,
        ResidueId,
        ELCache,
        Forcefield,
        RotamerLibrary,
        PlacementConfig,
        ProgressReporter<'static>,
    ) {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_a_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let res_b_id = system
            .add_residue(chain_id, 2, "LEU", Some(ResidueType::Leucine))
            .unwrap();

        let mut n1 = Atom::new(1, "N", res_a_id, Point3::new(0.0, 1.0, 0.0));
        n1.force_field_type = "N_BB".to_string();
        let mut ca1 = Atom::new(2, "CA", res_a_id, Point3::new(0.0, 0.0, 0.0));
        ca1.force_field_type = "C_BB".to_string();
        let mut c1 = Atom::new(3, "C", res_a_id, Point3::new(1.0, 0.0, 0.0));
        c1.force_field_type = "C_BB".to_string();
        system.add_atom_to_residue(res_a_id, n1).unwrap();
        system.add_atom_to_residue(res_a_id, ca1).unwrap();
        system.add_atom_to_residue(res_a_id, c1).unwrap();

        let mut n2 = Atom::new(4, "N", res_b_id, Point3::new(4.0, 1.0, 0.0));
        n2.force_field_type = "N_BB".to_string();
        let mut ca2 = Atom::new(5, "CA", res_b_id, Point3::new(4.0, 0.0, 0.0));
        ca2.force_field_type = "C_BB".to_string();
        let mut c2 = Atom::new(6, "C", res_b_id, Point3::new(5.0, 0.0, 0.0));
        c2.force_field_type = "C_BB".to_string();
        system.add_atom_to_residue(res_b_id, n2).unwrap();
        system.add_atom_to_residue(res_b_id, ca2).unwrap();
        system.add_atom_to_residue(res_b_id, c2).unwrap();

        let mut vdw = HashMap::new();
        vdw.insert(
            "N_BB".to_string(),
            VdwParam::LennardJones {
                radius: 3.5,
                well_depth: 0.15,
            },
        );
        vdw.insert(
            "C_BB".to_string(),
            VdwParam::LennardJones {
                radius: 3.8,
                well_depth: 0.18,
            },
        );
        vdw.insert(
            "C_CB".to_string(),
            VdwParam::LennardJones {
                radius: 4.0,
                well_depth: 1.0,
            },
        );

        let ff = Forcefield {
            non_bonded: NonBondedParams {
                globals: GlobalParams {
                    dielectric_constant: 1.0,
                    potential_function: "lennard-jones-12-6".to_string(),
                },
                vdw,
                hbond: HashMap::new(),
            },
            deltas: HashMap::new(),
        };

        let common_backbone_atoms = |residue_id: ResidueId| {
            let mut atoms = Vec::new();
            let mut n_atom = Atom::new(1001, "N", residue_id, Point3::new(0.0, 1.0, 0.0));
            n_atom.force_field_type = "N_BB".to_string();
            atoms.push(n_atom);

            let mut ca_atom = Atom::new(1002, "CA", residue_id, Point3::new(0.0, 0.0, 0.0));
            ca_atom.force_field_type = "C_BB".to_string();
            atoms.push(ca_atom);

            let mut c_atom = Atom::new(1003, "C", residue_id, Point3::new(1.0, 0.0, 0.0));
            c_atom.force_field_type = "C_BB".to_string();
            atoms.push(c_atom);
            atoms
        };

        let mut ala_rot0_atoms = common_backbone_atoms(res_a_id);
        let mut ala_cb0 = Atom::new(10, "CB", res_a_id, Point3::new(1.0, 0.0, 0.0));
        ala_cb0.force_field_type = "C_CB".to_string();
        ala_rot0_atoms.push(ala_cb0);
        let rotamer_a0 = Rotamer {
            atoms: ala_rot0_atoms,
        };

        let mut leu_rot0_atoms = common_backbone_atoms(res_b_id);
        let mut leu_cb0 = Atom::new(20, "CB", res_b_id, Point3::new(2.0, 0.0, 0.0));
        leu_cb0.force_field_type = "C_CB".to_string();
        leu_rot0_atoms.push(leu_cb0);
        let rotamer_b0 = Rotamer {
            atoms: leu_rot0_atoms,
        };

        let mut leu_rot1_atoms = common_backbone_atoms(res_b_id);
        let mut leu_cb1 = Atom::new(21, "CB", res_b_id, Point3::new(3.0, 0.0, 0.0));
        leu_cb1.force_field_type = "C_CB".to_string();
        leu_rot1_atoms.push(leu_cb1);
        let rotamer_b1 = Rotamer {
            atoms: leu_rot1_atoms,
        };

        let mut rot_lib_map = HashMap::new();
        rot_lib_map.insert(ResidueType::Alanine, vec![rotamer_a0]);
        rot_lib_map.insert(ResidueType::Leucine, vec![rotamer_b0, rotamer_b1]);

        let placement_info = PlacementInfo {
            anchor_atoms: vec!["N".to_string(), "CA".to_string(), "C".to_string()],
            sidechain_atoms: vec!["CB".to_string()],
            exact_match_atoms: vec![],
            connection_points: vec![],
        };

        let mut placement_map = HashMap::new();
        placement_map.insert(ResidueType::Alanine, placement_info.clone());
        placement_map.insert(ResidueType::Leucine, placement_info);

        let rot_lib = RotamerLibrary {
            rotamers: rot_lib_map,
            placement_info: placement_map,
        };

        let mut el_cache = ELCache::new();
        el_cache.insert(res_a_id, ResidueType::Alanine, 0, Default::default());
        el_cache.insert(res_b_id, ResidueType::Leucine, 0, Default::default());
        el_cache.insert(res_b_id, ResidueType::Leucine, 1, Default::default());

        let config = PlacementConfig {
            residues_to_optimize: ResidueSelection::All,
            scoring: Default::default(),
            optimization: Default::default(),
        };
        let reporter = ProgressReporter::new();

        (
            system, res_a_id, res_b_id, el_cache, ff, rot_lib, config, reporter,
        )
    }

    #[test]
    fn run_finds_unique_optimal_pair() {
        let (system, res_a_id, res_b_id, el_cache, ff, rot_lib, config, reporter) =
            setup_test_data();

        let context = Context::new(&system, &config, &ff, &rot_lib, &reporter);

        let result = run(res_a_id, res_b_id, &system, &el_cache, &context).unwrap();

        assert_eq!(
            result.rotamer_idx_a, 0,
            "Best rotamer for ALA should be index 0"
        );
        assert_eq!(
            result.rotamer_idx_b, 0,
            "Best rotamer for LEU should be index 0"
        );
        assert!(
            (result.best_local_energy - (-0.455568523264)).abs() < 1e-3,
            "Energy was {}, expected around -0.455568523264",
            result.best_local_energy
        );
    }

    #[test]
    fn run_handles_empty_rotamer_list() {
        let (system, res_a_id, res_b_id, el_cache, ff, mut rot_lib, config, reporter) =
            setup_test_data();

        rot_lib
            .rotamers
            .get_mut(&ResidueType::Alanine)
            .unwrap()
            .clear();

        let context = Context::new(&system, &config, &ff, &rot_lib, &reporter);

        let result = run(res_a_id, res_b_id, &system, &el_cache, &context);

        assert!(matches!(result, Err(EngineError::PhaseFailed { .. })));
    }
}

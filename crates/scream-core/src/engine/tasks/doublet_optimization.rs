use crate::core::forcefield::scoring::Scorer;
use crate::core::models::ids::{AtomId, ResidueId};
use crate::core::models::system::MolecularSystem;
use crate::engine::cache::ELCache;
use crate::engine::context::{OptimizationContext, ProvidesResidueSelections};
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
pub fn run<C>(
    res_a_id: ResidueId,
    res_b_id: ResidueId,
    system: &MolecularSystem,
    el_cache: &ELCache,
    context: &OptimizationContext<C>,
) -> Result<DoubletResult, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
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

    if rotamers_a.is_empty() || rotamers_b.is_empty() {
        return Err(EngineError::RotamerLibrary {
            residue_type: if rotamers_a.is_empty() {
                res_type_a.to_string()
            } else {
                res_type_b.to_string()
            },
            message: "Rotamer list is empty, cannot perform doublet optimization.".to_string(),
        });
    }

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
            placement::place_rotamer_on_system(
                &mut temp_system,
                res_a_id,
                rot_a,
                placement_info_a,
            )?;
            placement::place_rotamer_on_system(
                &mut temp_system,
                res_b_id,
                rot_b,
                placement_info_b,
            )?;

            let atoms_a_sidechain_ids: Vec<AtomId> = placement_info_a
                .sidechain_atoms
                .iter()
                .filter_map(|name| {
                    temp_system
                        .residue(res_a_id)
                        .unwrap()
                        .get_atom_id_by_name(name)
                })
                .collect();
            let atoms_b_sidechain_ids: Vec<AtomId> = placement_info_b
                .sidechain_atoms
                .iter()
                .filter_map(|name| {
                    temp_system
                        .residue(res_b_id)
                        .unwrap()
                        .get_atom_id_by_name(name)
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
        Some((best_energy, (best_idx_a, best_idx_b))) => Ok(DoubletResult {
            rotamer_idx_a: best_idx_a,
            rotamer_idx_b: best_idx_b,
            best_local_energy: best_energy,
        }),
        None => Err(EngineError::PhaseFailed {
            phase: "Doublet Optimization",
            reason: "No valid rotamer pairs could be evaluated.".to_string(),
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
    use crate::core::rotamers::library::RotamerLibrary;
    use crate::core::rotamers::placement::PlacementInfo;
    use crate::core::rotamers::rotamer::Rotamer;
    use crate::engine::config::{PlacementConfigBuilder, ResidueSelection};
    use crate::engine::context::OptimizationContext;
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
    ) {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);
        let res_a_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let res_b_id = system
            .add_residue(chain_id, 2, "LEU", Some(ResidueType::Leucine))
            .unwrap();

        let backbone_atoms_data = |residue_id, offset: f64| {
            vec![
                ("N", Point3::new(offset, 1.0, 0.0)),
                ("CA", Point3::new(offset, 0.0, 0.0)),
                ("C", Point3::new(offset + 1.0, 0.0, 0.0)),
            ]
        };

        for (name, pos) in backbone_atoms_data(res_a_id, 0.0) {
            let mut atom = Atom::new(system.atoms_iter().count() + 1, name, res_a_id, pos);
            atom.force_field_type = "BB".to_string();
            system.add_atom_to_residue(res_a_id, atom).unwrap();
        }
        for (name, pos) in backbone_atoms_data(res_b_id, 2.0) {
            let mut atom = Atom::new(system.atoms_iter().count() + 1, name, res_b_id, pos);
            atom.force_field_type = "BB".to_string();
            system.add_atom_to_residue(res_b_id, atom).unwrap();
        }

        let mut vdw = HashMap::new();
        vdw.insert(
            "BB".to_string(),
            VdwParam::LennardJones {
                radius: 1.0,
                well_depth: 0.0,
            },
        );
        vdw.insert(
            "C_SC".to_string(),
            VdwParam::LennardJones {
                radius: 3.8,
                well_depth: 0.1,
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

        let create_rotamer = |residue_id, cb_pos: Point3<f64>| -> Rotamer {
            let mut atoms = Vec::new();
            let mut n = Atom::new(1001, "N", residue_id, Point3::new(0.0, 1.0, 0.0));
            n.force_field_type = "BB".to_string();
            atoms.push(n);
            let mut ca = Atom::new(1002, "CA", residue_id, Point3::new(0.0, 0.0, 0.0));
            ca.force_field_type = "BB".to_string();
            atoms.push(ca);
            let mut c = Atom::new(1003, "C", residue_id, Point3::new(1.0, 0.0, 0.0));
            c.force_field_type = "BB".to_string();
            atoms.push(c);
            let mut cb = Atom::new(10, "CB", residue_id, cb_pos);
            cb.force_field_type = "C_SC".to_string();
            atoms.push(cb);
            Rotamer { atoms }
        };

        let rotamer_a0 = create_rotamer(res_a_id, Point3::new(-0.5, -0.8, 0.0));
        let rotamer_b0 = create_rotamer(res_b_id, Point3::new(0.5, 0.8, 0.0));
        let rotamer_b1 = create_rotamer(res_b_id, Point3::new(-0.5, -0.8, 0.0));

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

        let mut rot_lib = RotamerLibrary {
            rotamers: rot_lib_map,
            placement_info: placement_map,
        };

        let parameterizer =
            crate::core::forcefield::parameterization::Parameterizer::new(ff.clone(), 0.0);
        for rot_vec in rot_lib.rotamers.values_mut() {
            for rot in rot_vec {
                for atom in &mut rot.atoms {
                    parameterizer.parameterize_atom(atom, "ALA").unwrap();
                }
            }
        }

        let mut el_cache = ELCache::new();
        el_cache.insert(res_a_id, ResidueType::Alanine, 0, Default::default());
        el_cache.insert(res_b_id, ResidueType::Leucine, 0, Default::default());
        el_cache.insert(res_b_id, ResidueType::Leucine, 1, Default::default());

        (system, res_a_id, res_b_id, el_cache, ff, rot_lib)
    }

    #[test]
    fn run_finds_optimal_pair_with_less_clash() {
        let (system, res_a_id, res_b_id, el_cache, ff, rot_lib) = setup_test_data();

        let config = PlacementConfigBuilder::new()
            .forcefield_path("dummy.ff")
            .delta_params_path("dummy.delta")
            .s_factor(0.0)
            .rotamer_library_path("dummy.rotlib")
            .placement_registry_path("dummy.reg")
            .max_iterations(1)
            .convergence_threshold(0.1)
            .num_solutions(1)
            .residues_to_optimize(ResidueSelection::All)
            .build()
            .unwrap();
        let reporter = ProgressReporter::new();
        let context = OptimizationContext::new(&system, &ff, &reporter, &config, &rot_lib);

        let result = run(res_a_id, res_b_id, &system, &el_cache, &context).unwrap();

        assert_eq!(result.rotamer_idx_a, 0);
        assert_eq!(
            result.rotamer_idx_b, 0,
            "Expected LEU rotamer 0 to be selected as it results in less steric clash"
        );
    }

    #[test]
    fn run_handles_empty_rotamer_list() {
        let (system, res_a_id, res_b_id, el_cache, ff, mut rot_lib) = setup_test_data();
        rot_lib
            .rotamers
            .get_mut(&ResidueType::Alanine)
            .unwrap()
            .clear();

        let config = PlacementConfigBuilder::new()
            .forcefield_path("dummy.ff")
            .delta_params_path("dummy.delta")
            .s_factor(0.0)
            .rotamer_library_path("dummy.rotlib")
            .placement_registry_path("dummy.reg")
            .max_iterations(1)
            .convergence_threshold(0.1)
            .num_solutions(1)
            .residues_to_optimize(ResidueSelection::All)
            .build()
            .unwrap();

        let reporter = ProgressReporter::new();
        let context = OptimizationContext::new(&system, &ff, &reporter, &config, &rot_lib);

        let result = run(res_a_id, res_b_id, &system, &el_cache, &context);

        assert!(matches!(result, Err(EngineError::RotamerLibrary { .. })));
    }
}

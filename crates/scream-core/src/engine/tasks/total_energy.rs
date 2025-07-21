use crate::core::forcefield::params::Forcefield;
use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::ids::ResidueId;
use crate::core::models::system::MolecularSystem;
use crate::engine::cache::ELCache;
use crate::engine::error::EngineError;
use itertools::Itertools;
use std::collections::{HashMap, HashSet};
use std::iter::Sum;
use tracing::instrument;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

impl Sum for EnergyTerm {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::default(), |acc, term| acc + term)
    }
}

#[instrument(skip_all, name = "total_energy_task")]
pub fn run(
    system: &MolecularSystem,
    forcefield: &Forcefield,
    active_residues: &HashSet<ResidueId>,
    current_rotamers: &HashMap<ResidueId, usize>,
    el_cache: &ELCache,
) -> Result<EnergyTerm, EngineError> {
    let total_el_energy: EnergyTerm = active_residues
        .iter()
        .map(|&residue_id| {
            let residue = system.residue(residue_id).unwrap();
            let res_type = residue.res_type.unwrap();
            let rotamer_idx = current_rotamers.get(&residue_id).unwrap();

            el_cache
                .get(residue_id, res_type, *rotamer_idx)
                .copied()
                .unwrap_or_default()
        })
        .sum();

    let scorer = Scorer::new(system, forcefield);
    let residue_pairs: Vec<_> = active_residues.iter().combinations(2).collect();

    #[cfg(not(feature = "parallel"))]
    let iterator = residue_pairs.iter();

    #[cfg(feature = "parallel")]
    let iterator = residue_pairs.par_iter();

    let interaction_energy: EnergyTerm = iterator
        .map(|pair| {
            let res_a_id = *pair[0];
            let res_b_id = *pair[1];
            let atoms_a = system.residue(res_a_id).unwrap().atoms();
            let atoms_b = system.residue(res_b_id).unwrap().atoms();
            scorer
                .score_interaction(atoms_a, atoms_b)
                .map_err(EngineError::from)
        })
        .try_fold(EnergyTerm::default(), |acc, term| term.map(|t| acc + t))?;

    Ok(total_el_energy + interaction_energy)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::params::{GlobalParams, NonBondedParams, VdwParam};
    use crate::core::models::atom::Atom;
    use crate::core::models::chain::ChainType;
    use crate::core::models::residue::ResidueType;
    use nalgebra::Point3;

    fn setup_test_environment() -> (
        MolecularSystem,
        HashSet<ResidueId>,
        HashMap<ResidueId, usize>,
        ELCache,
        Forcefield,
        ResidueId,
        ResidueId,
    ) {
        let mut system = MolecularSystem::new();
        let chain_id = system.add_chain('A', ChainType::Protein);

        let res1_id = system
            .add_residue(chain_id, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let mut atom1 = Atom::new(1, "CA", res1_id, Point3::new(0.0, 0.0, 0.0));
        atom1.force_field_type = "C".to_string();
        system.add_atom_to_residue(res1_id, atom1).unwrap();

        let res2_id = system
            .add_residue(chain_id, 2, "GLY", Some(ResidueType::Glycine))
            .unwrap();
        let mut atom2 = Atom::new(2, "CA", res2_id, Point3::new(5.0, 0.0, 0.0));
        atom2.force_field_type = "C".to_string();
        system.add_atom_to_residue(res2_id, atom2).unwrap();

        let active_residues: HashSet<ResidueId> = [res1_id, res2_id].iter().cloned().collect();
        let current_rotamers: HashMap<ResidueId, usize> =
            [(res1_id, 0), (res2_id, 1)].iter().cloned().collect();

        let mut el_cache = ELCache::new();
        el_cache.insert(
            res1_id,
            ResidueType::Alanine,
            0,
            EnergyTerm::new(1.0, 2.0, 0.0),
        );
        el_cache.insert(
            res2_id,
            ResidueType::Glycine,
            1,
            EnergyTerm::new(3.0, 4.0, 0.0),
        );

        let mut vdw = HashMap::new();
        vdw.insert(
            "C".to_string(),
            VdwParam::LennardJones {
                radius: 4.0,
                well_depth: 0.1,
            },
        );
        let non_bonded = NonBondedParams {
            globals: GlobalParams {
                dielectric_constant: 1.0,
                potential_function: "lennard-jones-12-6".to_string(),
            },
            vdw,
            hbond: HashMap::new(),
        };
        let forcefield = Forcefield {
            non_bonded,
            deltas: HashMap::new(),
        };

        (
            system,
            active_residues,
            current_rotamers,
            el_cache,
            forcefield,
            res1_id,
            res2_id,
        )
    }

    #[test]
    fn calculates_total_energy_as_sum_of_el_and_interaction() {
        let (system, active_residues, current_rotamers, el_cache, forcefield, res1_id, res2_id) =
            setup_test_environment();

        let total_energy = run(
            &system,
            &forcefield,
            &active_residues,
            &current_rotamers,
            &el_cache,
        )
        .unwrap();

        let el_energy = EnergyTerm::new(1.0, 2.0, 0.0) + EnergyTerm::new(3.0, 4.0, 0.0);

        let scorer = Scorer::new(&system, &forcefield);
        let atoms1 = system.residue(res1_id).unwrap().atoms();
        let atoms2 = system.residue(res2_id).unwrap().atoms();
        let interaction = scorer.score_interaction(atoms1, atoms2).unwrap();

        assert_eq!(total_energy.vdw, el_energy.vdw + interaction.vdw);
        assert_eq!(
            total_energy.coulomb,
            el_energy.coulomb + interaction.coulomb
        );
        assert_eq!(total_energy.hbond, el_energy.hbond + interaction.hbond);
    }

    #[test]
    fn returns_only_el_energy_for_single_active_residue() {
        let (system, _, current_rotamers, el_cache, forcefield, res1_id, _) =
            setup_test_environment();

        let single_active: HashSet<ResidueId> = [res1_id].iter().cloned().collect();

        let total_energy = run(
            &system,
            &forcefield,
            &single_active,
            &current_rotamers,
            &el_cache,
        )
        .unwrap();

        let expected_energy = EnergyTerm::new(1.0, 2.0, 0.0);
        assert_eq!(total_energy, expected_energy);
    }

    #[test]
    fn returns_zero_energy_for_no_active_residues() {
        let (system, _, current_rotamers, el_cache, forcefield, _, _) = setup_test_environment();
        let no_active_residues = HashSet::new();

        let total_energy = run(
            &system,
            &forcefield,
            &no_active_residues,
            &current_rotamers,
            &el_cache,
        )
        .unwrap();

        assert_eq!(total_energy, EnergyTerm::default());
    }

    #[test]
    fn handles_missing_el_energy_in_cache_gracefully() {
        let (system, active_residues, current_rotamers, _, forcefield, res1_id, res2_id) =
            setup_test_environment();

        let empty_el_cache = ELCache::new();

        let total_energy = run(
            &system,
            &forcefield,
            &active_residues,
            &current_rotamers,
            &empty_el_cache,
        )
        .unwrap();

        let scorer = Scorer::new(&system, &forcefield);
        let atoms1 = system.residue(res1_id).unwrap().atoms();
        let atoms2 = system.residue(res2_id).unwrap().atoms();
        let interaction = scorer.score_interaction(atoms1, atoms2).unwrap();

        assert_eq!(total_energy, interaction);
    }

    #[test]
    fn propagates_scoring_error() {
        let (mut system, active_residues, current_rotamers, el_cache, forcefield, res1_id, _) =
            setup_test_environment();

        let atom1_id = system.residue(res1_id).unwrap().atoms()[0];
        system.atom_mut(atom1_id).unwrap().force_field_type = "UNKNOWN".to_string();

        let result = run(
            &system,
            &forcefield,
            &active_residues,
            &current_rotamers,
            &el_cache,
        );

        assert!(matches!(result, Err(EngineError::Scoring { .. })));
    }

    #[test]
    fn calculates_interaction_energy_for_multiple_pairs() {
        let (mut system, _, _, mut el_cache, forcefield, res1_id, res2_id) =
            setup_test_environment();

        let chain_id = system.find_chain_by_id('A').unwrap();
        let res3_id = system
            .add_residue(chain_id, 3, "LEU", Some(ResidueType::Leucine))
            .unwrap();
        let mut atom3 = Atom::new(3, "CA", res3_id, Point3::new(0.0, 5.0, 0.0));
        atom3.force_field_type = "C".to_string();
        system.add_atom_to_residue(res3_id, atom3).unwrap();

        let active_residues: HashSet<ResidueId> =
            [res1_id, res2_id, res3_id].iter().cloned().collect();
        let current_rotamers: HashMap<ResidueId, usize> =
            [(res1_id, 0), (res2_id, 1), (res3_id, 0)]
                .iter()
                .cloned()
                .collect();
        el_cache.insert(
            res3_id,
            ResidueType::Leucine,
            0,
            EnergyTerm::new(5.0, 6.0, 0.0),
        );

        let total_energy = run(
            &system,
            &forcefield,
            &active_residues,
            &current_rotamers,
            &el_cache,
        )
        .unwrap();

        let scorer = Scorer::new(&system, &forcefield);
        let atoms1 = system.residue(res1_id).unwrap().atoms();
        let atoms2 = system.residue(res2_id).unwrap().atoms();
        let atoms3 = system.residue(res3_id).unwrap().atoms();

        let interaction12 = scorer.score_interaction(atoms1, atoms2).unwrap();
        let interaction13 = scorer.score_interaction(atoms1, atoms3).unwrap();
        let interaction23 = scorer.score_interaction(atoms2, atoms3).unwrap();
        let total_interaction = interaction12 + interaction13 + interaction23;

        let el_energy = el_cache
            .get(res1_id, ResidueType::Alanine, 0)
            .unwrap()
            .clone()
            + el_cache
                .get(res2_id, ResidueType::Glycine, 1)
                .unwrap()
                .clone()
            + el_cache
                .get(res3_id, ResidueType::Leucine, 0)
                .unwrap()
                .clone();

        let expected_total_energy = el_energy + total_interaction;

        assert!((total_energy.vdw - expected_total_energy.vdw).abs() < 1e-9);
        assert!((total_energy.coulomb - expected_total_energy.coulomb).abs() < 1e-9);
        assert!((total_energy.hbond - expected_total_energy.hbond).abs() < 1e-9);
    }

    #[test]
    #[should_panic]
    fn panics_if_rotamer_index_is_missing_for_active_residue() {
        let (system, active_residues, _, el_cache, forcefield, _, _) = setup_test_environment();
        let mut incomplete_rotamers = HashMap::new();
        incomplete_rotamers.insert(*active_residues.iter().next().unwrap(), 0);

        let _ = run(
            &system,
            &forcefield,
            &active_residues,
            &incomplete_rotamers,
            &el_cache,
        );
    }
}

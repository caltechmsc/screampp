use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::atom::AtomRole;
use crate::core::models::ids::AtomId;
use crate::engine::context::{OptimizationContext, ProvidesResidueSelections};
use crate::engine::error::EngineError;
use tracing::{info, instrument};

struct FixedAtomSets {
    sidechain_atoms: Vec<AtomId>,
    non_sidechain_atoms: Vec<AtomId>,
}

#[instrument(skip_all, name = "fixed_energy_task")]
pub fn run<C>(context: &OptimizationContext<C>) -> Result<EnergyTerm, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    info!("Calculating constant energy offset for the fixed parts of the system.");

    let fixed_sets = collect_fixed_atom_sets(context)?;

    if fixed_sets.sidechain_atoms.is_empty() && fixed_sets.non_sidechain_atoms.is_empty() {
        info!("No fixed atoms found in the system; the energy offset is zero.");
        return Ok(EnergyTerm::default());
    }

    let scorer = Scorer::new(context.system, context.forcefield);

    let energy1 = scorer.score_group_internal(&fixed_sets.non_sidechain_atoms)?;

    let energy2 =
        scorer.score_interaction(&fixed_sets.sidechain_atoms, &fixed_sets.non_sidechain_atoms)?;

    // TODO: Implement summation of `E_internal_or_fallback` for `fixed_sets.sidechain_atoms`.
    let energy3 = EnergyTerm::default();

    let total_offset_energy = energy1 + energy2 + energy3;
    info!(
        energy = total_offset_energy.total(),
        vdw = total_offset_energy.vdw,
        coulomb = total_offset_energy.coulomb,
        hbond = total_offset_energy.hbond,
        "Fixed energy offset calculation complete."
    );

    Ok(total_offset_energy)
}

fn collect_fixed_atom_sets<C>(
    context: &OptimizationContext<C>,
) -> Result<FixedAtomSets, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    let active_residues = context.resolve_all_active_residues()?;

    let mut sidechain_atoms = Vec::new();
    let mut non_sidechain_atoms = Vec::new();

    for (atom_id, atom) in context.system.atoms_iter() {
        if active_residues.contains(&atom.residue_id) {
            if atom.role == AtomRole::Backbone {
                non_sidechain_atoms.push(atom_id);
            }
        } else {
            match atom.role {
                AtomRole::Sidechain => sidechain_atoms.push(atom_id),
                _ => non_sidechain_atoms.push(atom_id),
            }
        }
    }

    Ok(FixedAtomSets {
        sidechain_atoms,
        non_sidechain_atoms,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::parameterization::Parameterizer;
    use crate::core::forcefield::params::Forcefield;
    use crate::core::forcefield::scoring::Scorer;
    use crate::core::models::atom::Atom;
    use crate::core::models::chain::ChainType;
    use crate::core::models::ids::{AtomId, ResidueId};
    use crate::core::models::residue::ResidueType;
    use crate::core::models::system::MolecularSystem;
    use crate::core::models::topology::BondOrder;
    use crate::core::rotamers::library::RotamerLibrary;
    use crate::core::rotamers::rotamer::Rotamer;
    use crate::core::topology::registry::TopologyRegistry;
    use crate::engine::config::{
        ConvergenceConfig, PlacementConfig, PlacementConfigBuilder, ResidueSelection,
        ResidueSpecifier,
    };
    use crate::engine::progress::ProgressReporter;
    use nalgebra::Point3;
    use std::collections::HashSet;
    use std::fs::File;
    use std::io::Write;
    use tempfile::TempDir;

    struct TestSetup {
        _tmp: TempDir,
        pub system: MolecularSystem,
        pub forcefield: Forcefield,
        pub topology_registry: TopologyRegistry,
        pub rotamer_library: RotamerLibrary,
        pub config: PlacementConfig,
        pub reporter: ProgressReporter<'static>,
        pub res_ala_id: ResidueId,
        pub res_ser_id: ResidueId,
        pub res_hoh_id: ResidueId,
        pub res_lig_id: ResidueId,
    }

    fn add_atom_to_res(
        system: &mut MolecularSystem,
        res_id: ResidueId,
        name: &str,
        pos: [f64; 3],
    ) -> AtomId {
        let atom = Atom::new(name, res_id, Point3::new(pos[0], pos[1], pos[2]));
        system.add_atom_to_residue(res_id, atom).unwrap()
    }

    fn write_test_files(
        dir: &TempDir,
    ) -> (std::path::PathBuf, std::path::PathBuf, std::path::PathBuf) {
        let dirp = dir.path();

        let ff_path = dirp.join("ff.toml");
        let mut ff = File::create(&ff_path).unwrap();
        write!(
            ff,
            r#"
            [globals]
            dielectric_constant = 4.0
            potential_function = "lennard-jones-12-6"

            [vdw]
            N_R   = {{ radius = 1.6, well_depth = 0.10 }}
            C_BB  = {{ radius = 1.8, well_depth = 0.10 }}
            C_SC  = {{ radius = 1.9, well_depth = 0.12 }}
            C_R   = {{ radius = 1.7, well_depth = 0.09 }}
            O_2   = {{ radius = 1.6, well_depth = 0.15 }}
            O_HOH = {{ radius = 1.6, well_depth = 0.15 }}
            C_LIG = {{ radius = 1.75, well_depth = 0.11 }}

            [hbond]
            O_H-O_2 = {{ equilibrium_distance = 2.8, well_depth = 5.0 }}
            "#
        )
        .unwrap();

        let delta_path = dirp.join("delta.csv");
        let mut delta = File::create(&delta_path).unwrap();
        writeln!(delta, "residue_type,atom_name,mu,sigma").unwrap();

        let topo_path = dirp.join("topo.toml");
        let mut topo = File::create(&topo_path).unwrap();
        write!(
            topo,
            r#"
            [ALA]
            anchor_atoms = ["N", "CA", "C"]
            sidechain_atoms = ["CB"]

            [SER]
            anchor_atoms = ["N", "CA", "C"]
            sidechain_atoms = ["OG"]
            "#
        )
        .unwrap();

        (ff_path, delta_path, topo_path)
    }

    fn build_placement_config(
        ff_path: &std::path::Path,
        delta_path: &std::path::Path,
        topo_path: &std::path::Path,
    ) -> PlacementConfig {
        PlacementConfigBuilder::new()
            .forcefield_path(ff_path)
            .delta_params_path(delta_path)
            .rotamer_library_path(topo_path)
            .topology_registry_path(topo_path)
            .s_factor(0.0)
            .max_iterations(1)
            .num_solutions(1)
            .include_input_conformation(false)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.01,
                patience_iterations: 10,
            })
            .simulated_annealing_config(None)
            .final_refinement_iterations(0)
            .residues_to_optimize(ResidueSelection::List {
                include: vec![ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: 1,
                }],
                exclude: vec![],
            })
            .build()
            .unwrap()
    }

    fn setup() -> TestSetup {
        let tmp = tempfile::tempdir().unwrap();
        let (ff_path, delta_path, topo_path) = write_test_files(&tmp);
        let forcefield = Forcefield::load(&ff_path, &delta_path).unwrap();
        let topology_registry = TopologyRegistry::load(&topo_path).unwrap();

        let mut system = MolecularSystem::new();
        let chain_a = system.add_chain('A', ChainType::Protein);
        let res_ala_id = system
            .add_residue(chain_a, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();
        let res_ser_id = system
            .add_residue(chain_a, 2, "SER", Some(ResidueType::Serine))
            .unwrap();

        let ala_n = add_atom_to_res(&mut system, res_ala_id, "N", [-1.2, 0.5, 0.0]);
        let ala_ca = add_atom_to_res(&mut system, res_ala_id, "CA", [0.0, 0.0, 0.0]);
        let ala_c = add_atom_to_res(&mut system, res_ala_id, "C", [0.8, -0.8, 0.0]);
        let ala_cb = add_atom_to_res(&mut system, res_ala_id, "CB", [-0.5, -1.0, 1.2]);

        let ser_n = add_atom_to_res(&mut system, res_ser_id, "N", [5.0, 0.0, 0.0]);
        let ser_ca = add_atom_to_res(&mut system, res_ser_id, "CA", [6.4, 0.0, 0.0]);
        let ser_c = add_atom_to_res(&mut system, res_ser_id, "C", [7.8, 0.0, 0.0]);
        let ser_og = add_atom_to_res(&mut system, res_ser_id, "OG", [6.4, 1.3, 0.0]);

        let chain_w = system.add_chain('W', ChainType::Water);
        let res_hoh_id = system.add_residue(chain_w, 101, "HOH", None).unwrap();
        let hoh_o = add_atom_to_res(&mut system, res_hoh_id, "O", [12.0, 0.0, 0.0]);

        let chain_l = system.add_chain('L', ChainType::Ligand);
        let res_lig_id = system.add_residue(chain_l, 201, "LIG", None).unwrap();
        let lig_c1 = add_atom_to_res(&mut system, res_lig_id, "C1", [15.0, -1.0, 0.0]);

        let mut set = |id: crate::core::models::ids::AtomId, ff: &str, q: f64| {
            let a = system.atom_mut(id).unwrap();
            a.force_field_type = ff.to_string();
            a.partial_charge = q;
        };
        set(ala_n, "N_R", -0.3);
        set(ala_ca, "C_BB", 0.1);
        set(ala_c, "C_BB", 0.5);
        set(ala_cb, "C_SC", -0.1);

        set(ser_n, "N_R", -0.3);
        set(ser_ca, "C_BB", 0.1);
        set(ser_c, "C_BB", 0.4);
        set(ser_og, "O_2", -0.5);

        set(hoh_o, "O_HOH", -0.4);
        set(lig_c1, "C_LIG", -0.1);

        system.add_bond(ala_n, ala_ca, BondOrder::Single).unwrap();
        system.add_bond(ala_ca, ala_c, BondOrder::Single).unwrap();
        system.add_bond(ala_ca, ala_cb, BondOrder::Single).unwrap();

        system.add_bond(ser_n, ser_ca, BondOrder::Single).unwrap();
        system.add_bond(ser_ca, ser_c, BondOrder::Single).unwrap();
        system.add_bond(ser_ca, ser_og, BondOrder::Single).unwrap();

        let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 0.0);
        parameterizer.parameterize_system(&mut system).unwrap();

        let mut rotamer_library = RotamerLibrary::default();
        rotamer_library.rotamers.insert(
            ResidueType::Alanine,
            vec![Rotamer {
                atoms: vec![],
                bonds: vec![],
                energy: Default::default(),
            }],
        );

        let config = build_placement_config(&ff_path, &delta_path, &topo_path);
        let reporter = ProgressReporter::default();

        TestSetup {
            _tmp: tmp,
            system,
            forcefield,
            topology_registry,
            rotamer_library,
            config,
            reporter,
            res_ala_id,
            res_ser_id,
            res_hoh_id,
            res_lig_id,
        }
    }

    fn build_context<'a>(setup: &'a TestSetup) -> OptimizationContext<'a, PlacementConfig> {
        OptimizationContext::new(
            &setup.system,
            &setup.forcefield,
            &setup.reporter,
            &setup.config,
            &setup.rotamer_library,
            &setup.topology_registry,
        )
    }

    #[test]
    fn collects_correct_atoms_into_disjoint_sets() {
        let setup = setup();
        let context = build_context(&setup);

        let fixed = collect_fixed_atom_sets(&context).expect("collect_fixed_atom_sets should work");

        let res = |res_id| setup.system.residue(res_id).unwrap();
        let id_of = |res_id, name: &str| res(res_id).get_first_atom_id_by_name(name).unwrap();

        let expected_sc: HashSet<AtomId> =
            vec![id_of(setup.res_ser_id, "OG")].into_iter().collect();

        let mut expected_non_sc: HashSet<AtomId> = HashSet::new();
        for name in ["N", "CA", "C"].iter() {
            expected_non_sc.insert(id_of(setup.res_ala_id, name));
            expected_non_sc.insert(id_of(setup.res_ser_id, name));
        }
        expected_non_sc.insert(id_of(setup.res_hoh_id, "O"));
        expected_non_sc.insert(id_of(setup.res_lig_id, "C1"));

        let sc_set: HashSet<AtomId> = fixed.sidechain_atoms.iter().copied().collect();
        let non_sc_set: HashSet<AtomId> = fixed.non_sidechain_atoms.iter().copied().collect();

        assert_eq!(sc_set, expected_sc, "Sc_F should contain only SER.OG");
        assert_eq!(non_sc_set, expected_non_sc, "NonSc_F membership mismatch");

        let active = context.resolve_all_active_residues().unwrap();
        let mut expected_total = 0usize;
        for (_aid, atom) in setup.system.atoms_iter() {
            if active.contains(&atom.residue_id) {
                if atom.role == AtomRole::Backbone {
                    expected_total += 1;
                }
            } else {
                expected_total += 1;
            }
        }
        assert_eq!(
            fixed.sidechain_atoms.len() + fixed.non_sidechain_atoms.len(),
            expected_total
        );

        let intersection: HashSet<_> = sc_set.intersection(&non_sc_set).copied().collect();
        assert!(intersection.is_empty(), "Sets should be disjoint");
    }

    #[test]
    fn run_returns_backbone_internal_energy_when_only_active_ala_present() {
        let tmp = tempfile::tempdir().unwrap();
        let (ff_path, delta_path, topo_path) = write_test_files(&tmp);
        let forcefield = Forcefield::load(&ff_path, &delta_path).unwrap();
        let topology_registry = TopologyRegistry::load(&topo_path).unwrap();

        let mut system = MolecularSystem::new();
        let chain_a = system.add_chain('A', ChainType::Protein);
        let res_ala_id = system
            .add_residue(chain_a, 1, "ALA", Some(ResidueType::Alanine))
            .unwrap();

        let n = add_atom_to_res(&mut system, res_ala_id, "N", [-1.2, 0.5, 0.0]);
        let ca = add_atom_to_res(&mut system, res_ala_id, "CA", [0.0, 0.0, 0.0]);
        let c = add_atom_to_res(&mut system, res_ala_id, "C", [0.8, -0.8, 0.0]);
        let o = add_atom_to_res(&mut system, res_ala_id, "O", [1.5, -1.2, 0.0]);

        let mut set = |id: AtomId, ff: &str, q: f64| {
            let a = system.atom_mut(id).unwrap();
            a.force_field_type = ff.to_string();
            a.partial_charge = q;
        };
        set(n, "N_R", -0.3);
        set(ca, "C_BB", 0.1);
        set(c, "C_BB", 0.5);
        set(o, "O_2", -0.5);

        system.add_bond(n, ca, BondOrder::Single).unwrap();
        system.add_bond(ca, c, BondOrder::Single).unwrap();
        system.add_bond(c, o, BondOrder::Single).unwrap();

        let parameterizer = Parameterizer::new(&forcefield, &topology_registry, 0.0);
        parameterizer.parameterize_system(&mut system).unwrap();

        let mut rotamer_library = RotamerLibrary::default();
        rotamer_library.rotamers.insert(
            ResidueType::Alanine,
            vec![Rotamer {
                atoms: vec![],
                bonds: vec![],
                energy: Default::default(),
            }],
        );
        let reporter = ProgressReporter::default();
        let config = PlacementConfigBuilder::new()
            .forcefield_path(&ff_path)
            .delta_params_path(&delta_path)
            .rotamer_library_path(&topo_path)
            .topology_registry_path(&topo_path)
            .s_factor(0.0)
            .max_iterations(1)
            .num_solutions(1)
            .include_input_conformation(false)
            .convergence_config(ConvergenceConfig {
                energy_threshold: 0.01,
                patience_iterations: 10,
            })
            .final_refinement_iterations(0)
            .residues_to_optimize(ResidueSelection::List {
                include: vec![ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: 1,
                }],
                exclude: vec![],
            })
            .build()
            .unwrap();

        let context = OptimizationContext::new(
            &system,
            &forcefield,
            &reporter,
            &config,
            &rotamer_library,
            &topology_registry,
        );

        let fixed = collect_fixed_atom_sets(&context).unwrap();
        assert!(fixed.sidechain_atoms.is_empty());

        let scorer = Scorer::new(&system, &forcefield);
        let expected = scorer
            .score_group_internal(&fixed.non_sidechain_atoms)
            .unwrap();
        let total = run(&context).unwrap();

        assert!((total.total() - expected.total()).abs() < 1e-9);
        assert!(
            expected.total().abs() > 1e-9,
            "Expected non-zero internal energy for ALA backbone with O present"
        );
    }

    #[test]
    fn run_calculates_correct_energy_components() {
        let setup = setup();
        let context = build_context(&setup);
        let fixed = collect_fixed_atom_sets(&context).unwrap();

        let scorer = Scorer::new(&setup.system, &setup.forcefield);
        let expected_energy1 = scorer
            .score_group_internal(&fixed.non_sidechain_atoms)
            .unwrap();
        let expected_energy2 = scorer
            .score_interaction(&fixed.sidechain_atoms, &fixed.non_sidechain_atoms)
            .unwrap();

        let total = run(&context).unwrap();
        let expected_total = expected_energy1.total() + expected_energy2.total();
        assert!((total.total() - expected_total).abs() < 1e-9);
    }

    #[test]
    fn run_integration_with_real_values() {
        let setup = setup();
        let context = build_context(&setup);
        let fixed = collect_fixed_atom_sets(&context).unwrap();
        let mut union_ids = fixed.non_sidechain_atoms.clone();
        union_ids.extend(fixed.sidechain_atoms.iter().copied());

        let scorer = Scorer::new(&setup.system, &setup.forcefield);
        let expected_total_fixed_energy = scorer.score_group_internal(&union_ids).unwrap();

        let total = run(&context).unwrap();

        assert!((total.total() - expected_total_fixed_energy.total()).abs() < 1e-12);
    }
}

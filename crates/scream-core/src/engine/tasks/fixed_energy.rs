use crate::core::forcefield::scoring::Scorer;
use crate::core::forcefield::term::EnergyTerm;
use crate::core::models::atom::AtomRole;
use crate::core::models::ids::AtomId;
use crate::engine::context::{OptimizationContext, ProvidesResidueSelections};
use crate::engine::error::EngineError;
use tracing::{info, instrument};

#[instrument(skip_all, name = "fixed_energy_task")]
pub fn run<C>(context: &OptimizationContext<C>) -> Result<EnergyTerm, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    info!("Calculating constant energy offset for the fixed parts of the system.");

    // Step 1: Collect all fixed atoms into a single group.
    let fixed_atom_ids = collect_fixed_atom_ids(context)?;

    if fixed_atom_ids.is_empty() {
        info!("No fixed atoms found; the energy offset is zero.");
        return Ok(EnergyTerm::default());
    }

    let scorer = Scorer::new(context.system, context.forcefield);

    // Step 2: Calculate the total internal non-bonded energy of this single fixed group.
    let total_offset_energy = scorer.score_group_internal(&fixed_atom_ids)?;

    info!(
        "Fixed energy offset calculation complete. Energy = {:.4} kcal/mol (VDW: {:.4}, Coulomb: {:.4}, H-Bond: {:.4})",
        total_offset_energy.total(),
        total_offset_energy.vdw,
        total_offset_energy.coulomb,
        total_offset_energy.hbond
    );

    Ok(total_offset_energy)
}

fn collect_fixed_atom_ids<C>(context: &OptimizationContext<C>) -> Result<Vec<AtomId>, EngineError>
where
    C: ProvidesResidueSelections + Sync,
{
    let active_residues = context.resolve_all_active_residues()?;

    let fixed_atom_ids = context
        .system
        .atoms_iter()
        .filter_map(|(atom_id, atom)| {
            if !active_residues.contains(&atom.residue_id) || atom.role == AtomRole::Backbone {
                Some(atom_id)
            } else {
                None
            }
        })
        .collect();

    Ok(fixed_atom_ids)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::forcefield::parameterization::Parameterizer;
    use crate::core::forcefield::params::{EnergyWeights, Forcefield};
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
        let forcefield =
            Forcefield::load(&ff_path, &delta_path, &EnergyWeights::default()).unwrap();
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
    fn collects_correct_set_of_fixed_atoms() {
        let setup = setup();
        let context = build_context(&setup);

        let fixed_ids =
            collect_fixed_atom_ids(&context).expect("collect_fixed_atom_ids should work");

        let fixed_ids_set: HashSet<_> = fixed_ids.iter().cloned().collect();

        let res = |res_id| setup.system.residue(res_id).unwrap();
        let id_of = |res_id, name: &str| res(res_id).get_first_atom_id_by_name(name).unwrap();

        for atom_id in res(setup.res_ser_id).atoms() {
            assert!(fixed_ids_set.contains(atom_id));
        }
        for atom_id in res(setup.res_hoh_id).atoms() {
            assert!(fixed_ids_set.contains(atom_id));
        }
        for atom_id in res(setup.res_lig_id).atoms() {
            assert!(fixed_ids_set.contains(atom_id));
        }

        let ala_bb_ids: HashSet<_> = ["N", "CA", "C"]
            .iter()
            .map(|name| id_of(setup.res_ala_id, name))
            .collect();
        let ala_sc_id = id_of(setup.res_ala_id, "CB");

        for bb_id in &ala_bb_ids {
            assert!(fixed_ids_set.contains(bb_id));
        }
        assert!(!fixed_ids_set.contains(&ala_sc_id));

        let total_fixed_count = res(setup.res_ser_id).atoms().len()
            + res(setup.res_hoh_id).atoms().len()
            + res(setup.res_lig_id).atoms().len()
            + ala_bb_ids.len();
        assert_eq!(fixed_ids.len(), total_fixed_count);
    }

    #[test]
    fn run_calculates_correct_energy_for_all_fixed_atoms() {
        let setup = setup();
        let context = build_context(&setup);

        let fixed_atom_ids = collect_fixed_atom_ids(&context).unwrap();
        let scorer = Scorer::new(&setup.system, &setup.forcefield);
        let expected_energy = scorer.score_group_internal(&fixed_atom_ids).unwrap();

        let calculated_energy = run(&context).unwrap();

        assert!((calculated_energy.total() - expected_energy.total()).abs() < 1e-9);
        assert!(
            expected_energy.total().abs() > 1e-9,
            "Expected a non-zero energy for the fixed system"
        );
    }

    #[test]
    fn run_handles_system_with_only_active_residues() {
        let mut setup = setup();
        setup.config.residues_to_optimize = ResidueSelection::List {
            include: vec![
                ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: 1,
                },
                ResidueSpecifier {
                    chain_id: 'A',
                    residue_number: 2,
                },
            ],
            exclude: vec![],
        };
        let context = build_context(&setup);

        let fixed_atom_ids = collect_fixed_atom_ids(&context).unwrap();
        let scorer = Scorer::new(&setup.system, &setup.forcefield);
        let expected_energy = scorer.score_group_internal(&fixed_atom_ids).unwrap();

        let calculated_energy = run(&context).unwrap();

        assert!((calculated_energy.total() - expected_energy.total()).abs() < 1e-9);
    }
}

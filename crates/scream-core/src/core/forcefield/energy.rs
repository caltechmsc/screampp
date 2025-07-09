use crate::core::forcefield::params::VdwParam;
use crate::core::forcefield::potentials;
use crate::core::models::atom::Atom;
use std::f64::consts::PI;

pub struct EnergyCalculator;

impl EnergyCalculator {
    pub fn calculate_vdw(
        atom1: &Atom,
        atom2: &Atom,
        vdw_param1: &VdwParam,
        vdw_param2: &VdwParam,
    ) -> f64 {
        let dist = (atom1.position - atom2.position).norm();

        let total_delta = (atom1.delta.powi(2) + atom2.delta.powi(2)).sqrt();

        let (r_min1, well_depth1) = match vdw_param1 {
            VdwParam::LennardJones { radius, well_depth } => (*radius, *well_depth),
            VdwParam::Buckingham {
                radius, well_depth, ..
            } => (*radius, *well_depth),
        };
        let (r_min2, well_depth2) = match vdw_param2 {
            VdwParam::LennardJones { radius, well_depth } => (*radius, *well_depth),
            VdwParam::Buckingham {
                radius, well_depth, ..
            } => (*radius, *well_depth),
        };

        let r_min_combined = (r_min1 * r_min2).sqrt();
        let well_depth_combined = (well_depth1 * well_depth2).sqrt();

        let base_potential_fn = |d: f64| -> f64 {
            match (vdw_param1, vdw_param2) {
                (
                    VdwParam::Buckingham { scale: s1, .. },
                    VdwParam::Buckingham { scale: s2, .. },
                ) => potentials::buckingham_exp_6(
                    d,
                    r_min_combined,
                    well_depth_combined,
                    (s1 + s2) / 2.0,
                ),
                (VdwParam::Buckingham { scale, .. }, _)
                | (_, VdwParam::Buckingham { scale, .. }) => {
                    potentials::buckingham_exp_6(d, r_min_combined, well_depth_combined, *scale)
                }

                _ => potentials::lennard_jones_12_6(d, r_min_combined, well_depth_combined),
            }
        };

        potentials::apply_flat_bottom_vdw(dist, r_min_combined, total_delta, base_potential_fn)
    }

    pub fn calculate_coulomb(atom1: &Atom, atom2: &Atom, dielectric: f64) -> f64 {
        let dist = (atom1.position - atom2.position).norm();
        potentials::coulomb(dist, atom1.partial_charge, atom2.partial_charge, dielectric)
    }

    pub fn calculate_hbond(
        acceptor: &Atom,
        hydrogen: &Atom,
        donor: &Atom,
        r_hb: f64,
        d_hb: f64,
    ) -> f64 {
        let dist_ad = (acceptor.position - donor.position).norm();

        let v_ha = acceptor.position - hydrogen.position;
        let v_hd = donor.position - hydrogen.position;
        let angle_ahd_deg = (v_ha.angle(&v_hd)).to_degrees();

        let total_delta = (acceptor.delta.powi(2) + donor.delta.powi(2)).sqrt();

        let distance_potential_fn =
            |d: f64| -> f64 { potentials::dreiding_hbond_12_10(d, r_hb, d_hb) };

        let flat_bottom_distance_energy =
            potentials::apply_flat_bottom_hbond(dist_ad, r_hb, total_delta, distance_potential_fn);

        if angle_ahd_deg < 90.0 {
            return 0.0;
        }
        let cos_theta = (angle_ahd_deg * PI / 180.0).cos();
        let angular_term = cos_theta.powi(4);

        flat_bottom_distance_energy * angular_term
    }
}

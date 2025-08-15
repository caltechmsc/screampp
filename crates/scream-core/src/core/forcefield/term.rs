use std::ops::{Add, AddAssign};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct EnergyTerm {
    pub vdw: f64,
    pub coulomb: f64,
    pub hbond: f64,
    pub bond: f64,
    pub angle: f64,
    pub torsion: f64,
    pub inversion: f64,
}

impl EnergyTerm {
    pub fn new(
        vdw: f64,
        coulomb: f64,
        hbond: f64,
        bond: f64,
        angle: f64,
        torsion: f64,
        inversion: f64,
    ) -> Self {
        Self {
            vdw,
            coulomb,
            hbond,
            bond,
            angle,
            torsion,
            inversion,
        }
    }

    pub fn from_nonbonded(vdw: f64, coulomb: f64, hbond: f64) -> Self {
        Self {
            vdw,
            coulomb,
            hbond,
            bond: 0.0,
            angle: 0.0,
            torsion: 0.0,
            inversion: 0.0,
        }
    }

    #[inline]
    pub fn total(&self) -> f64 {
        self.vdw
            + self.coulomb
            + self.hbond
            + self.bond
            + self.angle
            + self.torsion
            + self.inversion
    }
}

impl Add for EnergyTerm {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            vdw: self.vdw + rhs.vdw,
            coulomb: self.coulomb + rhs.coulomb,
            hbond: self.hbond + rhs.hbond,
            bond: self.bond + rhs.bond,
            angle: self.angle + rhs.angle,
            torsion: self.torsion + rhs.torsion,
            inversion: self.inversion + rhs.inversion,
        }
    }
}

impl AddAssign for EnergyTerm {
    fn add_assign(&mut self, rhs: Self) {
        self.vdw += rhs.vdw;
        self.coulomb += rhs.coulomb;
        self.hbond += rhs.hbond;
        self.bond += rhs.bond;
        self.angle += rhs.angle;
        self.torsion += rhs.torsion;
        self.inversion += rhs.inversion;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_creates_energy_term_with_specified_values() {
        let term = EnergyTerm::from_nonbonded(1.0, 2.0, 3.0);
        assert_eq!(term.vdw, 1.0);
        assert_eq!(term.coulomb, 2.0);
        assert_eq!(term.hbond, 3.0);
    }

    #[test]
    fn total_returns_sum_of_all_terms() {
        let term = EnergyTerm::from_nonbonded(1.5, -2.0, 0.5);
        assert_eq!(term.total(), 0.0);
    }

    #[test]
    fn add_sums_each_field_correctly() {
        let a = EnergyTerm::from_nonbonded(1.0, 2.0, 3.0);
        let b = EnergyTerm::from_nonbonded(4.0, 5.0, 6.0);
        let result = a + b;
        assert_eq!(result, EnergyTerm::from_nonbonded(5.0, 7.0, 9.0));
    }

    #[test]
    fn add_assign_accumulates_each_field_correctly() {
        let mut a = EnergyTerm::from_nonbonded(1.0, 2.0, 3.0);
        let b = EnergyTerm::from_nonbonded(4.0, 5.0, 6.0);
        a += b;
        assert_eq!(a, EnergyTerm::from_nonbonded(5.0, 7.0, 9.0));
    }

    #[test]
    fn default_initializes_all_fields_to_zero() {
        let term = EnergyTerm::default();
        assert_eq!(term.vdw, 0.0);
        assert_eq!(term.coulomb, 0.0);
        assert_eq!(term.hbond, 0.0);
    }

    #[test]
    fn add_with_negative_values() {
        let a = EnergyTerm::from_nonbonded(-1.0, 2.0, -3.0);
        let b = EnergyTerm::from_nonbonded(4.0, -5.0, 6.0);
        let result = a + b;
        assert_eq!(result, EnergyTerm::from_nonbonded(3.0, -3.0, 3.0));
    }

    #[test]
    fn add_assign_with_zero_does_not_change_values() {
        let mut a = EnergyTerm::from_nonbonded(1.0, 2.0, 3.0);
        let b = EnergyTerm::default();
        a += b;
        assert_eq!(a, EnergyTerm::from_nonbonded(1.0, 2.0, 3.0));
    }
}

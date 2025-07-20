use std::ops::{Add, AddAssign};

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct EnergyTerm {
    pub vdw: f64,
    pub coulomb: f64,
    pub hbond: f64,
}

impl EnergyTerm {
    pub fn new(vdw: f64, coulomb: f64, hbond: f64) -> Self {
        Self {
            vdw,
            coulomb,
            hbond,
        }
    }

    #[inline]
    pub fn total(&self) -> f64 {
        self.vdw + self.coulomb + self.hbond
    }
}

impl Add for EnergyTerm {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            vdw: self.vdw + rhs.vdw,
            coulomb: self.coulomb + rhs.coulomb,
            hbond: self.hbond + rhs.hbond,
        }
    }
}

impl AddAssign for EnergyTerm {
    fn add_assign(&mut self, rhs: Self) {
        self.vdw += rhs.vdw;
        self.coulomb += rhs.coulomb;
        self.hbond += rhs.hbond;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_creates_energy_term_with_specified_values() {
        let term = EnergyTerm::new(1.0, 2.0, 3.0);
        assert_eq!(term.vdw, 1.0);
        assert_eq!(term.coulomb, 2.0);
        assert_eq!(term.hbond, 3.0);
    }

    #[test]
    fn total_returns_sum_of_all_terms() {
        let term = EnergyTerm::new(1.5, -2.0, 0.5);
        assert_eq!(term.total(), 0.0);
    }

    #[test]
    fn add_sums_each_field_correctly() {
        let a = EnergyTerm::new(1.0, 2.0, 3.0);
        let b = EnergyTerm::new(4.0, 5.0, 6.0);
        let result = a + b;
        assert_eq!(result, EnergyTerm::new(5.0, 7.0, 9.0));
    }

    #[test]
    fn add_assign_accumulates_each_field_correctly() {
        let mut a = EnergyTerm::new(1.0, 2.0, 3.0);
        let b = EnergyTerm::new(4.0, 5.0, 6.0);
        a += b;
        assert_eq!(a, EnergyTerm::new(5.0, 7.0, 9.0));
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
        let a = EnergyTerm::new(-1.0, 2.0, -3.0);
        let b = EnergyTerm::new(4.0, -5.0, 6.0);
        let result = a + b;
        assert_eq!(result, EnergyTerm::new(3.0, -3.0, 3.0));
    }

    #[test]
    fn add_assign_with_zero_does_not_change_values() {
        let mut a = EnergyTerm::new(1.0, 2.0, 3.0);
        let b = EnergyTerm::default();
        a += b;
        assert_eq!(a, EnergyTerm::new(1.0, 2.0, 3.0));
    }
}

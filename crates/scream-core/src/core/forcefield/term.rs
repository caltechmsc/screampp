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

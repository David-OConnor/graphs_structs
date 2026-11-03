use crate::{
    Transformable,
    chem::{Forbidden, Term},
};

impl Transformable for SpectroState {
    type Context = SpectroCtx;
    type Step = SpectroStep;
    type Invalid = Forbidden;

    fn steps(&self, _ctx: &Self::Context) -> Vec<Self::Step> {
        todo!()
    }

    fn instances(_ctx: &Self::Context) -> Vec<Self> {
        todo!()
    }

    fn transform(&self, _step: &Self::Step, _ctx: &Self::Context) -> Result<Self, Self::Invalid> {
        todo!()
    }
}

pub struct SpectroCtx {
    pub terms: Vec<Term>,
    /// Separated-atom limits, cm⁻¹. E.g. H(1s) + H(1s).
    pub asymptotes: Vec<f64>,
    pub homonuclear: bool,
    /// Nuclear spin × 2. H: 1. ¹²C: 0.
    pub nuc_spin_x2: u8,
    pub j_max: u8,
}

/// Stationary state as quantum numbers. The graph is an energy-level diagram.
/// Bond length is derived: R_v from B_v. No oscillation; eigenstates don't move.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SpectroState {
    Bound {
        /// Index into `SpectroCtx::terms`.
        term: u8,
        v: u8,
        j: u8,
    },
    /// Continuum collapsed to one vertex per asymptote.
    Dissociated { asymptote: u8 },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpectroStep {
    /// Absorb or emit, by sign of ΔE. Electric dipole rules.
    Photon {
        term: u8,
        dv: i8,
        dj: i8,
    },
    /// ΔJ = 0, ±2. The only vib-rot route for H2.
    Raman {
        dv: i8,
        dj: i8,
    },
    Collision {
        dv: i8,
        dj: i8,
    },
    Dissociate {
        asymptote: u8,
    },
}

//! Similar to the mol and node graphs we use in classical molecule representations. May be too
//! high level for what we're going for? And/or may provide insights to those applications,
//! even if outside of this project's more ambitious scope.

use crate::{
    Transformable,
    chem::{AtomNode, BondEdge, BondOrbital, Forbidden, Hybrid, Spin},
};

/// Atoms are nodes, bonds are edges; electrons live on either. A Lewis structure with MO labels
/// and spin. Generalizes past diatomics. No coordinates: geometry follows from bond lengths and
/// hybridization.
///
/// C2 ground state, roughly: 1 lone pair per C, π⁴, σ⁰. H2: σ with 1 up, 1 down.
/// Isomorphic states must hash equal; keep a canonical atom order.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ChemGraph {
    pub atoms: Vec<AtomNode>,
    pub bonds: Vec<BondEdge>,
}

impl Transformable for ChemGraph {
    type Context = ChemGraphCtx;
    type Step = ChemGraphStep;
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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChemGraphStep {
    /// Lone electron into a bond orbital.
    Share {
        atom: u8,
        bond: u8,
        orbital: BondOrbital,
        spin: Spin,
    },
    Unshare {
        atom: u8,
        bond: u8,
        orbital: BondOrbital,
        spin: Spin,
    },
    /// Within a bond; e.g. σ → σ*.
    Excite {
        bond: u8,
        from: BondOrbital,
        to: BondOrbital,
        spin: Spin,
    },
    Stretch {
        bond: u8,
        dr: i8,
    },
    Vibrate {
        bond: u8,
        dv: i8,
    },
    Rehybridize {
        atom: u8,
        to: Hybrid,
    },
    Form {
        atoms: (u8, u8),
    },
    /// Electrons return to atoms. Homolytic or heterolytic, by count.
    Break {
        bond: u8,
    },
}

pub struct ChemGraphCtx {
    /// Å
    pub r_min: f32,
    pub dr: f32,
    pub n_r: u8,
    pub v_max: u8,
}

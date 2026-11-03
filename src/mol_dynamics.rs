//! Classical MD. This probably won't work?
//! todo: Discrete steps, or distance weightings?

// use mol_defs::molecules::common::MoleculeCommon;
use lin_alg::f32::Vec3;

use crate::Transformable;

// todo: import from mol_defs?
/// This is the part of molecule state which is static; not part of our transformation
#[derive(Clone)]
struct MolCtx {
    atoms: Vec<Vec3>,
    bonds: Vec<(usize, usize)>,
}

// todo: multi-mol System instead of the mol?

#[derive(Clone)]
struct MolConformation {
    torsions: Vec<f32>, // todo: Discrete torsions for graph purposes?
}

#[derive(Clone)]
struct MolStep {}

enum MolInvalid {}

impl Transformable for MolConformation {
    type Context = MolCtx;
    type Step = MolStep;
    type Invalid = MolInvalid;

    fn steps(&self, ctx: &Self::Context) -> Vec<Self::Step> {
        todo!()
    }

    fn instances(ctx: &Self::Context) -> Vec<Self> {
        todo!()
    }

    fn transform(&self, step: &Self::Step, ctx: &Self::Context) -> Result<Self, Self::Invalid>
    where
        Self: Sized,
    {
        todo!()
    }
}

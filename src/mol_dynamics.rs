//! Classical MD. This probably won't work?
//! todo: Discrete steps, or distance weightings?

use mol_defs::molecules::common::MoleculeCommon;
use lin_alg::f32::Vec3;
use crate::Transformable;

// todo: import from mol_defs?
#[derive(Clone)]
pub struct Mol {
    atoms: Vec<Vec3>,
    bonds: Vec<(usize, usize)>
}

impl Transformable for Mol {
    type Context = ();
    type Step = ();
    type Invalid = ();

    fn steps(&self, ctx: &Self::Context) -> Vec<Self::Step> {
        todo!()
    }

    fn instances(ctx: &Self::Context) -> Vec<Self> {
        todo!()
    }

    fn transform(&self, step: &Self::Step, ctx: &Self::Context) -> Result<Self, Self::Invalid>
    where
        Self: Sized
    {
        todo!()
    }
}
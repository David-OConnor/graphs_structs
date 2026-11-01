//! Long shot

//! Classical MD. This probably won't work?
//! todo: Discrete steps, or distance weightings?

use crate::Transformable;

// todo: import from mol_defs?
#[derive(Clone)]
pub struct Atom {
    num_protons: usize,
    num_neutrons: usize,
    num_electrons: usize,
}

// todo: import from mol_defs?
#[derive(Clone)]
pub struct DiatomicMol {

}

impl Transformable for Atom {
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
//! Long shot

//! Classical MD. This probably won't work?
//! todo: Discrete steps, or distance weightings?

use crate::Transformable;

#[derive(Clone, Copy, PartialEq)]
enum EnergyLevel {
    Ground,
    _1,
    _2,
    _3,
}



// todo: import from mol_defs?
#[derive(Clone)]
pub struct Atom {
    num_protons: usize,
    num_neutrons: usize,
    num_electrons: usize,
    /// Len = num_electrons
    elec_energy_levels: Vec<EnergyLevel>,
}

// todo: import from mol_defs?
#[derive(Clone)]
pub struct DiatomicMol {
    atom_0: Atom,
    atm_1: Atom,
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
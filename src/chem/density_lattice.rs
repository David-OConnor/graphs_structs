use crate::{
    Transformable,
    chem::{Element, Forbidden, Spin},
};

/// Electron density as conserved quanta per voxel: a lattice gas of charge. Quanta are charge,
/// not density; voxel volume grows with ρ. Bond length: `nuc_z[1] - nuc_z[0]`.
///
/// Huge state space. Build `ConnectedOnly`, with steps pruned to an energy window.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct DensityLattice {
    /// Index: z · n_rho + ρ.
    pub up: Vec<u8>,
    pub down: Vec<u8>,
    pub nuc_z: [u16; 2],
}

impl Transformable for DensityLattice {
    type Context = DensityCtx;
    type Step = DensityStep;
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

pub struct DensityCtx {
    pub elements: [Element; 2],
    /// Bond axis is z. Axial symmetry reduces 3D to (z, ρ).
    pub n_z: u16,
    pub n_rho: u16,
    /// Å
    pub dz: f32,
    pub d_rho: f32,
    /// Resolution.
    pub quanta_per_e: u8,
    /// Max quanta per voxel per spin. Crude local Pauli.
    pub cap: u8,
    // todo: Energy functional. Orbital-free (TF alone won't bind; Teller), or learned.
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DensityStep {
    /// One quantum to an adjacent voxel.
    Flow {
        spin: Spin,
        from: u16,
        to: u16,
    },
    /// Needs spin-orbit coupling, or a field.
    Flip {
        from: Spin,
        voxel: u16,
    },
    MoveNucleus {
        which: u8,
        dz: i8,
    },
}

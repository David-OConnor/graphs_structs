use crate::{
    Transformable, chem,
    chem::{BOHR_PER_ANGSTROM, FS_PER_AU_TIME, Forbidden},
};

pub struct PhaseCellCtx {
    /// Å, at r = 0.
    pub r_min: f32,
    /// Å per r step.
    pub dr: f32,
    pub n_r: u16,
    pub p_max: i16,
    /// Per surface, per r. Lattice units; dt and μ folded in.
    pub force: Vec<Vec<i16>>,
    /// Per surface, per r. Hartree.
    pub potential: Vec<Vec<f64>>,
    /// Reduced mass, electron masses.
    pub mu: f64,
    /// fs per tick.
    pub dt: f64,
    /// Tried from every state. [Tick]: `ConnectedOnly` gives one trajectory.
    pub candidate_steps: Vec<PhaseCellStep>,
}

impl PhaseCellCtx {
    /// One surface per energy curve: (Å, Hartree), ascending; e.g. from `ground_curve`. Spans the
    /// r range they share. Force is linear between the curves' segment slopes.
    /// A unit of p moves 1 cell per tick. At most 65535 cells.
    pub fn from_curves(curves: &[Vec<(f32, f64)>], mu: f64, dr: f32, dt: f64, p_max: i16) -> Self {
        let r_min = curves.iter().map(|c| c[0].0).fold(f32::MIN, f32::max);
        let r_max = curves
            .iter()
            .map(|c| c[c.len() - 1].0)
            .fold(f32::MAX, f32::min);
        let n_r = u16::try_from(((r_max - r_min) / dr) as u32 + 1).expect("Over 65535 cells");
        let lattice: Vec<f64> = (0..n_r)
            .map(|k| f64::from(r_min) + f64::from(k) * f64::from(dr))
            .collect();

        // Hartree/Å → cells/tick, per tick.
        let dt_au = dt / FS_PER_AU_TIME;
        let to_lattice = dt_au.powi(2) / (mu * f64::from(dr) * BOHR_PER_ANGSTROM.powi(2));

        let mut force = Vec::new();
        let mut potential = Vec::new();
        for curve in curves {
            let curve: Vec<(f64, f64)> = curve.iter().map(|&(r, e)| (f64::from(r), e)).collect();
            let slopes: Vec<(f64, f64)> = curve
                .windows(2)
                .map(|w| {
                    let ((r0, e0), (r1, e1)) = (w[0], w[1]);

                    ((r0 + r1) / 2., -(e1 - e0) / (r1 - r0))
                })
                .collect();
            let f: Vec<f64> = lattice
                .iter()
                .map(|&r| chem::interpolate(&slopes, r))
                .collect();

            // Integrate −F, so V and F agree.
            let mut v = vec![chem::interpolate(&curve, lattice[0])];
            for w in f.windows(2) {
                v.push(v[v.len() - 1] - (w[0] + w[1]) / 2. * f64::from(dr));
            }

            force.push(f.iter().map(|x| (x * to_lattice).round() as i16).collect());
            potential.push(v);
        }

        Self {
            r_min,
            dr,
            n_r,
            p_max,
            force,
            potential,
            mu,
            dt,
            candidate_steps: vec![PhaseCellStep::Tick],
        }
    }

    /// At rest, at the cell nearest `r` (Å).
    pub fn at(&self, surface: u8, r: f32) -> PhaseCell {
        let r = ((r - self.r_min) / self.dr).round();

        PhaseCell {
            surface,
            r: r.clamp(0., f32::from(self.n_r - 1)) as u16,
            p: 0,
        }
    }

    /// `p`, if within ±`p_max`.
    fn momentum(&self, p: i32) -> Result<i16, Forbidden> {
        i16::try_from(p)
            .ok()
            .filter(|p| (-self.p_max..=self.p_max).contains(p))
            .ok_or(Forbidden::Unbound)
    }
}

/// Bond coordinate and momentum on an integer lattice, on one electronic surface.
///
/// `Tick` is integer leapfrog: p += F(r); r += p. Each half is a shear, so it's exactly
/// invertible, and bounded orbits are cycles. A cycle is the lattice state's exact recurrence; it
/// usually spans many vibrations, seen as r oscillating along it. Orbits leaving the grid
/// dissociate.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PhaseCell {
    /// Index into `PhaseCellCtx::force`. E.g. a `Term` or `Determinant` index.
    pub surface: u8,
    pub r: u16,
    pub p: i16,
}

impl PhaseCell {
    /// Å
    pub fn bond_length(self, ctx: &PhaseCellCtx) -> f32 {
        ctx.r_min + f32::from(self.r) * ctx.dr
    }

    /// Kinetic + potential, Hartree. p lags r by half a tick; adding F/2 aligns them.
    pub fn energy(self, ctx: &PhaseCellCtx) -> f64 {
        let (s, r) = (usize::from(self.surface), usize::from(self.r));
        let p = f64::from(self.p) + f64::from(ctx.force[s][r]) / 2.;
        let v = p * f64::from(ctx.dr) * BOHR_PER_ANGSTROM / (ctx.dt / FS_PER_AU_TIME);

        0.5 * ctx.mu * v * v + ctx.potential[s][r]
    }
}

impl Transformable for PhaseCell {
    type Context = PhaseCellCtx;
    type Step = PhaseCellStep;
    type Invalid = Forbidden;

    fn steps(&self, ctx: &Self::Context) -> Vec<Self::Step> {
        ctx.candidate_steps
            .iter()
            .copied()
            .filter(|step| self.transform(step, ctx).is_ok())
            .collect()
    }

    /// Surfaces × n_r × (2·p_max + 1): often too many to build.
    fn instances(ctx: &Self::Context) -> Vec<Self> {
        let mut instances = Vec::new();
        for surface in (0..=u8::MAX).take(ctx.force.len()) {
            for r in 0..ctx.n_r {
                for p in -ctx.p_max..=ctx.p_max {
                    instances.push(Self { surface, r, p });
                }
            }
        }

        instances
    }

    fn transform(&self, step: &Self::Step, ctx: &Self::Context) -> Result<Self, Self::Invalid> {
        let Some(force) = ctx.force.get(usize::from(self.surface)) else {
            return Err(Forbidden::Malformed);
        };
        if self.r >= ctx.n_r || !(-ctx.p_max..=ctx.p_max).contains(&self.p) {
            return Err(Forbidden::Malformed);
        }

        let mut next = *self;
        match *step {
            PhaseCellStep::Tick => {
                next.p = ctx.momentum(i32::from(self.p) + i32::from(force[usize::from(self.r)]))?;
                next.r = u16::try_from(i32::from(self.r) + i32::from(next.p))
                    .ok()
                    .filter(|&r| r < ctx.n_r)
                    .ok_or(Forbidden::Unbound)?;
            }
            PhaseCellStep::Hop { to } => {
                if to == self.surface {
                    return Err(Forbidden::Malformed);
                }
                if usize::from(to) >= ctx.force.len() {
                    return Err(Forbidden::Unbound);
                }

                next.surface = to;
            }
            PhaseCellStep::Kick(dp) => next.p = ctx.momentum(i32::from(self.p) + i32::from(dp))?,
        }

        Ok(next)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhaseCellStep {
    Tick,
    /// Electronic transition at fixed r: photon, or surface hop near a crossing.
    Hop {
        to: u8,
    },
    /// Momentum transfer, e.g. a collision.
    Kick(i16),
}

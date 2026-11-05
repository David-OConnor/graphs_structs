//! Electronic superpositions at fixed nuclei, in STO-3G. Atomic units except time (fs).
//! Matter waves obey Schrödinger evolution, not propagation at c. Density beats at energy
//! differences; an energy eigenstate has stationary density. Pauli exclusion and electron
//! repulsion enter through Slater determinants. Minimal-basis FCI, not a continuum model.

use std::f64::consts::TAU;

use crate::{
    Transformable,
    chem::{
        self, Element, FS_PER_AU_TIME, Forbidden, SpinCount,
        determinant::{DeterminantCtx, SlaterDeterminant},
        sto3g,
    },
};

/// |c_k|² = pop[k] / quanta. First occupied phase and all empty phases are zero.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct WavePackets {
    pub pop: Vec<u8>,
    pub phase: Vec<u16>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WavePacketsStep {
    Tick,
    /// Move one population quantum; state preparation, not free evolution.
    Transfer {
        from: usize,
        to: usize,
    },
    /// Select an occupied energy eigenstate; outcome probability is pop[to] / quanta.
    Collapse {
        to: usize,
    },
}

pub struct WavePacketsCtx {
    pub(crate) energies: Vec<f64>,
    pub(crate) vectors: Vec<Vec<f64>>,
    pub(crate) determinants: Vec<SlaterDeterminant>,
    aos: Vec<sto3g::Ao>,
    pub(crate) coefficients: Vec<Vec<f64>>,
    pub(crate) quanta: u8,
    pub(crate) phase_bins: u16,
    pub(crate) omega: Vec<u16>,
    pub(crate) dt: f64,
    /// Tick only by default. Other steps prepare or measure states.
    pub candidate_steps: Vec<WavePacketsStep>,
}

impl WavePacketsCtx {
    /// Lowest `levels` in a fixed (N↑, N↓) sector; nuclear positions in bohr, dt in fs.
    /// Dense FCI: at most 10 spatial orbitals and 512 determinants. H, C and H2 fit;
    /// low-spin C2 exceeds this limit. STO-3G gives isolated H only one spatial mode.
    /// `phase_bins` ≥ 2; phase advances round to this grid, so small gaps may freeze.
    pub fn new(
        atoms: &[(Element, [f64; 3])],
        electrons: SpinCount,
        levels: usize,
        quanta: u8,
        phase_bins: u16,
        dt: f64,
    ) -> Result<Self, Forbidden> {
        fn masks(n: usize, count: u8) -> Vec<u64> {
            (0..1u64 << n)
                .filter(|bits| bits.count_ones() == u32::from(count))
                .map(|bits| (0..n).map(|k| ((bits >> k) & 1) << (2 * k)).sum())
                .collect()
        }

        if atoms.is_empty()
            || levels == 0
            || quanta == 0
            || phase_bins < 2
            || !dt.is_finite()
            || dt <= 0.
            || atoms.iter().any(|(_, r)| r.iter().any(|x| !x.is_finite()))
        {
            return Err(Forbidden::Malformed);
        }

        let n: usize = atoms
            .iter()
            .map(|(e, _)| if *e == Element::H { 1 } else { 5 })
            .sum();
        if n > 10 || usize::from(electrons.up.max(electrons.down)) > n {
            return Err(Forbidden::Unbound);
        }

        for (i, (_, a)) in atoms.iter().enumerate() {
            for (_, b) in &atoms[i + 1..] {
                let r2: f64 = (0..3).map(|k| (a[k] - b[k]).powi(2)).sum();
                if !r2.is_finite() || r2 < 1e-12 {
                    return Err(Forbidden::Malformed);
                }
            }
        }

        let up = masks(n, electrons.up);
        let down = masks(n, electrons.down);
        let size = up.len() * down.len();
        if size > 512 || levels > size {
            return Err(Forbidden::Unbound);
        }

        let determinants: Vec<_> = up
            .iter()
            .flat_map(|&u| {
                down.iter().map(move |&d| SlaterDeterminant {
                    occ: u | (d << 1),
                    r: 0,
                })
            })
            .collect();
        let (integrals, overlap) = sto3g::molecule(atoms);
        if chem::eigh(overlap.clone())
            .iter()
            .any(|(e, _)| !e.is_finite() || *e < 1e-10)
        {
            return Err(Forbidden::Malformed);
        }

        let coefficients = chem::inverse_sqrt(&overlap);
        // Only integrals are read by SlaterDeterminant::hamiltonian.
        let det_ctx = DeterminantCtx {
            orbitals: Vec::new(),
            z_total: atoms.iter().map(|(e, _)| e.z()).sum(),
            charges: 0..=0,
            r_grid: vec![0.],
            integrals: vec![integrals.transformed(&coefficients)],
            spin_orbit: false,
            stretch: false,
        };
        let matrix: Vec<Vec<f64>> = determinants
            .iter()
            .map(|a| {
                determinants
                    .iter()
                    .map(|b| a.hamiltonian(*b, &det_ctx).unwrap())
                    .collect()
            })
            .collect();
        if matrix.iter().flatten().any(|x| !x.is_finite()) {
            return Err(Forbidden::Malformed);
        }

        let (energies, vectors): (Vec<_>, Vec<_>) =
            chem::eigh(matrix).into_iter().take(levels).unzip();
        let mut omega = Vec::with_capacity(levels);
        for &energy in &energies {
            let turns = (energy - energies[0]) * (dt / FS_PER_AU_TIME) / TAU;
            if !turns.is_finite() {
                return Err(Forbidden::Malformed);
            }

            omega.push(
                ((turns.rem_euclid(1.) * f64::from(phase_bins)).round() as u32
                    % u32::from(phase_bins)) as u16,
            );
        }

        Ok(Self {
            energies,
            vectors,
            determinants,
            aos: sto3g::basis(atoms),
            coefficients,
            quanta,
            phase_bins,
            omega,
            dt,
            candidate_steps: vec![WavePacketsStep::Tick],
        })
    }

    /// Hartree, ascending; includes nuclear repulsion.
    pub fn energies(&self) -> &[f64] {
        &self.energies
    }

    /// fs per tick.
    pub fn dt(&self) -> f64 {
        self.dt
    }
}

impl WavePackets {
    pub fn pure(level: usize, ctx: &WavePacketsCtx) -> Result<Self, Forbidden> {
        if level >= ctx.energies.len() {
            return Err(Forbidden::Unbound);
        }

        let mut pop = vec![0; ctx.energies.len()];
        pop[level] = ctx.quanta;

        Ok(Self {
            phase: vec![0; pop.len()],
            pop,
        })
    }

    /// Canonicalize a normalized population split and phases in grid units.
    pub fn new(pop: Vec<u8>, phase: Vec<u16>, ctx: &WavePacketsCtx) -> Result<Self, Forbidden> {
        if pop.len() != ctx.energies.len()
            || phase.len() != pop.len()
            || pop.iter().map(|&p| usize::from(p)).sum::<usize>() != usize::from(ctx.quanta)
        {
            return Err(Forbidden::Malformed);
        }

        Ok(Self { pop, phase }.canonical(ctx))
    }

    fn canonical(mut self, ctx: &WavePacketsCtx) -> Self {
        let reference = self
            .pop
            .iter()
            .position(|&p| p > 0)
            .map_or(0, |i| self.phase[i]);
        for (p, phase) in self.pop.iter().zip(&mut self.phase) {
            *phase = if *p == 0 {
                0
            } else {
                (i32::from(*phase) - i32::from(reference)).rem_euclid(i32::from(ctx.phase_bins))
                    as u16
            };
        }

        self
    }

    fn validate(&self, ctx: &WavePacketsCtx) -> Result<(), Forbidden> {
        if Self::new(self.pop.clone(), self.phase.clone(), ctx)? == *self {
            return Ok(());
        }

        Err(Forbidden::Malformed)
    }

    /// Expected total energy, Hartree.
    pub fn energy(&self, ctx: &WavePacketsCtx) -> Result<f64, Forbidden> {
        self.validate(ctx)?;

        Ok(self
            .pop
            .iter()
            .zip(&ctx.energies)
            .map(|(&p, e)| f64::from(p) * e / f64::from(ctx.quanta))
            .sum())
    }

    /// Real, spin-summed one-particle density matrix in the orthonormal AO basis.
    /// Its trace is the electron count. Imaginary coherences do not enter spatial density.
    pub fn density_matrix(&self, ctx: &WavePacketsCtx) -> Result<Vec<Vec<f64>>, Forbidden> {
        self.validate(ctx)?;

        let mut amplitudes = vec![[0.; 2]; ctx.determinants.len()];
        for (k, vector) in ctx.vectors.iter().enumerate() {
            let magnitude = (f64::from(self.pop[k]) / f64::from(ctx.quanta)).sqrt();
            let angle = TAU * f64::from(self.phase[k]) / f64::from(ctx.phase_bins);
            let (sin, cos) = angle.sin_cos();
            for (a, v) in amplitudes.iter_mut().zip(vector) {
                a[0] += magnitude * cos * v;
                a[1] += magnitude * sin * v;
            }
        }

        let mut density = vec![vec![0.; ctx.aos.len()]; ctx.aos.len()];
        for (i, det) in ctx.determinants.iter().enumerate() {
            let weight = amplitudes[i].iter().map(|a| a * a).sum::<f64>();
            for p in chem::bits(det.occ) {
                density[usize::from(p / 2)][usize::from(p / 2)] += weight;
            }

            for (j, other) in ctx.determinants.iter().enumerate().skip(i + 1) {
                let removed = det.occ & !other.occ;
                let added = other.occ & !det.occ;
                if removed.count_ones() != 1 || added.count_ones() != 1 {
                    continue;
                }

                let p = removed.trailing_zeros() as u8;
                let q = added.trailing_zeros() as u8;
                if p % 2 != q % 2 {
                    continue;
                }

                let coherence = chem::sign(det.occ, p, q)
                    * (amplitudes[i][0] * amplitudes[j][0] + amplitudes[i][1] * amplitudes[j][1]);
                density[usize::from(p / 2)][usize::from(q / 2)] += coherence;
                density[usize::from(q / 2)][usize::from(p / 2)] += coherence;
            }
        }

        Ok(density)
    }

    /// Electrons / bohr³ at `r` (bohr), including interference between energy states.
    pub fn density(&self, r: [f64; 3], ctx: &WavePacketsCtx) -> Result<f64, Forbidden> {
        if r.iter().any(|x| !x.is_finite()) {
            return Err(Forbidden::Malformed);
        }

        let density = self.density_matrix(ctx)?;
        let raw: Vec<_> = ctx.aos.iter().map(|ao| ao.value(r)).collect();
        let orbitals: Vec<f64> = (0..raw.len())
            .map(|p| {
                raw.iter()
                    .enumerate()
                    .map(|(a, v)| v * ctx.coefficients[a][p])
                    .sum()
            })
            .collect();
        let value: f64 = chem::indices::<2>(raw.len())
            .map(|[p, q]| density[p][q] * orbitals[p] * orbitals[q])
            .sum();

        Ok(value.max(0.))
    }
}

impl Transformable for WavePackets {
    type Context = WavePacketsCtx;
    type Step = WavePacketsStep;
    type Invalid = Forbidden;

    fn steps(&self, ctx: &Self::Context) -> Vec<Self::Step> {
        ctx.candidate_steps
            .iter()
            .copied()
            .filter(|step| self.transform(step, ctx).is_ok())
            .collect()
    }

    /// All normalized splits and relative phases; exponential. Prefer ConnectedOnly.
    fn instances(ctx: &Self::Context) -> Vec<Self> {
        fn enumerate(
            state: &mut WavePackets,
            k: usize,
            remaining: u8,
            occupied: bool,
            ctx: &WavePacketsCtx,
            result: &mut Vec<WavePackets>,
        ) {
            if k == state.pop.len() {
                if remaining == 0 {
                    result.push(state.clone());
                }

                return;
            }

            let minimum = if k + 1 == state.pop.len() {
                remaining
            } else {
                0
            };
            for pop in minimum..=remaining {
                state.pop[k] = pop;
                let phases = if occupied && pop > 0 {
                    ctx.phase_bins
                } else {
                    1
                };
                for phase in 0..phases {
                    state.phase[k] = phase;
                    enumerate(
                        state,
                        k + 1,
                        remaining - pop,
                        occupied || pop > 0,
                        ctx,
                        result,
                    );
                }
            }
        }

        let mut result = Vec::new();
        let mut state = Self {
            pop: vec![0; ctx.energies.len()],
            phase: vec![0; ctx.energies.len()],
        };
        enumerate(&mut state, 0, ctx.quanta, false, ctx, &mut result);

        result
    }

    fn transform(&self, step: &Self::Step, ctx: &Self::Context) -> Result<Self, Self::Invalid> {
        self.validate(ctx)?;

        let mut next = self.clone();
        match *step {
            WavePacketsStep::Tick => {
                for (phase, &omega) in next.phase.iter_mut().zip(&ctx.omega) {
                    *phase = (i32::from(*phase) - i32::from(omega))
                        .rem_euclid(i32::from(ctx.phase_bins)) as u16;
                }
            }
            WavePacketsStep::Transfer { from, to } => {
                if from >= self.pop.len() || to >= self.pop.len() {
                    return Err(Forbidden::Unbound);
                }
                if from == to {
                    return Err(Forbidden::Malformed);
                }
                if self.pop[from] == 0 {
                    return Err(Forbidden::Empty);
                }

                if self.pop[to] == 0 {
                    next.phase[to] = self.phase[from];
                }
                next.pop[from] -= 1;
                next.pop[to] += 1;
            }
            WavePacketsStep::Collapse { to } => {
                if to >= self.pop.len() {
                    return Err(Forbidden::Unbound);
                }
                if self.pop[to] == 0 {
                    return Err(Forbidden::Empty);
                }

                return Self::pure(to, ctx);
            }
        }

        Ok(next.canonical(ctx))
    }
}

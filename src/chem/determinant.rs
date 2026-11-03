//! Slater determinants.

use DeterminantStep::*;

use crate::{
    Transformable, chem,
    chem::{Center, Element, Forbidden, Lambda, Orbital, OrbitalKind, Parity, UP_BITS, sto3g},
};

/// Indices are spin-orbitals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DeterminantStep {
    Excite {
        from: u8,
        to: u8,
    },
    DoubleExcite {
        from: [u8; 2],
        to: [u8; 2],
    },
    Ionize {
        from: u8,
    },
    Attach {
        to: u8,
    },
    /// ±1 bin.
    Stretch(i8),
}

/// One Slater determinant, and a bond-length bin. Charge, Ms, bond order are derived.
///
/// Spin-conserving `Excite` and `DoubleExcite` edges cover every nonzero CI Hamiltonian element
/// (Slater-Condon). `hamiltonian` gives vertex and edge weights: the graph *is* H.
///
/// H2, minimal basis: 4 spin-orbitals, 6 dets. With `Center::Atom`, these are the covalent and
/// ionic structures (Hubbard dimer). C2 valence: 16 spin-orbitals, 8 e-, 12870 dets; 4900 with
/// Ms = 0.
/// Dets aren't spin eigenstates; CSFs are (GUGA's Shavitt graph enumerates those).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SlaterDeterminant {
    /// Bit 2k: orbital k, up. Bit 2k + 1: down.
    pub occ: u64,
    /// Index into `DeterminantCtx::r_grid`.
    pub r: u8,
}

impl SlaterDeterminant {
    pub fn n_elec(self) -> u8 {
        self.occ.count_ones() as u8
    }

    pub fn charge(self, ctx: &DeterminantCtx) -> i8 {
        ctx.z_total as i8 - self.n_elec() as i8
    }

    /// 2·Ms
    pub fn ms_x2(self) -> i8 {
        (self.occ & UP_BITS).count_ones() as i8 - (self.occ & !UP_BITS).count_ones() as i8
    }

    /// (Bonding − antibonding electrons) / 2. MO basis only.
    pub fn bond_order(self, ctx: &DeterminantCtx) -> f32 {
        let mut net = 0;
        for (k, orbital) in ctx.orbitals.iter().enumerate() {
            let n = ((self.occ >> (2 * k)) & 0b11).count_ones() as i8;
            match orbital.kind {
                OrbitalKind::Bonding => net += n,
                OrbitalKind::Antibonding => net -= n,
                OrbitalKind::Core | OrbitalKind::Nonbonding => (),
            }
        }

        f32::from(net) / 2.
    }

    /// ⟨self|H|other⟩, Hartree, by Slater-Condon. Vertex weight if `other == self`, else edge
    /// weight. Zero across r bins or electron counts. None without integrals.
    pub fn hamiltonian(self, other: Self, ctx: &DeterminantCtx) -> Option<f64> {
        let ints = ctx.integrals.get(usize::from(self.r))?;

        if self.r != other.r || self.n_elec() != other.n_elec() {
            return Some(0.);
        }

        let occupied = chem::bits(self.occ);
        let removed = chem::bits(self.occ & !other.occ);
        let added = chem::bits(other.occ & !self.occ);

        let element = match (&removed[..], &added[..]) {
            ([], []) => {
                let one: f64 = occupied.iter().map(|&i| ints.h_so(i, i)).sum();
                let two: f64 = occupied
                    .iter()
                    .flat_map(|&i| occupied.iter().map(move |&j| ints.anti(i, j, i, j)))
                    .sum();

                ints.e_nuc + one + two / 2.
            }
            (&[i], &[a]) => {
                let two: f64 = occupied.iter().map(|&j| ints.anti(i, j, a, j)).sum();

                chem::sign(self.occ, i, a) * (ints.h_so(i, a) + two)
            }
            (&[i, j], &[a, b]) => {
                let after_ia = self.occ ^ (1u64 << i) ^ (1u64 << a);

                chem::sign(self.occ, i, a) * chem::sign(after_ia, j, b) * ints.anti(i, j, a, b)
            }
            _ => 0.,
        };

        Some(element)
    }

    /// Empties `from` and fills `to`; they must be full and empty.
    fn moved(self, from: &[u8], to: &[u8], ctx: &DeterminantCtx) -> Result<u64, Forbidden> {
        let from = ctx.mask(from)?;
        let to = ctx.mask(to)?;

        if (self.occ & from) != from {
            return Err(Forbidden::Empty);
        }
        if (self.occ & to) != 0 {
            return Err(Forbidden::Pauli);
        }

        Ok((self.occ & !from) | to)
    }
}

impl Transformable for SlaterDeterminant {
    type Context = DeterminantCtx;
    type Step = DeterminantStep;
    type Invalid = Forbidden;

    fn steps(&self, ctx: &Self::Context) -> Vec<Self::Step> {
        let pairs = |items: &[u8]| {
            let mut pairs = Vec::new();
            for (i, &a) in items.iter().enumerate() {
                for &b in &items[i + 1..] {
                    pairs.push([a, b]);
                }
            }

            pairs
        };

        let (occupied, empty): (Vec<u8>, Vec<u8>) =
            (0..ctx.n_spin_orbitals()).partition(|&i| (self.occ & 1u64 << i) != 0);

        let mut steps = if ctx.stretch {
            vec![Stretch(-1), Stretch(1)]
        } else {
            Vec::new()
        };
        for &from in &occupied {
            steps.push(Ionize { from });
            steps.extend(empty.iter().map(|&to| Excite { from, to }));
        }
        steps.extend(empty.iter().map(|&to| Attach { to }));

        let empty_pairs = pairs(&empty);
        for from in pairs(&occupied) {
            steps.extend(empty_pairs.iter().map(|&to| DoubleExcite { from, to }));
        }

        steps
            .into_iter()
            .filter(|step| self.transform(step, ctx).is_ok())
            .collect()
    }

    fn instances(ctx: &Self::Context) -> Vec<Self> {
        let mut instances = Vec::new();
        for r in (0..=u8::MAX).take(ctx.r_grid.len()) {
            for occ in 0..=ctx.all_bits() {
                let det = Self { occ, r };
                if ctx.charges.contains(&det.charge(ctx)) {
                    instances.push(det);
                }
            }
        }

        instances
    }

    fn transform(&self, step: &Self::Step, ctx: &Self::Context) -> Result<Self, Self::Invalid> {
        if (self.occ & !ctx.all_bits()) != 0
            || usize::from(self.r) >= ctx.r_grid.len()
            || !ctx.charges.contains(&self.charge(ctx))
        {
            return Err(Forbidden::Malformed);
        }

        let mut next = *self;
        match *step {
            DeterminantStep::Excite { from, to } => {
                next.occ = self.moved(&[from], &[to], ctx)?;

                if next.ms_x2() != self.ms_x2() && !ctx.spin_orbit {
                    return Err(Forbidden::Spin);
                }
            }
            DeterminantStep::DoubleExcite { from, to } => {
                next.occ = self.moved(&from, &to, ctx)?;

                // Spin-orbit is 1-electron; it can't change Ms in a double.
                if next.ms_x2() != self.ms_x2() {
                    return Err(Forbidden::Spin);
                }
            }
            DeterminantStep::Ionize { from } => next.occ = self.moved(&[from], &[], ctx)?,
            DeterminantStep::Attach { to } => next.occ = self.moved(&[], &[to], ctx)?,
            DeterminantStep::Stretch(dr) => {
                if !matches!(dr, -1 | 1) {
                    return Err(Forbidden::NotAdjacent);
                }

                next.r = self
                    .r
                    .checked_add_signed(dr)
                    .filter(|&r| usize::from(r) < ctx.r_grid.len())
                    .ok_or(Forbidden::Unbound)?;
            }
        }

        if !ctx.charges.contains(&next.charge(ctx)) {
            return Err(Forbidden::Unbound);
        }

        Ok(next)
    }
}

pub struct DeterminantCtx {
    /// At most 32; `occ` is a u64.
    pub orbitals: Vec<Orbital>,
    /// Sum of nuclear charges, less frozen core electrons.
    pub z_total: u8,
    /// Allowed molecular charges. Bounds `Ionize` and `Attach`.
    pub charges: std::ops::RangeInclusive<i8>,
    /// Å, per bond-length bin. At most 256; `r` is a u8.
    pub r_grid: Vec<f32>,
    /// Per r bin, over `orbitals`. Empty: no energies.
    pub integrals: Vec<Integrals>,
    /// Allows Ms-changing `Excite`.
    pub spin_orbit: bool,
    /// Generate `Stretch` steps. Off: the graph stays on one bin; `lowest_curve` reweights it.
    pub stretch: bool,
}

impl DeterminantCtx {
    /// H2, STO-3G, as MOs. H2+, H2, H2-.
    pub fn h2_mo(r_grid: Vec<f32>) -> Self {
        let sigma = |label, parity, kind| Orbital {
            label,
            lambda: Lambda::Sigma,
            parity: Some(parity),
            center: Center::Delocalized,
            kind,
        };

        // σg, σu ∝ 1sA ± 1sB.
        let coeffs = |s: f64| {
            let g = 1. / (2. * (1. + s)).sqrt();
            let u = 1. / (2. * (1. - s)).sqrt();

            vec![vec![g, u], vec![g, -u]]
        };

        let orbitals = vec![
            sigma("σg1s", Parity::Gerade, OrbitalKind::Bonding),
            sigma("σu*1s", Parity::Ungerade, OrbitalKind::Antibonding),
        ];

        Self::h2(orbitals, r_grid, coeffs)
    }

    /// H2, STO-3G, as Löwdin-orthogonalized 1s orbitals: a site basis. Neutral dets: covalent H·H
    /// and ionic H+H-.
    pub fn h2_ao(r_grid: Vec<f32>) -> Self {
        let s1 = |label, atom| Orbital {
            label,
            lambda: Lambda::Sigma,
            parity: None,
            center: Center::Atom(atom),
            kind: OrbitalKind::Nonbonding,
        };

        // S^(-1/2)
        let coeffs = |s: f64| {
            let g = 1. / (1. + s).sqrt();
            let u = 1. / (1. - s).sqrt();

            vec![
                vec![(g + u) / 2., (g - u) / 2.],
                vec![(g - u) / 2., (g + u) / 2.],
            ]
        };

        Self::h2(vec![s1("1sA", 0), s1("1sB", 1)], r_grid, coeffs)
    }

    /// `coeffs`: AO overlap → AO coefficients, [μ][p].
    fn h2(orbitals: Vec<Orbital>, r_grid: Vec<f32>, coeffs: impl Fn(f64) -> Vec<Vec<f64>>) -> Self {
        let integrals = r_grid
            .iter()
            .map(|&r| {
                let (ao, s) = sto3g::molecule(&chem::diatomic(Element::H, r));
                ao.transformed(&coeffs(s[0][1]))
            })
            .collect();

        Self {
            orbitals,
            z_total: 2,
            charges: -1..=1,
            r_grid,
            integrals,
            spin_orbit: false,
            stretch: true,
        }
    }

    /// C2, STO-3G, 1s cores frozen. Neutral; no `Stretch`. Orbitals: symmetric and antisymmetric
    /// combinations of Löwdin-orthogonalized AOs. Not SCF orbitals, but FCI doesn't care.
    pub fn c2(r_grid: Vec<f32>) -> Self {
        use std::f64::consts::FRAC_1_SQRT_2;

        use crate::{
            chem,
            chem::{
                Lambda::{Pi, Sigma},
                OrbitalKind::{Antibonding, Bonding},
                Parity::{Gerade, Ungerade},
            },
        };

        let orbital = |label, lambda, parity, kind| Orbital {
            label,
            lambda,
            parity: Some(parity),
            center: Center::Delocalized,
            kind,
        };
        let orbitals = vec![
            orbital("σg2s", Sigma, Gerade, Bonding),
            orbital("σu*2s", Sigma, Ungerade, Antibonding),
            orbital("πu2px", Pi, Ungerade, Bonding),
            orbital("πu2py", Pi, Ungerade, Bonding),
            orbital("σg2p", Sigma, Gerade, Bonding),
            orbital("πg*2px", Pi, Gerade, Antibonding),
            orbital("πg*2py", Pi, Gerade, Antibonding),
            orbital("σu*2p", Sigma, Ungerade, Antibonding),
        ];

        // (AO per atom: 1s, 2s, 2px, 2py, 2pz; sign of atom B's). Core σg1s, σu1s, then
        // `orbitals`. pz points +z on both atoms, so σ bonding is A − B.
        let combos = [
            (0, 1.),
            (0, -1.),
            (1, 1.),
            (1, -1.),
            (2, 1.),
            (3, 1.),
            (4, -1.),
            (2, -1.),
            (3, -1.),
            (4, 1.),
        ];
        let mut u = vec![vec![0.; 10]; 10];
        for (p, &(k, sign)) in combos.iter().enumerate() {
            u[k][p] = FRAC_1_SQRT_2;
            u[k + 5][p] = sign * FRAC_1_SQRT_2;
        }

        let integrals = r_grid
            .iter()
            .map(|&r| {
                let (ao, s) = sto3g::molecule(&chem::diatomic(Element::C, r));
                let lowdin = chem::inverse_sqrt(&s);
                let c: Vec<Vec<f64>> = (0..10)
                    .map(|m| {
                        (0..10)
                            .map(|p| (0..10).map(|v| lowdin[m][v] * u[v][p]).sum())
                            .collect()
                    })
                    .collect();

                ao.transformed(&c).frozen_core(2)
            })
            .collect();

        Self {
            orbitals,
            z_total: 8,
            charges: 0..=0,
            r_grid,
            integrals,
            spin_orbit: false,
            stretch: false,
        }
    }

    fn n_spin_orbitals(&self) -> u8 {
        2 * self.orbitals.len() as u8
    }

    /// Every spin-orbital's bit.
    fn all_bits(&self) -> u64 {
        u64::MAX
            .checked_shr(64 - u32::from(self.n_spin_orbitals()))
            .unwrap_or(0)
    }

    /// Bitmask of distinct, in-range spin-orbitals.
    fn mask(&self, indices: &[u8]) -> Result<u64, Forbidden> {
        let mut mask = 0;
        for &i in indices {
            if i >= self.n_spin_orbitals() {
                return Err(Forbidden::Unbound);
            }

            let bit = 1u64 << i;
            if mask & bit != 0 {
                return Err(Forbidden::Malformed);
            }

            mask |= bit;
        }

        Ok(mask)
    }
}

// Also used by sto3g.
/// Spatial-orbital integrals, Hartree. Orbitals are orthonormal, unless straight from `sto3g`.
#[derive(Clone, Debug)]
pub struct Integrals {
    pub n: usize,
    /// Nuclear repulsion.
    pub e_nuc: f64,
    /// h_pq: kinetic + nuclear attraction. Row-major.
    pub h: Vec<f64>,
    /// (pq|rs), chemists' notation. Row-major.
    pub eri: Vec<f64>,
}

impl Integrals {
    fn h(&self, p: usize, q: usize) -> f64 {
        self.h[p * self.n + q]
    }

    fn eri(&self, p: usize, q: usize, r: usize, s: usize) -> f64 {
        self.eri[((p * self.n + q) * self.n + r) * self.n + s]
    }

    /// ⟨i|h|j⟩, over spin-orbitals.
    fn h_so(&self, i: u8, j: u8) -> f64 {
        if i % 2 != j % 2 {
            return 0.;
        }

        self.h(usize::from(i / 2), usize::from(j / 2))
    }

    /// ⟨ij||kl⟩ = ⟨ij|kl⟩ − ⟨ij|lk⟩, over spin-orbitals.
    fn anti(&self, i: u8, j: u8, k: u8, l: u8) -> f64 {
        // ⟨ij|kl⟩ = (ik|jl), if spins match.
        let phys = |i: u8, j: u8, k: u8, l: u8| {
            if i % 2 != k % 2 || j % 2 != l % 2 {
                return 0.;
            }

            let [i, j, k, l] = [i, j, k, l].map(|x| usize::from(x / 2));
            self.eri(i, k, j, l)
        };

        phys(i, j, k, l) - phys(i, j, l, k)
    }

    /// Into orbitals φ_p = Σ_μ c[μ][p]·χ_μ; `c` is square. One index at a time: O(n⁵).
    pub fn transformed(&self, c: &[Vec<f64>]) -> Self {
        let n = self.n;
        let flat = |[p, q, r, s]: [usize; 4]| ((p * n + q) * n + r) * n + s;

        let h = chem::indices::<2>(n)
            .map(|[p, q]| {
                chem::indices::<2>(n)
                    .map(|[m, v]| c[m][p] * c[v][q] * self.h(m, v))
                    .sum()
            })
            .collect();

        let mut eri = self.eri.clone();
        for k in 0..4 {
            eri = chem::indices::<4>(n)
                .map(|index| {
                    (0..n)
                        .map(|m| {
                            let mut from = index;
                            from[k] = m;
                            c[m][index[k]] * eri[flat(from)]
                        })
                        .sum()
                })
                .collect();
        }

        Self {
            n,
            e_nuc: self.e_nuc,
            h,
            eri,
        }
    }

    /// Folds orbitals 0..n_core, doubly occupied, into `h` and `e_nuc`, and drops them.
    pub fn frozen_core(&self, n_core: usize) -> Self {
        let n = self.n - n_core;

        let h = chem::indices::<2>(n)
            .map(|[p, q]| {
                let (p, q) = (p + n_core, q + n_core);
                let core: f64 = (0..n_core)
                    .map(|c| 2. * self.eri(p, q, c, c) - self.eri(p, c, c, q))
                    .sum();

                self.h(p, q) + core
            })
            .collect();
        let eri = chem::indices::<4>(n)
            .map(|[p, q, r, s]| self.eri(p + n_core, q + n_core, r + n_core, s + n_core))
            .collect();
        let e_core = (0..n_core).map(|c| 2. * self.h(c, c)).sum::<f64>()
            + chem::indices::<2>(n_core)
                .map(|[c, d]| 2. * self.eri(c, c, d, d) - self.eri(c, d, d, c))
                .sum::<f64>();

        Self {
            n,
            e_nuc: self.e_nuc + e_core,
            h,
            eri,
        }
    }
}

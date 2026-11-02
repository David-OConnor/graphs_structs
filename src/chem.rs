//! Long shot

//! Classical MD. This probably won't work?
//! todo: Discrete steps, or distance weightings?

mod sto3g;

use std::f64::consts::{FRAC_1_SQRT_2, TAU};

use crate::{Transformable, graph::Graph};

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

// ---------- Candidate states ----------
// State: small ints and enums; finite and `Eq + Hash`, so each is a graph vertex.
// Physics tables (energies, forces, integrals) live in `Context`.

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Spin {
    Up,
    Down,
}

/// |Orbital angular momentum projection on the bond axis|. Σ, Π, Δ for states; σ, π, δ for orbitals.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Lambda {
    Sigma,
    Pi,
    Delta,
}

/// Inversion symmetry. Homonuclear only.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Parity {
    Gerade,
    Ungerade,
}

/// Σ states only.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Reflection {
    Plus,
    Minus,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Element {
    H,
    C,
    N,
    O,
}

impl Element {
    pub fn z(self) -> u8 {
        match self {
            Self::H => 1,
            Self::C => 6,
            Self::N => 7,
            Self::O => 8,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct SpinCount {
    pub up: u8,
    pub down: u8,
}

/// Why a step is rejected. Mostly selection rules.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Forbidden {
    /// Target orbital or cell is full.
    Pauli,
    /// Nothing to move.
    Empty,
    /// ΔS ≠ 0, or an Ms change without spin-orbit coupling.
    Spin,
    /// g↔g or u↔u.
    Parity,
    /// |ΔΛ| > 1, or Σ+↔Σ−.
    Lambda,
    /// ΔJ outside the rule for this transition.
    Rotation,
    /// No transition dipole; e.g. vib-rot within a homonuclear state.
    Dipole,
    /// Ortho↔para, or a J level absent for identical nuclei.
    NuclearSpin,
    /// Above dissociation, or outside the context's grids, tables or limits.
    Unbound,
    NotAdjacent,
    /// State or step violates its invariants; e.g. a repeated index.
    Malformed,
}

// ---------- 1: Spectroscopic quantum numbers ----------

/// Electronic term ^(2S+1)Λ_(g/u)^(±), with constants for E(v, J).
#[derive(Clone, Debug)]
pub struct Term {
    /// E.g. 'X', 'a', 'B'.
    pub label: char,
    pub charge: i8,
    /// 2S + 1
    pub multiplicity: u8,
    pub lambda: Lambda,
    pub parity: Option<Parity>,
    pub reflection: Option<Reflection>,
    /// Index into `SpectroCtx::asymptotes`.
    pub asymptote: u8,
    /// Highest bound v.
    pub v_max: u8,
    // cm⁻¹. E = T_e + ω_e(v+½) − ω_e·x_e(v+½)² + B_v·J(J+1); B_v = B_e − α_e(v+½)
    pub t_e: f64,
    pub omega_e: f64,
    pub omega_e_x_e: f64,
    pub b_e: f64,
    pub alpha_e: f64,
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

// ---------- 2: Occupation bitmask (Fock space) ----------

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OrbitalKind {
    Core,
    Bonding,
    Antibonding,
    Nonbonding,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Center {
    /// Atomic orbital: valence bond, Hubbard.
    Atom(u8),
    /// Molecular orbital.
    Delocalized,
}

/// Spatial orbital; 2 spin-orbitals.
#[derive(Clone, Debug)]
pub struct Orbital {
    /// E.g. "σg2s", "πu2px".
    pub label: &'static str,
    pub lambda: Lambda,
    pub parity: Option<Parity>,
    pub center: Center,
    pub kind: OrbitalKind,
}

const BOHR_PER_ANGSTROM: f64 = 1. / 0.529_177_210_544;

/// Two atoms on the z axis, `r` Å apart; bohr.
fn diatomic(element: Element, r: f32) -> [(Element, [f64; 3]); 2] {
    [
        (element, [0.; 3]),
        (element, [0., 0., f64::from(r) * BOHR_PER_ANGSTROM]),
    ]
}

/// Every [i, j, ...] over 0..n, row-major.
fn indices<const D: usize>(n: usize) -> impl Iterator<Item = [usize; D]> {
    (0..n.pow(D as u32)).map(move |mut flat| {
        let mut index = [0; D];
        for slot in index.iter_mut().rev() {
            *slot = flat % n;
            flat /= n;
        }

        index
    })
}

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

        let h = indices::<2>(n)
            .map(|[p, q]| {
                indices::<2>(n)
                    .map(|[m, v]| c[m][p] * c[v][q] * self.h(m, v))
                    .sum()
            })
            .collect();

        let mut eri = self.eri.clone();
        for k in 0..4 {
            eri = indices::<4>(n)
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

        let h = indices::<2>(n)
            .map(|[p, q]| {
                let (p, q) = (p + n_core, q + n_core);
                let core: f64 = (0..n_core)
                    .map(|c| 2. * self.eri(p, q, c, c) - self.eri(p, c, c, q))
                    .sum();

                self.h(p, q) + core
            })
            .collect();
        let eri = indices::<4>(n)
            .map(|[p, q, r, s]| self.eri(p + n_core, q + n_core, r + n_core, s + n_core))
            .collect();
        let e_core = (0..n_core).map(|c| 2. * self.h(c, c)).sum::<f64>()
            + indices::<2>(n_core)
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
                let (ao, s) = sto3g::molecule(&diatomic(Element::H, r));
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
        use Lambda::{Pi, Sigma};
        use OrbitalKind::{Antibonding, Bonding};
        use Parity::{Gerade, Ungerade};

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
                let (ao, s) = sto3g::molecule(&diatomic(Element::C, r));
                let lowdin = inverse_sqrt(&s);
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
pub struct Determinant {
    /// Bit 2k: orbital k, up. Bit 2k + 1: down.
    pub occ: u64,
    /// Index into `DeterminantCtx::r_grid`.
    pub r: u8,
}

/// Even bits: spin up.
const UP_BITS: u64 = 0x5555_5555_5555_5555;

/// Set bit indices, ascending.
fn bits(mut occ: u64) -> Vec<u8> {
    let mut bits = Vec::new();
    while occ != 0 {
        bits.push(occ.trailing_zeros() as u8);
        occ &= occ - 1;
    }

    bits
}

/// Fermion sign of moving an electron i → a: −1 per occupied spin-orbital between them.
fn sign(occ: u64, i: u8, a: u8) -> f64 {
    let (lo, hi) = (i.min(a), i.max(a));
    let between = (1u64 << hi) - (1u64 << (lo + 1));

    (-1f64).powi((occ & between).count_ones() as i32)
}

impl Determinant {
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

        let occupied = bits(self.occ);
        let removed = bits(self.occ & !other.occ);
        let added = bits(other.occ & !self.occ);

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

                sign(self.occ, i, a) * (ints.h_so(i, a) + two)
            }
            (&[i, j], &[a, b]) => {
                let after_ia = self.occ ^ (1u64 << i) ^ (1u64 << a);

                sign(self.occ, i, a) * sign(after_ia, j, b) * ints.anti(i, j, a, b)
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

impl Transformable for Determinant {
    type Context = DeterminantCtx;
    type Step = DeterminantStep;
    type Invalid = Forbidden;

    fn steps(&self, ctx: &Self::Context) -> Vec<Self::Step> {
        use DeterminantStep::{Attach, DoubleExcite, Excite, Ionize, Stretch};

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

/// An H eigenstate on a block of the graph.
#[derive(Clone, Debug)]
pub struct Eigenstate {
    /// Hartree.
    pub energy: f64,
    pub r: u8,
    pub charge: i8,
    pub ms_x2: i8,
    /// (Graph vertex index, coefficient). Normalized.
    pub amplitudes: Vec<(usize, f64)>,
}

/// H's eigenstates on the graph, ascending. None without integrals.
///
/// Blocks are connected components over nonzero-weight edges; these split r, charge, Ms and
/// spatial symmetry with no labels. Each is diagonalized densely: fine for H2, not C2.
/// Missing vertices give truncated CI; energies are then upper bounds.
pub fn eigenstates(graph: &Graph<Determinant>, ctx: &DeterminantCtx) -> Option<Vec<Eigenstate>> {
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }

        i
    }

    let n = graph.vertices.len();
    let states: Vec<Determinant> = graph.vertices.iter().map(|v| v.instance).collect();

    let mut couplings = Vec::new();
    let mut parent: Vec<usize> = (0..n).collect();
    for edge in &graph.edges {
        let weight = states[edge.from].hamiltonian(states[edge.to], ctx)?;
        if weight.abs() > 1e-12 {
            couplings.push((edge.from, edge.to, weight));

            let (a, b) = (root(&mut parent, edge.from), root(&mut parent, edge.to));
            parent[a] = b;
        }
    }
    let block_of: Vec<usize> = (0..n).map(|i| root(&mut parent, i)).collect();

    // Indexed by block root.
    let mut members = vec![Vec::new(); n];
    let mut local = vec![0; n];
    for (i, &block) in block_of.iter().enumerate() {
        local[i] = members[block].len();
        members[block].push(i);
    }

    let mut matrices: Vec<Vec<Vec<f64>>> = members
        .iter()
        .map(|m| vec![vec![0.; m.len()]; m.len()])
        .collect();
    for (i, state) in states.iter().enumerate() {
        matrices[block_of[i]][local[i]][local[i]] = state.hamiltonian(*state, ctx)?;
    }
    for (from, to, weight) in couplings {
        let matrix = &mut matrices[block_of[from]];
        matrix[local[from]][local[to]] = weight;
        matrix[local[to]][local[from]] = weight;
    }

    let mut result = Vec::new();
    for (block, matrix) in members.iter().zip(matrices) {
        let Some(&first) = block.first() else {
            continue;
        };

        for (energy, vector) in eigh(matrix) {
            result.push(Eigenstate {
                energy,
                r: states[first].r,
                charge: states[first].charge(ctx),
                ms_x2: states[first].ms_x2(),
                amplitudes: block.iter().copied().zip(vector).collect(),
            });
        }
    }
    result.sort_by(|a, b| a.energy.total_cmp(&b.energy));

    Some(result)
}

/// Lowest energy per r bin, at `charge`: (Å, Hartree).
pub fn ground_curve(states: &[Eigenstate], ctx: &DeterminantCtx, charge: i8) -> Vec<(f32, f64)> {
    ctx.r_grid
        .iter()
        .enumerate()
        .filter_map(|(r, &r_angstrom)| {
            states
                .iter()
                .find(|s| usize::from(s.r) == r && s.charge == charge)
                .map(|s| (r_angstrom, s.energy))
        })
        .collect()
}

/// Lowest eigenstate in `start`'s block, by Lanczos: each H·v is one pass over the graph's
/// edges. For blocks too big for `eigenstates`. Weights use bin `r` for every vertex, so a graph
/// built on one bin (`stretch` off) serves them all. None without integrals.
pub fn lowest_state(
    graph: &Graph<Determinant>,
    ctx: &DeterminantCtx,
    start: usize,
    r: u8,
) -> Option<Eigenstate> {
    let states: Vec<Determinant> = graph
        .vertices
        .iter()
        .map(|v| Determinant { r, ..v.instance })
        .collect();
    let n = states.len();

    let diag: Vec<f64> = states
        .iter()
        .map(|s| s.hamiltonian(*s, ctx))
        .collect::<Option<_>>()?;
    let mut couplings = Vec::new();
    for edge in graph.edges.iter().filter(|e| e.from < e.to) {
        let weight = states[edge.from].hamiltonian(states[edge.to], ctx)?;
        if weight.abs() > 1e-12 {
            couplings.push((edge.from, edge.to, weight));
        }
    }

    let apply = |x: &[f64]| {
        let mut y: Vec<f64> = diag.iter().zip(x).map(|(d, x)| d * x).collect();
        for &(a, b, w) in &couplings {
            y[a] += w * x[b];
            y[b] += w * x[a];
        }

        y
    };
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();

    let mut basis = vec![vec![0.; n]];
    basis[0][start] = 1.;
    let (mut alpha, mut beta) = (Vec::new(), Vec::new());
    let (mut energy, mut ritz) = (0., Vec::new());
    for _ in 0..n.min(300) {
        let q = &basis[basis.len() - 1];
        let mut w = apply(q);
        alpha.push(dot(&w, q));

        // Full reorthogonalization; twice is enough.
        for _ in 0..2 {
            for q in &basis {
                let c = dot(&w, q);
                w.iter_mut().zip(q).for_each(|(w, q)| *w -= c * q);
            }
        }
        let b = dot(&w, &w).sqrt();

        energy = bisect(&alpha, &beta, 0);
        ritz = inverse_iterate(&alpha, &beta, energy);
        if (b * ritz[ritz.len() - 1]).abs() < 1e-9 {
            break;
        }

        beta.push(b);
        basis.push(w.iter().map(|w| w / b).collect());
    }

    let vector: Vec<f64> = (0..n)
        .map(|i| basis.iter().zip(&ritz).map(|(q, y)| q[i] * y).sum())
        .collect();

    Some(Eigenstate {
        energy,
        r,
        charge: states[start].charge(ctx),
        ms_x2: states[start].ms_x2(),
        amplitudes: vector
            .into_iter()
            .enumerate()
            .filter(|&(_, c)| c != 0.)
            .collect(),
    })
}

/// `lowest_state` per r bin: (Å, Hartree).
pub fn lowest_curve(
    graph: &Graph<Determinant>,
    ctx: &DeterminantCtx,
    start: usize,
) -> Vec<(f32, f64)> {
    (0..=u8::MAX)
        .zip(&ctx.r_grid)
        .filter_map(|(r, &r_angstrom)| {
            lowest_state(graph, ctx, start, r).map(|s| (r_angstrom, s.energy))
        })
        .collect()
}

/// Symmetric eigensolver; cyclic Jacobi. (Eigenvalue, eigenvector), ascending. O(n³) per sweep.
fn eigh(mut a: Vec<Vec<f64>>) -> Vec<(f64, Vec<f64>)> {
    let n = a.len();
    let scale: f64 = a.iter().flatten().map(|x| x * x).sum();

    let mut v = vec![vec![0.; n]; n];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.;
    }

    for _ in 0..50 {
        let off: f64 = indices::<2>(n)
            .filter(|[p, q]| p != q)
            .map(|[p, q]| a[p][q].powi(2))
            .sum();
        if off <= 1e-30 * scale {
            break;
        }

        for [p, q] in indices::<2>(n).filter(|[p, q]| p < q) {
            if a[p][q] == 0. {
                continue;
            }

            let theta = (a[q][q] - a[p][p]) / (2. * a[p][q]);
            let t = theta.signum() / (theta.abs() + theta.hypot(1.));
            let c = 1. / t.hypot(1.);
            let s = t * c;

            // A ← JᵀAJ; V ← VJ.
            for row in a.iter_mut().chain(v.iter_mut()) {
                (row[p], row[q]) = (c * row[p] - s * row[q], s * row[p] + c * row[q]);
            }
            let (upper, lower) = a.split_at_mut(q);
            for (x, y) in upper[p].iter_mut().zip(lower[0].iter_mut()) {
                (*x, *y) = (c * *x - s * *y, s * *x + c * *y);
            }
        }
    }

    let mut pairs: Vec<(f64, Vec<f64>)> = (0..n)
        .map(|i| (a[i][i], v.iter().map(|row| row[i]).collect()))
        .collect();
    pairs.sort_by(|x, y| x.0.total_cmp(&y.0));

    pairs
}

/// S^(−1/2): Löwdin orthogonalization.
fn inverse_sqrt(s: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let pairs = eigh(s.to_vec());

    (0..s.len())
        .map(|i| {
            (0..s.len())
                .map(|j| pairs.iter().map(|(l, v)| v[i] * v[j] / l.sqrt()).sum())
                .collect()
        })
        .collect()
}

// ---------- 3: Integer phase-space lattice (nuclear motion) ----------

/// ¹H₂ reduced mass, electron masses.
pub const MU_H2: f64 = 1.007_825_032 / 2. * 1_822.888_486;

/// ¹²C₂ reduced mass, electron masses.
pub const MU_C2: f64 = 6. * 1_822.888_486;

const FS_PER_AU_TIME: f64 = 0.024_188_843_265_864;

/// Linear through (x, y) points, ascending x; clamped at the ends.
fn interpolate(points: &[(f64, f64)], x: f64) -> f64 {
    let i = points.partition_point(|p| p.0 < x);
    if i == 0 {
        return points[0].1;
    }
    if i == points.len() {
        return points[i - 1].1;
    }

    let (a, b) = (points[i - 1], points[i]);

    a.1 + (b.1 - a.1) * (x - a.0) / (b.0 - a.0)
}

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
            let f: Vec<f64> = lattice.iter().map(|&r| interpolate(&slopes, r)).collect();

            // Integrate −F, so V and F agree.
            let mut v = vec![interpolate(&curve, lattice[0])];
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

// ---------- 4: Quantized superposition ----------

/// Eigenvalues below `x` of the tridiagonal matrix with `diag`, and `off[i]` coupling i and i + 1.
fn sturm_count(diag: &[f64], off: &[f64], x: f64) -> usize {
    let mut count = 0;
    let mut q = 0.;
    for (i, &d) in diag.iter().enumerate() {
        q = d - x - if i == 0 { 0. } else { off[i - 1].powi(2) / q };
        if q == 0. {
            q = 1e-300;
        }
        if q < 0. {
            count += 1;
        }
    }

    count
}

/// The k-th lowest eigenvalue, from 0; bisection.
fn bisect(diag: &[f64], off: &[f64], k: usize) -> f64 {
    let reach = 2. * off.iter().fold(0., |m: f64, o| m.max(o.abs()));
    let mut lo = diag.iter().copied().fold(f64::MAX, f64::min) - reach;
    let mut hi = diag.iter().copied().fold(f64::MIN, f64::max) + reach;
    for _ in 0..100 {
        let mid = (lo + hi) / 2.;
        if sturm_count(diag, off, mid) > k {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    (lo + hi) / 2.
}

/// Eigenvector for `e`; inverse iteration. Σ χ² = 1; first lobe positive.
fn inverse_iterate(diag: &[f64], off: &[f64], e: f64) -> Vec<f64> {
    let mut x = vec![1.; diag.len()];
    for _ in 0..3 {
        // Solve (T − e)·y = x; Thomas algorithm.
        let (mut sup, mut y) = (Vec::new(), Vec::new());
        let (mut prev_sup, mut prev_y) = (0., 0.);
        for (i, (&d, &xi)) in diag.iter().zip(&x).enumerate() {
            let below = if i == 0 { 0. } else { off[i - 1] };
            let pivot = d - e - below * prev_sup;
            let pivot = if pivot.abs() < 1e-300 { 1e-300 } else { pivot };
            prev_sup = off.get(i).copied().unwrap_or(0.) / pivot;
            prev_y = (xi - below * prev_y) / pivot;
            sup.push(prev_sup);
            y.push(prev_y);
        }
        for i in (0..y.len() - 1).rev() {
            y[i] -= sup[i] * y[i + 1];
        }

        // Scale before squaring; y can be near overflow.
        let peak = y.iter().fold(0., |m: f64, v| m.max(v.abs()));
        let norm = y.iter().map(|v| (v / peak).powi(2)).sum::<f64>().sqrt() * peak;
        x = y.iter().map(|v| v / norm).collect();
    }

    let peak = x.iter().fold(0., |m: f64, v| m.max(v.abs()));
    if x.iter()
        .find(|v| v.abs() > 0.01 * peak)
        .is_some_and(|&v| v < 0.)
    {
        x.iter_mut().for_each(|v| *v = -*v);
    }

    x
}

/// Bound vibrational states of one energy curve; J = 0. On a uniform grid the nuclear H is
/// tridiagonal: a path graph.
#[derive(Clone, Debug)]
pub struct VibLevels {
    /// Hartree.
    pub energies: Vec<f64>,
    /// Å
    pub grid: Vec<f64>,
    /// Per level, on `grid`. Σ χ² = 1.
    pub wavefunctions: Vec<Vec<f64>>,
    /// ⟨v|R|v'⟩, Å.
    pub r_matrix: Vec<Vec<f64>>,
}

impl VibLevels {
    /// 3-point finite differences, on `n` points spanning `curve`: (Å, Hartree), ascending; e.g.
    /// from `ground_curve`. `mu`: electron masses. Up to `max` levels, below the curve's last E.
    pub fn new(curve: &[(f32, f64)], mu: f64, n: usize, max: usize) -> Self {
        let points: Vec<(f64, f64)> = curve.iter().map(|&(r, e)| (f64::from(r), e)).collect();
        let (first, last) = (points[0], points[points.len() - 1]);
        let grid: Vec<f64> = (0..n)
            .map(|i| first.0 + (last.0 - first.0) * i as f64 / (n - 1) as f64)
            .collect();

        // −1/(2μ)·d²/dR² + V
        let h = (grid[1] - grid[0]) * BOHR_PER_ANGSTROM;
        let off = -1. / (2. * mu * h * h);
        let diag: Vec<f64> = grid
            .iter()
            .map(|&r| interpolate(&points, r) - 2. * off)
            .collect();

        let off = vec![off; n - 1];
        let n_levels = sturm_count(&diag, &off, last.1).min(max);
        let energies: Vec<f64> = (0..n_levels).map(|k| bisect(&diag, &off, k)).collect();
        let wavefunctions: Vec<Vec<f64>> = energies
            .iter()
            .map(|&e| inverse_iterate(&diag, &off, e))
            .collect();
        let r_matrix = wavefunctions
            .iter()
            .map(|a| {
                wavefunctions
                    .iter()
                    .map(|b| (0..n).map(|i| a[i] * b[i] * grid[i]).sum())
                    .collect()
            })
            .collect();

        Self {
            energies,
            grid,
            wavefunctions,
            r_matrix,
        }
    }
}

/// Phase unit is 2π/65536.
const PHASE_TURN: f64 = 65_536.;

#[derive(Clone, Copy, Debug)]
pub enum PacketBasis {
    Spectro(SpectroState),
    Determinant(Determinant),
}

pub struct WavepacketCtx {
    /// Sum of `pop`.
    pub quanta: u8,
    pub basis: Vec<PacketBasis>,
    /// Hartree.
    pub energies: Vec<f64>,
    /// fs per tick.
    pub dt: f64,
    /// Phase turned per tick: (E_k − E_0)·dt/ħ, rounded.
    pub omega: Vec<u16>,
    /// ⟨k|R|k'⟩, Å. Off-diagonals make ⟨R⟩ oscillate.
    pub r_matrix: Vec<Vec<f64>>,
    /// Tried from every state. [Tick]: `ConnectedOnly` gives free evolution.
    pub candidate_steps: Vec<WavepacketStep>,
}

impl WavepacketCtx {
    /// Over levels v = 0..k, as `SpectroState::Bound { term, v, j: 0 }`.
    pub fn vibrational(levels: &VibLevels, term: u8, k: usize, quanta: u8, dt: f64) -> Self {
        let dt_au = dt / FS_PER_AU_TIME;
        let energies = levels.energies[..k].to_vec();
        let omega = energies
            .iter()
            .map(|e| {
                let turns = ((e - energies[0]) * dt_au / TAU).rem_euclid(1.);
                ((turns * PHASE_TURN).round() as u32 % 65_536) as u16
            })
            .collect();

        Self {
            quanta,
            basis: (0..k)
                .map(|v| {
                    PacketBasis::Spectro(SpectroState::Bound {
                        term,
                        v: v as u8,
                        j: 0,
                    })
                })
                .collect(),
            energies,
            dt,
            omega,
            r_matrix: levels.r_matrix[..k]
                .iter()
                .map(|row| row[..k].to_vec())
                .collect(),
            candidate_steps: vec![WavepacketStep::Tick],
        }
    }
}

/// Superposition of K basis states. |c_k|² = pop_k / quanta; arg c_k = phase_k · 2π/65536.
///
/// u16 phases wrap, so `Tick` is a permutation: free evolution is a cycle, spanning a whole number
/// of beats. v = 0, 1 superposition: ⟨R⟩ oscillates at (E₁ − E₀)/ħ (bond vibration).
/// Two atom-centered dets on H2+: the electron tunnels between nuclei.
///
/// Canonical, so equal states hash equal: first populated phase = 0 (global phase); empty
/// components' phase = 0.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Wavepacket<const K: usize> {
    pub pop: [u8; K],
    pub phase: [u16; K],
}

impl<const K: usize> Wavepacket<K> {
    /// All population in component `k`.
    pub fn pure(k: usize, ctx: &WavepacketCtx) -> Self {
        let mut pop = [0; K];
        pop[k] = ctx.quanta;

        Self { pop, phase: [0; K] }
    }

    /// ⟨R⟩, Å.
    pub fn mean_r(&self, ctx: &WavepacketCtx) -> f64 {
        let amp = |k: usize| (f64::from(self.pop[k]) / f64::from(ctx.quanta)).sqrt();
        let angle = |k: usize| f64::from(self.phase[k]) * TAU / PHASE_TURN;

        indices::<2>(K)
            .map(|[a, b]| amp(a) * amp(b) * (angle(b) - angle(a)).cos() * ctx.r_matrix[a][b])
            .sum()
    }

    /// ⟨H⟩, Hartree.
    pub fn energy(&self, ctx: &WavepacketCtx) -> f64 {
        self.pop
            .iter()
            .zip(&ctx.energies)
            .map(|(&p, e)| f64::from(p) / f64::from(ctx.quanta) * e)
            .sum()
    }

    fn canonical(mut self) -> Self {
        let reference = self
            .pop
            .iter()
            .position(|&p| p > 0)
            .map_or(0, |k| self.phase[k]);
        for (phase, &pop) in self.phase.iter_mut().zip(&self.pop) {
            *phase = if pop == 0 {
                0
            } else {
                phase.wrapping_sub(reference)
            };
        }

        self
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WavepacketStep {
    /// c_k ← c_k·e^(−iE_k·dt/ħ)
    Tick,
    /// One population quantum; a pulse or coupling. A newly populated target takes the source's
    /// phase.
    Transfer { from: u8, to: u8 },
    /// Measurement or decoherence. Non-unitary.
    Collapse { to: u8 },
}

impl<const K: usize> Transformable for Wavepacket<K> {
    type Context = WavepacketCtx;
    type Step = WavepacketStep;
    type Invalid = Forbidden;

    fn steps(&self, ctx: &Self::Context) -> Vec<Self::Step> {
        ctx.candidate_steps
            .iter()
            .copied()
            .filter(|step| self.transform(step, ctx).is_ok())
            .collect()
    }

    /// 65536^(populated − 1) per population split: only K ≤ 2 is practical.
    fn instances(ctx: &Self::Context) -> Vec<Self> {
        let mut instances = Vec::new();
        for split in indices::<K>(usize::from(ctx.quanta) + 1) {
            if split.iter().sum::<usize>() != usize::from(ctx.quanta) {
                continue;
            }

            let pop = split.map(|p| p as u8);
            let free: Vec<usize> = (0..K).filter(|&k| pop[k] > 0).skip(1).collect();
            let count = 1u64
                .checked_shl(16 * free.len() as u32)
                .expect("Too many instances");
            for phases in 0..count {
                let mut phase = [0; K];
                for (i, &k) in free.iter().enumerate() {
                    phase[k] = (phases >> (16 * i)) as u16;
                }

                instances.push(Self { pop, phase });
            }
        }

        instances
    }

    fn transform(&self, step: &Self::Step, ctx: &Self::Context) -> Result<Self, Self::Invalid> {
        let total: u32 = self.pop.iter().map(|&p| u32::from(p)).sum();
        if ctx.omega.len() != K || total != u32::from(ctx.quanta) || self.canonical() != *self {
            return Err(Forbidden::Malformed);
        }

        let mut next = *self;
        match *step {
            WavepacketStep::Tick => {
                for (phase, &omega) in next.phase.iter_mut().zip(&ctx.omega) {
                    *phase = phase.wrapping_sub(omega);
                }
            }
            WavepacketStep::Transfer { from, to } => {
                let (from, to) = (usize::from(from), usize::from(to));
                if from >= K || to >= K {
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
            WavepacketStep::Collapse { to } => {
                let to = usize::from(to);
                if to >= K {
                    return Err(Forbidden::Unbound);
                }
                if self.pop[to] == 0 {
                    return Err(Forbidden::Empty);
                }

                next = Self::pure(to, ctx);
            }
        }

        Ok(next.canonical())
    }
}

// ---------- 5: Density quanta on a grid ----------

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

// ---------- 6: Molecule as a graph, with electron bookkeeping ----------

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Hybrid {
    None,
    Sp,
    Sp2,
    Sp3,
}

/// Local 2-center orbitals. Index into `BondEdge::occ`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondOrbital {
    Sigma,
    SigmaStar,
    /// 2 degenerate orbitals; up to 4 e-.
    Pi,
    PiStar,
}

/// Core electrons are implied by `element`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct AtomNode {
    pub element: Element,
    pub hybrid: Hybrid,
    /// Nonbonding valence electrons.
    pub lone: SpinCount,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct BondEdge {
    /// Indices into `ChemGraph::atoms`. a < b.
    pub atoms: (u8, u8),
    /// Indexed by `BondOrbital`. Bond order = (σ + π − σ* − π*) / 2.
    pub occ: [SpinCount; 4],
    /// Length bin.
    pub r: u8,
    /// Vibrational quantum number.
    pub v: u8,
}

pub struct ChemGraphCtx {
    /// Å
    pub r_min: f32,
    pub dr: f32,
    pub n_r: u8,
    pub v_max: u8,
}

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

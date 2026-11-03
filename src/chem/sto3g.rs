//! STO-3G integrals over Cartesian s and p Gaussians, by McMurchie-Davidson. Hartree atomic
//! units; bohr. After Helgaker, Jørgensen & Olsen, ch. 9.

use std::f64::consts::PI;
use crate::chem::determinant::Integrals;
use super::{Element, indices};

/// Contraction coefficients; the same for every element.
const C_1S: [f64; 3] = [0.154_328_97, 0.535_328_14, 0.444_634_54];
const C_2S: [f64; 3] = [-0.099_967_23, 0.399_512_83, 0.700_115_47];
const C_2P: [f64; 3] = [0.155_916_27, 0.607_683_72, 0.391_957_39];

/// 1s exponents; 2sp exponents past He.
fn exponents(element: Element) -> ([f64; 3], Option<[f64; 3]>) {
    match element {
        Element::H => ([3.425_250_91, 0.623_913_73, 0.168_855_40], None),
        Element::C => (
            [71.616_837, 13.045_096, 3.530_512_2],
            Some([2.941_249_4, 0.683_483_1, 0.222_289_9]),
        ),
        Element::N => (
            [99.106_169, 18.052_312, 4.885_660_2],
            Some([3.780_455_9, 0.878_496_6, 0.285_714_4]),
        ),
        Element::O => (
            [130.709_32, 23.808_861, 6.443_608_3],
            Some([5.033_151_3, 1.169_596_1, 0.380_389]),
        ),
    }
}

/// Contracted Σ c·x^l·y^m·z^n·exp(−α·|r − center|²).
struct Ao {
    center: [f64; 3],
    lmn: [i32; 3],
    /// (α, c); c includes normalization.
    prims: Vec<(f64, f64)>,
}

/// Per atom: 1s; then 2s, 2px, 2py, 2pz past He.
fn basis(atoms: &[(Element, [f64; 3])]) -> Vec<Ao> {
    let mut aos = Vec::new();
    for &(element, center) in atoms {
        let (core, valence) = exponents(element);

        let mut shells = vec![([0, 0, 0], core, C_1S)];
        if let Some(valence) = valence {
            shells.push(([0, 0, 0], valence, C_2S));
            shells.extend([[1, 0, 0], [0, 1, 0], [0, 0, 1]].map(|lmn| (lmn, valence, C_2P)));
        }

        for (lmn, exps, coeffs) in shells {
            let l = f64::from(lmn.iter().sum::<i32>());
            let prims = exps
                .iter()
                .zip(coeffs)
                .map(|(&a, c)| (a, c * (2. * a / PI).powf(0.75) * (4. * a).powf(l / 2.)))
                .collect();

            let mut ao = Ao { center, lmn, prims };
            let norm = overlap(&ao, &ao).sqrt();
            ao.prims.iter_mut().for_each(|p| p.1 /= norm);
            aos.push(ao);
        }
    }

    aos
}

/// Over the raw AOs of `atoms` (bohr); not orthogonal. Also returns the AO overlap matrix.
pub fn molecule(atoms: &[(Element, [f64; 3])]) -> (Integrals, Vec<Vec<f64>>) {
    let aos = basis(atoms);
    let n = aos.len();

    let s = aos
        .iter()
        .map(|a| aos.iter().map(|b| overlap(a, b)).collect())
        .collect();
    let h = indices::<2>(n)
        .map(|[a, b]| {
            let attraction: f64 = atoms
                .iter()
                .map(|&(element, c)| f64::from(element.z()) * nuclear(&aos[a], &aos[b], c))
                .sum();

            kinetic(&aos[a], &aos[b]) - attraction
        })
        .collect();

    // Unique (ab|cd) only: a ≥ b, c ≥ d, ab ≥ cd.
    let pairs: Vec<Vec<PrimPair>> = indices::<2>(n)
        .map(|[a, b]| prim_pairs(&aos[a], &aos[b]))
        .collect();
    let mut eri = vec![0.; n.pow(4)];
    for [a, b, c, d] in indices::<4>(n) {
        if b > a || d > c || (c, d) > (a, b) {
            continue;
        }

        let value = repulsion(&pairs[a * n + b], &pairs[c * n + d]);
        for [p, q, r, s] in [
            [a, b, c, d],
            [b, a, c, d],
            [a, b, d, c],
            [b, a, d, c],
            [c, d, a, b],
            [d, c, a, b],
            [c, d, b, a],
            [d, c, b, a],
        ] {
            eri[((p * n + q) * n + r) * n + s] = value;
        }
    }

    let e_nuc = indices::<2>(atoms.len())
        .filter(|[i, j]| i < j)
        .map(|[i, j]| {
            let ((a, ra), (b, rb)) = (atoms[i], atoms[j]);
            f64::from(a.z()) * f64::from(b.z()) / distance(ra, rb)
        })
        .sum();

    (Integrals { n, e_nuc, h, eri }, s)
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>().sqrt()
}

/// Hermite expansion coefficient E_t^ij, on one axis. q = A − B.
fn hermite(i: i32, j: i32, t: i32, q: f64, a: f64, b: f64) -> f64 {
    let p = a + b;
    if i < 0 || j < 0 || t < 0 || t > i + j {
        return 0.;
    }
    if i == 0 && j == 0 {
        return (-a * b / p * q * q).exp();
    }

    let e = |i, j, t| hermite(i, j, t, q, a, b);
    if j == 0 {
        e(i - 1, j, t - 1) / (2. * p) - b * q / p * e(i - 1, j, t)
            + f64::from(t + 1) * e(i - 1, j, t + 1)
    } else {
        e(i, j - 1, t - 1) / (2. * p)
            + a * q / p * e(i, j - 1, t)
            + f64::from(t + 1) * e(i, j - 1, t + 1)
    }
}

/// Σ over primitive pairs of c_a·c_b·f.
fn contract(a: &Ao, b: &Ao, f: impl Fn(f64, f64) -> f64) -> f64 {
    let mut sum = 0.;
    for &(alpha, ca) in &a.prims {
        for &(beta, cb) in &b.prims {
            sum += ca * cb * f(alpha, beta);
        }
    }

    sum
}

/// Primitive overlap, with b's powers `lb`.
fn overlap_prim(a: &Ao, alpha: f64, b: &Ao, beta: f64, lb: [i32; 3]) -> f64 {
    let p = alpha + beta;
    let e: f64 = (0..3)
        .map(|k| hermite(a.lmn[k], lb[k], 0, a.center[k] - b.center[k], alpha, beta))
        .product();

    e * (PI / p).powf(1.5)
}

fn overlap(a: &Ao, b: &Ao) -> f64 {
    contract(a, b, |alpha, beta| overlap_prim(a, alpha, b, beta, b.lmn))
}

fn kinetic(a: &Ao, b: &Ao) -> f64 {
    contract(a, b, |alpha, beta| {
        let shifted = |k: usize, d: i32| {
            let mut lb = b.lmn;
            lb[k] += d;
            overlap_prim(a, alpha, b, beta, lb)
        };
        let l = b.lmn.iter().sum::<i32>();

        beta * f64::from(2 * l + 3) * shifted(0, 0)
            - 2. * beta * beta * (0..3).map(|k| shifted(k, 2)).sum::<f64>()
            - 0.5
                * (0..3)
                    .map(|k| f64::from(b.lmn[k] * (b.lmn[k] - 1)) * shifted(k, -2))
                    .sum::<f64>()
    })
}

/// A primitive pair's product: a Hermite expansion about `center`.
struct PrimPair {
    p: f64,
    center: [f64; 3],
    /// c_a·c_b
    coeff: f64,
    /// Total angular momentum.
    l: usize,
    /// ([t, u, v], E_t·E_u·E_v)
    terms: Vec<([usize; 3], f64)>,
}

fn prim_pairs(a: &Ao, b: &Ao) -> Vec<PrimPair> {
    let mut pairs = Vec::new();
    for &(alpha, ca) in &a.prims {
        for &(beta, cb) in &b.prims {
            let p = alpha + beta;
            let tops: [usize; 3] = std::array::from_fn(|k| (a.lmn[k] + b.lmn[k]) as usize);
            let e = |k: usize, t: usize| {
                hermite(
                    a.lmn[k],
                    b.lmn[k],
                    t as i32,
                    a.center[k] - b.center[k],
                    alpha,
                    beta,
                )
            };

            pairs.push(PrimPair {
                p,
                center: std::array::from_fn(|k| (alpha * a.center[k] + beta * b.center[k]) / p),
                coeff: ca * cb,
                l: tops.iter().sum(),
                terms: indices::<3>(tops.iter().max().unwrap() + 1)
                    .filter(|tuv| (0..3).all(|k| tuv[k] <= tops[k]))
                    .map(|tuv| (tuv, e(0, tuv[0]) * e(1, tuv[1]) * e(2, tuv[2])))
                    .collect(),
            });
        }
    }

    pairs
}

/// ∫ a·b / |r − c|
fn nuclear(a: &Ao, b: &Ao, c: [f64; 3]) -> f64 {
    prim_pairs(a, b)
        .iter()
        .map(|pair| {
            let r = hermite_coulomb(
                pair.l,
                pair.p,
                std::array::from_fn(|k| pair.center[k] - c[k]),
            );
            let s = pair.l + 1;
            let sum: f64 = pair
                .terms
                .iter()
                .map(|&([t, u, v], e)| e * r[(t * s + u) * s + v])
                .sum();

            pair.coeff * 2. * PI / pair.p * sum
        })
        .sum()
}

/// (ab|cd), over each pair's primitive products.
fn repulsion(ab: &[PrimPair], cd: &[PrimPair]) -> f64 {
    let mut sum = 0.;
    for x in ab {
        for y in cd {
            let (p, q) = (x.p, y.p);
            let l = x.l + y.l;
            let s = l + 1;
            let r = hermite_coulomb(
                l,
                p * q / (p + q),
                std::array::from_fn(|k| x.center[k] - y.center[k]),
            );

            let mut inner = 0.;
            for &([t, u, v], ex) in &x.terms {
                for &([tau, nu, phi], ey) in &y.terms {
                    let sign = if (tau + nu + phi) % 2 == 0 { 1. } else { -1. };
                    inner += ex * ey * sign * r[((t + tau) * s + u + nu) * s + v + phi];
                }
            }

            sum += x.coeff * y.coeff * 2. * PI.powf(2.5) / (p * q * (p + q).sqrt()) * inner;
        }
    }

    sum
}

/// Hermite Coulomb integrals R_tuv, for t + u + v ≤ l; index (t·s + u)·s + v, s = l + 1.
/// `pc`: from the charge (or Hermite) center to the other center.
fn hermite_coulomb(l: usize, p: f64, pc: [f64; 3]) -> Vec<f64> {
    let s = l + 1;
    let at = |n: usize, t: usize, u: usize, v: usize| ((n * s + t) * s + u) * s + v;
    let boys = boys(l, p * pc.iter().map(|x| x * x).sum::<f64>());

    // R^n_tuv from R^(n+1), n descending.
    let mut r = vec![0.; s.pow(4)];
    for n in (0..=l).rev() {
        for [t, u, v] in indices::<3>(s) {
            if t + u + v > l - n {
                continue;
            }

            let step = |k: usize, m: usize, back: [usize; 3], back2: [usize; 3]| {
                let two_back = if m > 1 {
                    (m - 1) as f64 * r[at(n + 1, back2[0], back2[1], back2[2])]
                } else {
                    0.
                };

                two_back + pc[k] * r[at(n + 1, back[0], back[1], back[2])]
            };
            r[at(n, t, u, v)] = if t > 0 {
                step(0, t, [t - 1, u, v], [t.saturating_sub(2), u, v])
            } else if u > 0 {
                step(1, u, [t, u - 1, v], [t, u.saturating_sub(2), v])
            } else if v > 0 {
                step(2, v, [t, u, v - 1], [t, u, v.saturating_sub(2)])
            } else {
                (-2. * p).powi(n as i32) * boys[n]
            };
        }
    }
    r.truncate(s.pow(3));

    r
}

/// Boys functions F_0..=F_l at t. Top by series or asymptote, then downward recursion.
fn boys(l: usize, t: f64) -> Vec<f64> {
    let top = if t > 50. {
        let double_factorial: f64 = (1..=l).map(|k| (2 * k - 1) as f64).product();
        double_factorial / 2f64.powi(l as i32 + 1) * (PI / t.powi(2 * l as i32 + 1)).sqrt()
    } else {
        // e^(−t)·Σ (2t)^k / ((2l+1)(2l+3)···(2l+2k+1)); terms are all positive.
        let mut term = 1. / (2 * l + 1) as f64;
        let mut sum = term;
        for k in 1.. {
            term *= 2. * t / (2 * l + 2 * k + 1) as f64;
            sum += term;
            if term <= sum * 1e-17 {
                break;
            }
        }

        (-t).exp() * sum
    };

    let mut f = vec![top; l + 1];
    for n in (0..l).rev() {
        f[n] = (2. * t * f[n + 1] + (-t).exp()) / (2 * n + 1) as f64;
    }

    f
}

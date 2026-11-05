use std::{fmt::Debug, hash::Hash};

use crate::{
    Transformable,
    chem::{
        self, Element, SpinCount, VibLevels, Wavepacket, WavepacketCtx, WavepacketStep,
        determinant::{DeterminantCtx, SlaterDeterminant},
        phase_cell::{PhaseCellCtx, PhaseCellStep},
        wave_packets::{WavePackets, WavePacketsCtx, WavePacketsStep},
    },
    ex_sinking_star::{SokobanContext, SokobanEx0},
    graph::Graph,
};

#[derive(Clone, Debug)]
pub struct Field {
    pub name: String,
    pub value: String,
}

impl Field {
    pub fn new(name: impl Into<String>, value: impl Debug) -> Self {
        Self {
            name: name.into(),
            value: format!("{value:#?}"),
        }
    }
}

macro_rules! fields {
    ($value:ident; $($field:ident),* $(,)?) => {
        vec![$(Field::new(stringify!($field), &$value.$field)),*]
    };
}

pub struct StateView {
    pub label: String,
    pub category: String,
    pub fields: Vec<Field>,
    pub observables: Vec<Field>,
}

pub struct TransitionView {
    pub from: usize,
    pub to: usize,
    pub step: String,
}

pub struct Snapshot {
    pub states: Vec<StateView>,
    pub transitions: Vec<TransitionView>,
    pub context: Vec<Field>,
    pub truncated: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Model {
    Sokoban,
    HydrogenMo,
    HydrogenAo,
    Carbon,
    Nuclear,
    Wavepacket,
    Electronic,
}

impl Model {
    pub const ALL: [Self; 7] = [
        Self::Sokoban,
        Self::HydrogenMo,
        Self::HydrogenAo,
        Self::Carbon,
        Self::Nuclear,
        Self::Wavepacket,
        Self::Electronic,
    ];

    pub fn title(self) -> &'static str {
        match self {
            Self::Sokoban => "Sinking Star",
            Self::HydrogenMo => "H₂ · molecular orbitals",
            Self::HydrogenAo => "H₂ · atomic orbitals",
            Self::Carbon => "C₂ · determinants",
            Self::Nuclear => "Nuclear phase space",
            Self::Wavepacket => "Vibrational wavepacket",
            Self::Electronic => "Electronic wavepackets",
        }
    }

    pub fn description(self) -> &'static str {
        match self {
            Self::Sokoban => {
                "Move the monk, transform the dragon, and push its crystal form on a 5 × 5 board."
            }
            Self::HydrogenMo => {
                "STO-3G determinants: excitation, ionization, attachment, and bond stretching."
            }
            Self::HydrogenAo => {
                "The same H₂ system in an atom-centered basis. Inspect covalent and ionic occupations."
            }
            Self::Carbon => {
                "Frozen-core C₂ determinants at 1.24 Å. Explore spin-conserving electronic transitions."
            }
            Self::Nuclear => {
                "Bond position and momentum on an integer lattice, using an H₂ ground-state curve."
            }
            Self::Wavepacket => {
                "Two H₂ vibrational levels. Follow phase ticks, population transfer, collapse, and mean bond length."
            }
            Self::Electronic => {
                "Two H₂ energy eigenstates from minimal-basis FCI. Inspect populations, phases, and density."
            }
        }
    }

    pub fn build(self, limit: usize) -> Result<Snapshot, String> {
        match self {
            Self::Sokoban => {
                let ctx = SokobanContext {
                    grid_size: (5, 5),
                    grid_water_tiles: Vec::new(),
                };

                Ok(project(
                    SokobanEx0::new((2, 2), (2, 3), false),
                    &ctx,
                    limit,
                    fields!(ctx; grid_size, grid_water_tiles),
                    |s, _| StateView {
                        label: format!("M{},{} D{},{}", s.monk.0, s.monk.1, s.dragon.0, s.dragon.1),
                        category: if s.dragon_is_xtal {
                            "Crystal"
                        } else {
                            "Dragon"
                        }
                        .into(),
                        fields: fields!(s; monk, dragon, dragon_is_xtal),
                        observables: Vec::new(),
                    },
                ))
            }
            Self::HydrogenMo | Self::HydrogenAo | Self::Carbon => {
                let ctx = match self {
                    Self::HydrogenMo => DeterminantCtx::h2_mo(vec![0.6, 0.74, 1., 1.5, 2.]),
                    Self::HydrogenAo => DeterminantCtx::h2_ao(vec![0.6, 0.74, 1., 1.5, 2.]),
                    _ => DeterminantCtx::c2(vec![1.24]),
                };
                let start = SlaterDeterminant {
                    occ: if self == Self::Carbon { 0xff } else { 3 },
                    r: 0,
                };

                Ok(project(
                    start,
                    &ctx,
                    limit,
                    determinant_context(&ctx),
                    |s, c| {
                        let mut state = StateView {
                            label: format!("{:x} · r{}", s.occ, s.r),
                            category: format!("Charge {:+}", s.charge(c)),
                            fields: vec![
                                Field::new("occ (bit mask)", format!("{:#018b}", s.occ)),
                                Field::new("r", s.r),
                            ],
                            observables: vec![
                                Field::new("electrons", s.n_elec()),
                                Field::new("charge", s.charge(c)),
                                Field::new("2 × Ms", s.ms_x2()),
                                Field::new("R (Å)", c.r_grid[s.r as usize]),
                                Field::new("Hᵢᵢ (Hartree)", s.hamiltonian(*s, c)),
                            ],
                        };
                        for (k, orbital) in c.orbitals.iter().enumerate() {
                            state.observables.push(Field::new(
                                format!("{} (up, down)", orbital.label),
                                ((s.occ >> (2 * k)) & 1, (s.occ >> (2 * k + 1)) & 1),
                            ));
                        }
                        if self != Self::HydrogenAo {
                            state
                                .observables
                                .push(Field::new("bond order", s.bond_order(c)));
                        }

                        state
                    },
                ))
            }
            Self::Nuclear => {
                let curve = hydrogen_curve();
                let mut ctx = PhaseCellCtx::from_curves(&[curve], chem::MU_H2, 0.02, 0.1, 12);
                ctx.candidate_steps
                    .extend([PhaseCellStep::Kick(-1), PhaseCellStep::Kick(1)]);

                Ok(project(
                    ctx.at(0, 0.9),
                    &ctx,
                    limit,
                    fields!(ctx; r_min, dr, n_r, p_max, force, potential, mu, dt, candidate_steps),
                    |s, c| StateView {
                        label: format!("r{} p{:+}", s.r, s.p),
                        category: format!("Surface {}", s.surface),
                        fields: fields!(s; surface, r, p),
                        observables: vec![
                            Field::new("R (Å)", s.bond_length(c)),
                            Field::new("energy (Hartree)", s.energy(c)),
                        ],
                    },
                ))
            }
            Self::Wavepacket => {
                let levels = VibLevels::new(&hydrogen_curve(), chem::MU_H2, 128, 4);
                if levels.energies.len() < 2 {
                    return Err("The H₂ curve did not yield two bound vibrational levels.".into());
                }

                let mut ctx = WavepacketCtx::vibrational(&levels, 0, 2, 4, 0.1);
                ctx.candidate_steps.extend([
                    WavepacketStep::Transfer { from: 0, to: 1 },
                    WavepacketStep::Transfer { from: 1, to: 0 },
                    WavepacketStep::Collapse { to: 0 },
                    WavepacketStep::Collapse { to: 1 },
                ]);
                let context =
                    fields!(ctx; quanta, basis, energies, dt, omega, r_matrix, candidate_steps);

                Ok(project(
                    Wavepacket::<2> {
                        pop: [2, 2],
                        phase: [0, 0],
                    },
                    &ctx,
                    limit,
                    context,
                    |s, c| StateView {
                        label: format!("{:?} · {}", s.pop, s.phase[1]),
                        category: format!("Population {}", s.pop[0]),
                        fields: fields!(s; pop, phase),
                        observables: vec![
                            Field::new("energy (Hartree)", s.energy(c)),
                            Field::new("mean R (Å)", s.mean_r(c)),
                        ],
                    },
                ))
            }
            Self::Electronic => {
                let atoms = chem::diatomic(Element::H, 0.74);
                let electrons = SpinCount { up: 1, down: 1 };
                let mut ctx = WavePacketsCtx::new(&atoms, electrons, 2, 4, 32, 0.1)
                    .map_err(|e| format!("{e:?}"))?;
                ctx.candidate_steps.extend([
                    WavePacketsStep::Transfer { from: 0, to: 1 },
                    WavePacketsStep::Transfer { from: 1, to: 0 },
                    WavePacketsStep::Collapse { to: 0 },
                    WavePacketsStep::Collapse { to: 1 },
                ]);
                let mut context = fields!(ctx; energies, vectors, determinants, coefficients, quanta, phase_bins, omega, dt, candidate_steps);
                context.extend([
                    Field::new("atoms (bohr)", &atoms),
                    Field::new("electrons", electrons),
                    Field::new("basis", "STO-3G"),
                ]);
                let start =
                    WavePackets::new(vec![2, 2], vec![0, 0], &ctx).map_err(|e| format!("{e:?}"))?;

                Ok(project(start, &ctx, limit, context, |s, c| StateView {
                    label: format!("{:?} · {}", s.pop, s.phase[1]),
                    category: format!("Population {}", s.pop[0]),
                    fields: fields!(s; pop, phase),
                    observables: vec![
                        Field::new("energy (Hartree)", s.energy(c)),
                        Field::new("density at midpoint (e/bohr³)", s.density([0., 0., 0.], c)),
                        Field::new("density matrix", s.density_matrix(c)),
                    ],
                }))
            }
        }
    }
}

fn hydrogen_curve() -> Vec<(f32, f64)> {
    let mut ctx = DeterminantCtx::h2_mo((0..31).map(|i| 0.4 + i as f32 * 0.06).collect());
    ctx.stretch = false;
    ctx.charges = 0..=0;
    let graph = Graph::build(
        crate::graph::AllOrConnected::ConnectedOnly(SlaterDeterminant { occ: 3, r: 0 }),
        &ctx,
    );

    chem::lowest_curve(&graph, &ctx, 0)
}

fn determinant_context(ctx: &DeterminantCtx) -> Vec<Field> {
    let mut fields = fields!(ctx; orbitals, z_total, charges, r_grid, spin_orbit, stretch);
    for (r, integrals) in ctx.integrals.iter().enumerate() {
        fields.extend([
            Field::new(format!("integrals[{r}].n"), integrals.n),
            Field::new(format!("integrals[{r}].e_nuc"), integrals.e_nuc),
            Field::new(format!("integrals[{r}].h"), &integrals.h),
            Field::new(format!("integrals[{r}].eri"), &integrals.eri),
        ]);
    }

    fields
}

fn project<T>(
    start: T,
    ctx: &T::Context,
    limit: usize,
    context: Vec<Field>,
    inspect: impl Fn(&T, &T::Context) -> StateView,
) -> Snapshot
where
    T: Transformable + Eq + Hash,
    T::Step: Debug,
{
    let (graph, truncated) = Graph::build_limited(start, ctx, limit, limit.saturating_mul(12));
    let states = graph
        .vertices
        .iter()
        .map(|v| inspect(&v.instance, ctx))
        .collect();
    let transitions = graph
        .edges
        .into_iter()
        .map(|e| TransitionView {
            from: e.from,
            to: e.to,
            step: format!("{:?}", e.step),
        })
        .collect();

    Snapshot {
        states,
        transitions,
        context,
        truncated,
    }
}

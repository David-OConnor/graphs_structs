mod ex_sinking_star;
mod graph;

use std::hash::Hash;

use crate::{
    ex_sinking_star::{SokobanContext, SokobanEx0},
    graph::{AllOrConnected, Graph},
};

// /// Describes a transformation between one instance of a struct and another
// struct TransformationRule {}

/// Will this become the core part of this project?
pub trait Transformable: Clone + Eq + Hash {
    // /// Enum rule for returning if a result was valid, or what the result was.
    // type Validity;

    // todo: What is this?
    type Context; // static level data; `()` for systems that don't need it

    /// Describes a single transformation step.
    type Step: Clone;
    type Invalid;

    /// Candidate steps; `transform` rejects invalid ones.
    fn steps(&self, ctx: &Self::Context) -> Vec<Self::Step>;

    /// All valid states for a finite context.
    fn instances(ctx: &Self::Context) -> Vec<Self>;

    /// Applies one step, or reports why it is invalid.
    // todo this sized trait??
    fn transform(&self, step: &Self::Step, ctx: &Self::Context) -> Result<Self, Self::Invalid>
    where
        Self: Sized;

    // /// Compares raw states. Do we want this approach of comparing arbitrary
    // /// steps, or should we focus on curating the limited set of (always valid) steps?
    // fn is_valid(&self, other: &Self) -> bool;

    // Should we use this, or the plain bool type?
    // fn is_valid_details(&self, other: &Self) -> Self::Validity;
}

fn main() {
    println!("Hello, rules?");

    // An avenue of approach for building a robust system:
    // - Try implementing this for various Sinking Star levels
    // - Try implementing this for various examples and exercises in the Graphs book (e.g. early chapters)

    // Init
    let state_0 = SokobanEx0::new((2, 2), (2, 3), false);

    // 0 -> 1: Valid: Transform
    let state_1 = SokobanEx0::new((2, 2), (2, 3), true);

    // 0 -> 2 or 1 -> 2: Invalid (Dragon can't move on its own)
    let state_2 = SokobanEx0::new((2, 2), (2, 4), false);

    // 0 -> 3: Invalid: Only crystal form can be pushed
    let state_3 = SokobanEx0::new((2, 3), (2, 4), false);

    // 0 -> 4: Invalid: Can't push and transform in the same step
    // 1 -> 4: Valid: Monk pushes crystal
    let state_4 = SokobanEx0::new((2, 3), (2, 4), true);

    // 0 -> 5: Invalid: Monk can't move 2 squares, dragon can't move on its own, or transform here
    // 4 -> 5: Valid: Monk moves.
    let state_5 = SokobanEx0::new((3, 3), (2, 4), true);

    let ctx = SokobanContext {
        grid_size: (5, 5),
        grid_water_tiles: Vec::new(),
    };
    let connected = Graph::build(AllOrConnected::ConnectedOnly(state_0.clone()), &ctx);
    let all = Graph::<SokobanEx0>::build(AllOrConnected::AllPossibleInstances, &ctx);

    for (from, to) in [
        (&state_0, &state_1),
        (&state_0, &state_2),
        (&state_1, &state_2),
        (&state_0, &state_3),
        (&state_0, &state_4),
        (&state_1, &state_4),
        (&state_0, &state_5),
        (&state_4, &state_5),
    ] {
        let step = all
            .edges
            .iter()
            .find(|edge| {
                all.vertices[edge.from].instance == *from && all.vertices[edge.to].instance == *to
            })
            .map(|edge| &edge.step);

        println!("{from:?} -> {to:?}: {step:?}");
    }

    println!(
        "Reachable: {} states, {} steps",
        connected.vertices.len(),
        connected.edges.len()
    );
    println!(
        "All: {} states, {} steps",
        all.vertices.len(),
        all.edges.len()
    );
}

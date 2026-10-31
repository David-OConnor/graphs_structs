mod ex_sinking_star;
mod graph;

use graph::{Graph, Vertex};


/// Describes a transformation between one instance of a struct and another
// struct TransformationRule {}


/// Will this become the core part of this project?
trait Transformable {
    /// Enum rule for returning if a result was valid, or what the result was.
    type Validity;

    /// Describes a singnle transformation step. Enum or struct?
    type Step;

    /// Instead of an option, should we use a StepResult type?
    // todo this sized trait??
    fn transform(&self, step: &Self::Step) -> Option<Self>
    where
        Self: Sized;

    /// Compares raw states. Do we want this approach of comparing arbitrary
    /// steps, or should we focus on curating the limited set of (always valid) steps?
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

    // 0 -> 3: Valid: Monk push dragon
    let state_3 = SokobanEx0::new((2, 3), (2, 4), false);

    // 0 -> 4: Invalid: Can't push and transform in the same step
    let state_4 = SokobanEx0::new((2, 3), (2, 4), true);

    // 0 -> 5: Invalid: Monk can't move 2 squares, dragon can't move on its own, or transform here
    // 4 -> 5: Valid: Monk moves.
    let state_5 = SokobanEx0::new((3, 3), (2, 4), true);
}

mod graph;

use graph::{Graph, Vertex};

/// Describes a transformation between one instance of a struct and another
struct TransformationRule {}

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
    fn is_valid(&self, other: &Self) -> bool;

    // Should we use this, or the plain bool type?
    fn is_valid_details(&self, other: &Self) -> Self::Validity;

    /// todo: Or add to graph?
    fn create_vertex(&self) -> Vertex;
}

/// Example use
type Posit = (u8, u8);

fn adjacent(own: Posit, step: Posit) -> bool {
    let x0 = step.0 as i16;
    let x1 = step.0 as i16;
    let y0 = step.1 as i16;
    let y1 = step.1 as i16;

    (x0 == x1 && (y0 - y1).abs() == 1) || (y0 == y1 && (x0 - x1).abs() == 1)
}

/// Example
struct SokobanEx0 {
    monk: Posit,
    dragon: Posit,
    dragon_is_xtal: bool,
    // Alternative: Grid struct, with each tile having state.
    grid_size: (usize, usize),
    grid_water_tiles: Vec<(usize, usize)>,
}

impl SokobanEx0 {
    fn new(monk: Posit, dragon: Posit, dragon_is_xtal: bool) -> Self {
        Self {
            monk,
            dragon,
            dragon_is_xtal,
            grid_size: (0, 0),
            grid_water_tiles: Vec::new(),
        }
    }
}

#[derive(Clone, Copy, PartialEq)]
enum MoveValidity {
    Valid,
    InvalidMonkMovement,
    InvalidPushAdjacency,
    InvalidXformAdjacency,
    InvalidOther,
}

#[derive(Clone, Copy, PartialEq)]
enum Direction {
    Up,
    Down,
    Left,
    Right,
}

enum SokobanStep {
    /// Also handles a push
    MoveMonk(Direction),
    Transform(Direction),
}

impl Transformable for SokobanEx0 {
    type Validity = MoveValidity;
    type Step = SokobanStep;

    fn transform(&self, step: &SokobanStep) -> Option<Self> {
        unimplemented!()
    }

    fn is_valid(&self, next: &Self) -> bool {
        // Adjacent to transform
        if !adjacent(self.monk, self.dragon) && (self.dragon_is_xtal != next.dragon_is_xtal) {
            return false;
        }

        // Monk can only move one step
        if !adjacent(self.monk, next.monk) {
            return false;
        }

        true
    }

    fn is_valid_details(&self, next: &Self) -> Self::Validity {
        unimplemented!()
    }

    fn create_vertex(&self) -> Vertex {
        unimplemented!()
    }
}

fn main() {
    println!("Hello, rules?");

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

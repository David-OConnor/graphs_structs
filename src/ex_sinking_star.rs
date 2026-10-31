//! An example set, using *Order of the Sinking Star* levels.

use crate::{Transformable, graph::Graph};


type Posit = (u8, u8);

fn adjacent(own: Posit, step: Posit) -> bool {
    let x0 = own.0 as i16;
    let x1 = step.0 as i16;
    let y0 = own.1 as i16;
    let y1 = step.1 as i16;

    (x0 == x1 && (y0 - y1).abs() == 1) || (y0 == y1 && (x0 - x1).abs() == 1)
}



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

    // fn is_valid(&self, next: &Self) -> bool {
        // Adjacent to transform
    //    if !adjacent(self.monk, self.dragon) && (self.dragon_is_xtal != next.dragon_is_xtal) {
    //        return false;
     //   }

        // Monk can only move one step
   //     if !adjacent(self.monk, next.monk) {
    //        return false;
    //    }

   //     true
    // }

    // fn is_valid_details(&self, next: &Self) -> Self::Validity {
    //     unimplemented!()
    // }
}
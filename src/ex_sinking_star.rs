//! An example set, using *Order of the Sinking Star* levels.

use crate::Transformable;

type Posit = (u8, u8);

fn adjacent(own: Posit, step: Posit) -> bool {
    let x0 = own.0 as i16;
    let x1 = step.0 as i16;
    let y0 = own.1 as i16;
    let y1 = step.1 as i16;

    (x0 == x1 && (y0 - y1).abs() == 1) || (y0 == y1 && (x0 - x1).abs() == 1)
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct SokobanEx0 {
    pub monk: Posit,
    pub dragon: Posit,
    pub dragon_is_xtal: bool,
}

pub struct SokobanContext {
    // Alternative: Grid struct, with each tile having state.
    /// Width and height; coordinates start at zero.
    pub grid_size: (usize, usize),
    /// Impassable for both characters.
    pub grid_water_tiles: Vec<(usize, usize)>,
}

impl SokobanContext {
    fn is_open(&self, position: Posit) -> bool {
        let position = (usize::from(position.0), usize::from(position.1));

        position.0 < self.grid_size.0
            && position.1 < self.grid_size.1
            && !self.grid_water_tiles.contains(&position)
    }
}

impl SokobanEx0 {
    pub fn new(monk: Posit, dragon: Posit, dragon_is_xtal: bool) -> Self {
        Self {
            monk,
            dragon,
            dragon_is_xtal,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MoveValidity {
    InvalidMonkMovement,
    InvalidPushAdjacency,
    InvalidXformAdjacency,
    DragonNotCrystal,
    MalformedState,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Direction {
    Up,
    Down,
    Left,
    Right,
}

impl Direction {
    fn apply(self, (x, y): Posit) -> Option<Posit> {
        match self {
            Self::Up => y.checked_sub(1).map(|y| (x, y)),
            Self::Down => y.checked_add(1).map(|y| (x, y)),
            Self::Left => x.checked_sub(1).map(|x| (x, y)),
            Self::Right => x.checked_add(1).map(|x| (x, y)),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SokobanStep {
    /// Also handles a crystal push.
    MoveMonk(Direction),
    Transform(Direction),
}

impl Transformable for SokobanEx0 {
    type Context = SokobanContext;
    type Step = SokobanStep;
    type Invalid = MoveValidity;

    fn steps(&self, ctx: &Self::Context) -> Vec<Self::Step> {
        use Direction::{Down, Left, Right, Up};

        [Up, Down, Left, Right]
            .into_iter()
            .flat_map(|direction| {
                [
                    SokobanStep::MoveMonk(direction),
                    SokobanStep::Transform(direction),
                ]
            })
            .filter(|step| self.transform(step, ctx).is_ok())
            .collect()
    }

    fn instances(ctx: &Self::Context) -> Vec<Self> {
        let mut positions = Vec::new();
        for x in 0..ctx.grid_size.0.min(usize::from(u8::MAX) + 1) {
            for y in 0..ctx.grid_size.1.min(usize::from(u8::MAX) + 1) {
                let position = (x as u8, y as u8);
                if ctx.is_open(position) {
                    positions.push(position);
                }
            }
        }

        let mut instances = Vec::new();
        for &monk in &positions {
            for &dragon in &positions {
                if monk == dragon {
                    continue;
                }

                for dragon_is_xtal in [false, true] {
                    instances.push(Self::new(monk, dragon, dragon_is_xtal));
                }
            }
        }

        instances
    }

    fn transform(&self, step: &SokobanStep, ctx: &Self::Context) -> Result<Self, Self::Invalid> {
        if self.monk == self.dragon || !ctx.is_open(self.monk) || !ctx.is_open(self.dragon) {
            return Err(MoveValidity::MalformedState);
        }

        let mut next = self.clone();
        match *step {
            SokobanStep::MoveMonk(direction) => {
                let monk = direction
                    .apply(self.monk)
                    .filter(|&position| ctx.is_open(position))
                    .ok_or(MoveValidity::InvalidMonkMovement)?;

                if monk == self.dragon {
                    if !self.dragon_is_xtal {
                        return Err(MoveValidity::DragonNotCrystal);
                    }

                    next.dragon = direction
                        .apply(self.dragon)
                        .filter(|&position| ctx.is_open(position))
                        .ok_or(MoveValidity::InvalidPushAdjacency)?;
                }

                next.monk = monk;
            }
            SokobanStep::Transform(direction) => {
                if !adjacent(self.monk, self.dragon)
                    || direction.apply(self.monk) != Some(self.dragon)
                {
                    return Err(MoveValidity::InvalidXformAdjacency);
                }

                next.dragon_is_xtal = !self.dragon_is_xtal;
            }
        }

        Ok(next)
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

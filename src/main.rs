/// Describes a transformation between one instance of a struct and another
struct TransformationRule {

}


trait Transformable {
    fn transformation_rule();

    fn is_valid(&self, other: &Self) -> bool;
}


/// Example use
type posit = (u8, u8);



/// Example
struct SokobanEx0 {
	monk: posit,
    dragon: posit,
    dragon_is_xtal: bool,
    // Alternative: Grid struct, with each tile having state.
    grid_size: (usize, usize),
    grid_water_tiles: Vec<(usize, usize)>,
}

impl SokobanEx0 {
    fn new(monk: posit, dragon: posit, dragon_is_xtal: bool) -> Self {
        Self {
            monk,
            dragon,
            dragon_is_xtal,
            grid_size: (0, 0),
            grid_water_tiles: Vec::new(),
        }
    }
}


impl Transformable for SokobanEx0 {
    fn transform(&self) -> Option<Self> {

    }

    fn is_valid(&self, other: &Self) -> bool {

    }
}

fn main() {
    println!("Hello, world!");

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



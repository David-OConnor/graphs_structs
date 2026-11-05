#![cfg_attr(target_os = "windows", windows_subsystem = "windows")]

mod chem;
mod ex_sinking_star;
mod explorer;
mod graph;
mod gui;
mod mol_dynamics;

// /// Describes a transformation between one instance of a struct and another
// struct TransformationRule {}

/// Will this become the core part of this project?
// pub trait Transformable: Clone + Eq + Hash {
pub trait Transformable: Clone {
    // /// Enum rule for returning if a result was valid, or what the result was.
    // type Validity;

    /// Static data; doesn't transform that's part of the implementing struct.
    /// Use `()` for systems that don't need it.
    type Context;

    /// Describes a single transformation step.
    type Step: Clone;
    type Invalid;

    /// Candidate steps; `transform` rejects invalid ones.
    fn steps(&self, ctx: &Self::Context) -> Vec<Self::Step>;

    /// All valid states for a finite context.
    fn instances(ctx: &Self::Context) -> Vec<Self>;

    /// Applies one step, or reports why it is invalid.
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
    gui::run();
}

use std::collections::HashMap;

use crate::Transformable;

/// Placeholder? Rename??
pub enum AllOrConnected<T: Transformable> {
    ConnectedOnly(T),
    AllPossibleInstances,
}

/// This graph represents all possible states of the system, and how they can
/// transform to each other.
///
/// Each vertex corresponds to a given instance.
/// Each edge corresponds to a transformation (step?) between one instance and the next
///
/// There is a vertex for every instance possibility, although perhaps we limit this to ones
/// that are not disconnected, given a starting instance (vertex), and the transformation
/// rules. (?)
pub struct Graph<T: Transformable> {
    pub vertices: Vec<Vertex<T>>,
    pub edges: Vec<Edge<T::Step>>,
}

impl<T: Transformable> Graph<T> {
    /// `all` determines if we only create vertices that are connected to a starting one, or if we
    /// include all possible instances. (Warning: For instances which include non-discrete values,
    /// the latter may not be viable due to having an astounding number of values for a given field,
    /// and this compounding!) Both modes require a finite state space.
    pub fn build(all: AllOrConnected<T>, ctx: &T::Context) -> Self {
        let mut graph = Self {
            vertices: Vec::new(),
            edges: Vec::new(),
        };
        let mut indices = HashMap::new();
        let initial = match all {
            AllOrConnected::ConnectedOnly(start) => vec![start],
            AllOrConnected::AllPossibleInstances => T::instances(ctx),
        };

        for state in initial {
            indices
                .entry(state.clone())
                .or_insert_with(|| graph.add_vertex(state));
        }

        let mut from = 0;
        while from < graph.vertices.len() {
            let state = graph.vertices[from].instance.clone();

            for step in state.steps(ctx) {
                let Ok(next) = state.transform(&step, ctx) else {
                    continue;
                };

                let to = *indices
                    .entry(next.clone())
                    .or_insert_with(|| graph.add_vertex(next));

                graph.edges.push(Edge { from, to, step });
            }

            from += 1;
        }

        graph
    }

    /// Returns the new vertex's index, for building edges later.
    pub fn add_vertex(&mut self, state: T) -> usize {
        self.vertices.push(Vertex { instance: state });

        self.vertices.len() - 1
    }
}

/// Each instance of a `Transformable` is a Vertex
pub struct Vertex<T: Transformable> {
    pub instance: T,
}

/// Each transition from one specific `Transformable` instance to another is an edge.
pub struct Edge<S> {
    pub from: usize,
    pub to: usize,
    pub step: S,
}

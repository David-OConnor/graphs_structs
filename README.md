# Graphs and Structs

Lab: Mapping between graphs and structs, and defining rules for how struct instances can and cannot transform.

Keyword: *State-space graph*?

Provides the `Transformable` trait, which defines the rules for how a given rust type (e.g. struct) transforms. Once implemented, we can create a graph where each vertex is an instance, and each edge is a transform between instances.

We can create a graph of all possible instances, or only ones reachable from a given starting point, and rules. (i.e. no disconnected vertices)

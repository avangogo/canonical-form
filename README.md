# canonical-form

Algorithm to reduce combinatorial structures modulo isomorphism.

This can typically be used to to test if two graphs are isomorphic.

The algorithm manipulates its input as a black box by
the action of permutations
and by testing equallity with element of its orbit,
plus some user-defined functions
that help to break symmetries.

```rust
use canonical_form::Canonize;

// Simple Graph implementation as adjacency lists
#[derive(Ord, PartialOrd, PartialEq, Eq, Clone, Debug)]
struct Graph {
      adj: Vec<Vec<usize>>,
}


impl Graph {
   fn new(n: usize, edges: &[(usize, usize)]) -> Self {
       let mut adj = vec![Vec::new(); n];
       for &(u, v) in edges {
           adj[u].push(v);
           adj[v].push(u);
       }
       for list in &mut adj {
           list.sort() // Necessary to make the derived `==` correct
       }
       Graph { adj }
   }
}

// The Canonize trait allows to use the canonial form algorithms
impl Canonize for Graph {
   fn size(&self) -> usize {
       self.adj.len()
   }
   fn apply_morphism(&self, perm: &[usize]) -> Self {
       let mut adj = vec![Vec::new(); self.size()];
       for (i, nbrs) in self.adj.iter().enumerate() {
           adj[perm[i]] = nbrs.iter().map(|&u| perm[u]).collect();
           adj[perm[i]].sort();
       }
       Graph { adj }
   }
   fn invariant_neighborhood(&self, u: usize) -> Vec<Vec<usize>> {
       vec![self.adj[u].clone()]
   }
}

// Usage of library functions
// Two isomorphic graphs
let c5 = Graph::new(5, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]);
let other_c5 = Graph::new(5, &[(0, 2), (2, 1), (1, 4), (4, 3), (3, 0)]);
assert_eq!(c5.canonical(), other_c5.canonical());

// Non-isomorphic graphs
let p5 = Graph::new(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
assert!(c5.canonical() != p5.canonical());

// Recovering the permutation that gives the canonical form
let p = c5.morphism_to_canonical();
assert_eq!(c5.apply_morphism(&p), c5.canonical());

// Enumerating automorphisms
assert_eq!(c5.canonical().automorphisms().count(), 10)
```

License: MIT

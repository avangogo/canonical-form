# canonical-form

Algorithms to reduce combinatorial structures modulo isomorphism.

This can typically be use to to test if two graphs are isomorphic.
```rust
use canonical_form::*;

// Simple Graph implementation as adjacency list
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
let c5 = Graph::new(5, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]);
let other_c5 = Graph::new(5, &[(0, 2), (2, 1), (1, 4), (4, 3), (3, 0)]);
assert_eq!(canonical_form(&c5), canonical_form(&other_c5));

let p5 = Graph::new(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
assert!(canonical_form(&c5) != canonical_form(&p5));

let p = canonical_form_morphism(&c5);
assert_eq!(c5.apply_morphism(&p), canonical_form(&c5));
```

License: MIT

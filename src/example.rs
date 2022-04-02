//! Structure used in doc examples
use crate::Canonize;

/// Simple Graph implementation as adjacency lists
#[derive(Ord, PartialOrd, PartialEq, Eq, Clone, Debug)]
pub struct Graph {
    adj: Vec<Vec<usize>>,
}

impl Graph {
    /// Create a graph of size `n` with the edges in `edges`
    pub fn new(n: usize, edges: &[(usize, usize)]) -> Self {
        let mut adj = vec![Vec::new(); n];
        for &(u, v) in edges {
            adj[u].push(v);
            adj[v].push(u);
        }
        for nbrs in &mut adj {
            nbrs.sort_unstable() // Necessary to make the derived `==` correct
        }
        Graph { adj }
    }
}

impl Canonize for Graph {
    fn size(&self) -> usize {
        self.adj.len()
    }
    fn apply_morphism(&self, perm: &[usize]) -> Self {
        let mut adj = vec![Vec::new(); self.size()];
        for (i, nbrs) in self.adj.iter().enumerate() {
            adj[perm[i]] = nbrs.iter().map(|&u| perm[u]).collect();
            adj[perm[i]].sort_unstable();
        }
        Graph { adj }
    }
    fn invariant_neighborhood(&self, u: usize) -> Vec<Vec<usize>> {
        vec![self.adj[u].clone()]
    }
}

//! Coverage for colored `invariant_neighborhood` hints.
//!
//! Colors are opaque labels, so two hints of color `a` must never be confused
//! with one hint of color `2a`. The refiner fingerprints them with a hash, so
//! these tests check both that the canonical form stays isomorphism invariant
//! and that the hint colors are actually used to tell structures apart.

mod common;

use canonical_form::Canonize;
use common::{check_contract, check_sampled};
use rand::{RngExt, SeedableRng, rngs::StdRng};

/// A graph whose edges carry a color, exposed as `(neighbor, color)` hints.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct ColoredGraph {
    n: usize,
    adj: Vec<Vec<(usize, u64)>>,
}

impl ColoredGraph {
    fn new(n: usize, edges: &[(usize, usize, u64)]) -> Self {
        let mut adj = vec![Vec::new(); n];
        for &(u, v, c) in edges {
            adj[u].push((v, c));
            adj[v].push((u, c));
        }
        for l in &mut adj {
            l.sort_unstable();
            l.dedup();
        }
        ColoredGraph { n, adj }
    }

    fn edges(&self) -> Vec<(usize, usize, u64)> {
        let mut e = Vec::new();
        for (u, l) in self.adj.iter().enumerate() {
            for &(v, c) in l {
                e.push((u, v, c));
            }
        }
        e
    }
}

impl Canonize for ColoredGraph {
    fn size(&self) -> usize {
        self.n
    }
    fn apply_morphism(&self, p: &[usize]) -> Self {
        let edges: Vec<_> = self
            .edges()
            .into_iter()
            .map(|(u, v, c)| (p[u], p[v], c))
            .collect();
        ColoredGraph::new(self.n, &edges)
    }
    fn invariant_neighborhood(&self, u: usize) -> impl Iterator<Item = (usize, u64)> {
        self.adj[u].iter().copied()
    }
}

#[test]
fn invariant_on_colored_circulants() {
    let mut rng = StdRng::seed_from_u64(2024);
    for &(n, steps) in &[
        (12usize, &[1usize, 3][..]),
        (14, &[1, 4]),
        (16, &[1, 2, 7]),
        (20, &[1, 5]),
    ] {
        // Each step gets its own color, so every vertex looks the same and
        // only individualization can break the symmetry.
        let edges: Vec<_> = (0..n)
            .flat_map(|u| {
                steps
                    .iter()
                    .enumerate()
                    .map(move |(i, &s)| (u, (u + s) % n, i as u64))
            })
            .collect();
        check_sampled(
            &format!("circulant {n}"),
            &ColoredGraph::new(n, &edges),
            50,
            &mut rng,
        );
    }
}

#[test]
fn colors_are_opaque_labels() {
    // Same underlying graph, different edge colors: the colors alone must make
    // the two structures inequivalent.
    let path = [(0usize, 1usize), (1, 2), (2, 3)];
    let a = ColoredGraph::new(4, &path.map(|(u, v)| (u, v, 1)));
    let b = ColoredGraph::new(4, &[(0, 1, 1), (1, 2, 2), (2, 3, 1)]);
    assert_ne!(a.canonical(), b.canonical());

    // Two hints of color 21 must not look like one hint of color 42.
    let c = ColoredGraph::new(3, &[(0, 1, 21), (0, 2, 21)]);
    let d = ColoredGraph::new(3, &[(0, 1, 42), (0, 2, 21)]);
    assert_ne!(c.canonical(), d.canonical());
}

#[test]
fn the_contract_holds_on_small_colored_graphs() {
    let mut rng = StdRng::seed_from_u64(31337);
    for trial in 0..60 {
        let n = 4 + trial % 2;
        let ncolors = 2 + (trial % 5) as u64;
        let mut edges = Vec::new();
        for u in 0..n {
            for v in 0..u {
                if rng.random_bool(0.4) {
                    // Colors are spread out so that small sums of colors are
                    // not colors themselves.
                    edges.push((u, v, 7 * rng.random_range(0..ncolors)));
                }
            }
        }
        check_contract(&ColoredGraph::new(n, &edges));
    }
}

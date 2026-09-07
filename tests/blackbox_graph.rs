//! Canonization with no invariants at all.
//!
//! [`NaiveGraph`] implements only `size` and `apply_morphism`, leaving
//! `invariant_neighborhood` and `invariant_color` at their defaults. The
//! refiner therefore never splits a cell: the search individualizes vertices
//! one at a time down to a discrete partition, and the *only* thing that keeps
//! it from walking the whole symmetric group is the automorphism pruning.
//!
//! That inverts the usual cost, and it is what sets the sizes below. The work
//! is about `n! / |Aut|`, so:
//!
//! | graph | leaves at n = 10 |
//! | --- | --- |
//! | clique, empty (`|Aut| = n!`) | 18 |
//! | cycle (`|Aut| = 2n`) | 181 445 |
//! | asymmetric (`|Aut| = 1`) | 3 628 800 |
//!
//! So the symmetric graphs are the cheap ones here and can be tested large,
//! while anything with a small automorphism group has to stay small: an
//! asymmetric graph fills the leaf map with one image per permutation, 720 at
//! n = 6 and 40 320 at n = 8.

mod common;

use canonical_form::Canonize;
use common::{check_contract, permutations, relabel};
use rand::rngs::StdRng;
use rand::{Rng, RngExt, SeedableRng};

/// A graph that gives the algorithm nothing beyond the group action.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct NaiveGraph {
    n: usize,
    /// Sorted, so that `==` and `Ord` compare labelled graphs.
    adj: Vec<Vec<usize>>,
}

impl NaiveGraph {
    fn new(n: usize, edges: &[(usize, usize)]) -> Self {
        let mut adj = vec![Vec::new(); n];
        for &(u, v) in edges {
            adj[u].push(v);
            adj[v].push(u);
        }
        for list in &mut adj {
            list.sort_unstable();
            list.dedup();
        }
        NaiveGraph { n, adj }
    }
}

impl Canonize for NaiveGraph {
    fn size(&self) -> usize {
        self.n
    }

    fn apply_morphism(&self, p: &[usize]) -> Self {
        let mut adj = vec![Vec::new(); self.n];
        for (u, neighbors) in self.adj.iter().enumerate() {
            adj[p[u]] = neighbors.iter().map(|&v| p[v]).collect();
            adj[p[u]].sort_unstable();
        }
        NaiveGraph { n: self.n, adj }
    }
}

fn gnp<R: Rng>(n: usize, p: f64, rng: &mut R) -> NaiveGraph {
    let mut edges = Vec::new();
    for i in 0..n {
        for j in 0..i {
            if rng.random_bool(p) {
                edges.push((i, j));
            }
        }
    }
    NaiveGraph::new(n, &edges)
}

fn clique(n: usize) -> NaiveGraph {
    let mut edges = Vec::new();
    for i in 0..n {
        for j in 0..i {
            edges.push((i, j));
        }
    }
    NaiveGraph::new(n, &edges)
}

fn cycle(n: usize) -> NaiveGraph {
    NaiveGraph::new(n, &(0..n).map(|i| (i, (i + 1) % n)).collect::<Vec<_>>())
}

fn empty(n: usize) -> NaiveGraph {
    NaiveGraph::new(n, &[])
}

#[test]
fn a_large_automorphism_group_keeps_it_tractable() {
    // These would be 16! leaves without the pruning; they are a few dozen.
    let mut rng = StdRng::seed_from_u64(7);
    for n in [1, 2, 6, 10, 16] {
        for g in [clique(n), empty(n)] {
            let expected = g.canonical();
            for _ in 0..4 {
                assert_eq!(relabel(&g, &mut rng).canonical(), expected, "n = {n}");
            }
        }
    }
    // A cycle only has 2n automorphisms, so it grows as n!/2n: n = 8 is
    // already 2 525 leaves and n = 12 takes a minute.
    for n in [3, 4, 5, 6, 7] {
        let g = cycle(n);
        let expected = g.canonical();
        for _ in 0..4 {
            assert_eq!(relabel(&g, &mut rng).canonical(), expected, "cycle {n}");
        }
    }
}

#[test]
fn non_isomorphic_graphs_keep_different_canonical_forms() {
    // The pruning is the only thing cutting the search here, so this is where
    // it would show if it cut too much.
    let mut rng = StdRng::seed_from_u64(31415);
    for n in 3..=5 {
        let mut graphs: Vec<NaiveGraph> = (0..12).map(|_| gnp(n, 0.5, &mut rng)).collect();
        graphs.push(clique(n));
        graphs.push(empty(n));
        graphs.push(cycle(n));
        let perms = permutations(n);
        for (i, a) in graphs.iter().enumerate() {
            for b in &graphs[..i] {
                let isomorphic = perms.iter().any(|p| a.apply_morphism(p) == *b);
                assert_eq!(a.canonical() == b.canonical(), isomorphic, "n = {n}");
            }
        }
    }
}

#[test]
fn the_contract_holds_with_no_invariants() {
    let mut rng = StdRng::seed_from_u64(20260906);
    for n in 0..=5 {
        for _ in 0..3 {
            check_contract(&gnp(n, 0.5, &mut rng));
        }
        check_contract(&clique(n));
        check_contract(&empty(n));
    }
    for n in 3..=5 {
        check_contract(&cycle(n));
    }
    // One pass at n = 6, where an asymmetric graph really does put all 720
    // permutations in the leaf map: checking the whole contract there would
    // canonize every one of them.
    let g = gnp(6, 0.5, &mut rng);
    assert_eq!(relabel(&g, &mut rng).canonical(), g.canonical());
    // The typed forms, which nothing else here exercises.
    for sigma in 0..=2 {
        let g = gnp(5, 0.5, &mut rng);
        let phi = g.morphism_to_canonical_typed(sigma);
        assert_eq!(g.apply_morphism(&phi), g.canonical_typed(sigma));
    }
}

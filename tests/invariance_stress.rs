//! Stress test for the property the whole algorithm rests on:
//! `canonical` must return the *exact same* structure on isomorphic inputs.
//!
//! The graphs here are chosen to make the target-cell selection matter:
//! regular graphs (where the initial refinement does nothing) and disjoint
//! unions of identical components (where equal-sized cells tie constantly).

mod common;

use canonical_form::Canonize;
use canonical_form::example::Graph;
use common::check_sampled;
use rand::{SeedableRng, prelude::SliceRandom, rngs::StdRng};

/// Check that `g` and 50 random relabellings of it have the same canonical form.
fn check(name: &str, g: &Graph, rng: &mut StdRng) {
    check_sampled(name, g, 50, rng);
}

fn cycle(n: usize) -> Vec<(usize, usize)> {
    (0..n).map(|i| (i, (i + 1) % n)).collect()
}

/// `k` disjoint copies of the graph given by `edges` on `n` vertices.
fn disjoint_copies(n: usize, edges: &[(usize, usize)], k: usize) -> Graph {
    let mut all = Vec::new();
    for c in 0..k {
        all.extend(edges.iter().map(|&(u, v)| (u + c * n, v + c * n)));
    }
    Graph::new(n * k, &all)
}

fn hypercube(d: usize) -> Graph {
    let n = 1 << d;
    let mut edges = Vec::new();
    for u in 0..n {
        for b in 0..d {
            let v = u ^ (1 << b);
            if u < v {
                edges.push((u, v));
            }
        }
    }
    Graph::new(n, &edges)
}

/// Random `d`-regular-ish graph: a union of `d` random perfect matchings /
/// hamiltonian cycles, which 1-WL cannot split at all.
fn random_circulant(n: usize, steps: &[usize]) -> Graph {
    let mut edges = Vec::new();
    for u in 0..n {
        for &s in steps {
            let v = (u + s) % n;
            if u != v {
                edges.push((u, v));
            }
        }
    }
    Graph::new(n, &edges)
}

/// 4x4 rook's graph and the Shrikhande graph: the classic pair of
/// non-isomorphic strongly regular graphs that 1-WL cannot tell apart.
fn rook4() -> Graph {
    let mut edges = Vec::new();
    for u in 0..16 {
        for v in 0..u {
            if u / 4 == v / 4 || u % 4 == v % 4 {
                edges.push((u, v));
            }
        }
    }
    Graph::new(16, &edges)
}

fn shrikhande() -> Graph {
    let idx = |x: usize, y: usize| (x % 4) * 4 + (y % 4);
    let mut edges = Vec::new();
    for x in 0..4 {
        for y in 0..4 {
            for &(dx, dy) in &[(1, 0), (0, 1), (1, 1)] {
                let (u, v) = (idx(x, y), idx(x + dx, y + dy));
                if u < v {
                    edges.push((u, v));
                } else {
                    edges.push((v, u));
                }
            }
        }
    }
    edges.sort_unstable();
    edges.dedup();
    Graph::new(16, &edges)
}

#[test]
fn invariant_on_symmetric_graphs() {
    let mut rng = StdRng::seed_from_u64(0xca7_f00d);

    check(
        "petersen",
        &{
            let mut edges: Vec<_> = cycle(5);
            edges.extend((0..5).map(|i| (i + 5, (i + 2) % 5 + 5)));
            edges.extend((0..5).map(|i| (i, i + 5)));
            Graph::new(10, &edges)
        },
        &mut rng,
    );

    check("q3", &hypercube(3), &mut rng);
    check("q4", &hypercube(4), &mut rng);
    check("rook4", &rook4(), &mut rng);
    check("shrikhande", &shrikhande(), &mut rng);

    // Non-isomorphic 1-WL-equivalent pair must stay distinguishable.
    assert_ne!(rook4().canonical(), shrikhande().canonical());

    for &(n, k) in &[
        (3, 2),
        (3, 4),
        (4, 2),
        (4, 3),
        (5, 2),
        (5, 3),
        (5, 4),
        (6, 3),
        (7, 3),
    ] {
        check(
            &format!("{k} x C{n}"),
            &disjoint_copies(n, &cycle(n), k),
            &mut rng,
        );
    }

    for &(n, steps) in &[
        (12usize, &[1usize, 3][..]),
        (12, &[1, 2, 5]),
        (14, &[1, 4]),
        (16, &[1, 2, 7]),
        (18, &[1, 6]),
        (20, &[1, 5]),
    ] {
        check(
            &format!("circulant {n} {steps:?}"),
            &random_circulant(n, steps),
            &mut rng,
        );
    }
}

#[test]
fn invariant_on_random_regular_graphs() {
    let mut rng = StdRng::seed_from_u64(12345);
    // Unions of random hamiltonian cycles: regular, so refinement stalls
    // immediately and every split comes from individualization.
    for n in [8usize, 10, 12, 14, 16] {
        for d in [2usize, 3] {
            for _ in 0..5 {
                let mut edges = Vec::new();
                for _ in 0..d {
                    let mut perm: Vec<_> = (0..n).collect();
                    perm.shuffle(&mut rng);
                    for i in 0..n {
                        let (u, v) = (perm[i], perm[(i + 1) % n]);
                        edges.push(if u < v { (u, v) } else { (v, u) });
                    }
                }
                edges.sort_unstable();
                edges.dedup();
                let g = Graph::new(n, &edges);
                check(&format!("random {d}-regular on {n}"), &g, &mut rng);
            }
        }
    }
}

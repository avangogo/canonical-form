extern crate canonical_form;

use canonical_form::example::Graph;
use canonical_form::Canonize;
use rand::{prelude::SliceRandom, random, thread_rng};

fn random_graph(n: usize) -> Graph {
    let mut edges = Vec::new();
    for u in 0..n {
        for v in 0..u {
            if random() {
                edges.push((u, v));
            }
        }
    }
    Graph::new(n, &edges)
}

fn random_permutation(n: usize) -> Vec<usize> {
    let mut perm: Vec<_> = (0..n).collect();
    perm.shuffle(&mut thread_rng());
    perm
}

#[test]
fn canonical() {
    for n in 0..20 {
        let g = random_graph(n);
        let permuted_g = g.apply_morphism(&random_permutation(n));
        assert_eq!(g.canonical(), permuted_g.canonical());
    }
}

#[test]
fn morphism_to_canonical() {
    for n in 0..20 {
        let g = random_graph(n);
        let morphism = g.morphism_to_canonical();
        assert_eq!(g.canonical(), g.apply_morphism(&morphism));
    }
}

#[test]
fn automorphisms() {
    for n in 0..20 {
        let g = random_graph(n).canonical();
        for morphism in g.automorphisms() {
            assert_eq!(g.apply_morphism(&morphism), g);
        }
    }
}

#[test]
fn count_automorphisms() {
    let p4 = Graph::new(4, &[(0, 1), (1, 2), (2, 3)]).canonical();
    assert_eq!(p4.automorphisms().count(), 2);
    let three_edges = Graph::new(6, &[(0, 1), (2, 3), (4, 5)]).canonical();
    assert_eq!(three_edges.automorphisms().count(), 48);
}

#[test]
fn all_automorphisms() {
    let k3 = Graph::new(3, &[(0, 1), (1, 2), (2, 0)]).canonical();
    let mut automorphisms: Vec<_> = k3.automorphisms().collect();
    automorphisms.sort_unstable();
    assert_eq!(
        automorphisms,
        vec![
            vec![0, 1, 2],
            vec![0, 2, 1],
            vec![1, 0, 2],
            vec![1, 2, 0],
            vec![2, 0, 1],
            vec![2, 1, 0]
        ]
    );
}

#[test]
fn single_automorphisms() {
    let non_symetric_graph =
        Graph::new(6, &[(0, 1), (1, 2), (2, 3), (1, 3), (3, 4), (4, 5)]).canonical();
    let automorphisms: Vec<_> = non_symetric_graph.automorphisms().collect();
    assert_eq!(automorphisms, vec![vec![0, 1, 2, 3, 4, 5]]);
}

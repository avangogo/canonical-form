extern crate canonical_form;

use canonical_form::Canonize;
use rand::{prelude::SliceRandom, rng};

/// A permutation of `0..n`, canonicalized under conjugation by `S_n`:
/// `p` acts on `pi` as `p . pi . p^-1`.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct Perm {
    image: Vec<usize>,
}

impl Perm {
    fn from_cycles(n: usize, cycles: &[&[usize]]) -> Self {
        let mut image: Vec<usize> = (0..n).collect();
        for cycle in cycles {
            for w in 0..cycle.len() {
                image[cycle[w]] = cycle[(w + 1) % cycle.len()];
            }
        }
        Perm { image }
    }
}

impl Canonize for Perm {
    fn size(&self) -> usize {
        self.image.len()
    }

    fn apply_morphism(&self, p: &[usize]) -> Self {
        let n = self.size();
        let mut image = vec![0; n];
        for i in 0..n {
            image[p[i]] = p[self.image[i]];
        }
        Perm { image }
    }

    fn invariant_neighborhood(&self, u: usize) -> Vec<Vec<usize>> {
        vec![vec![self.image[u]]]
    }
}

fn random_permutation(n: usize) -> Vec<usize> {
    let mut perm: Vec<_> = (0..n).collect();
    perm.shuffle(&mut rng());
    perm
}

fn permutations(n: usize) -> Vec<Vec<usize>> {
    fn helper(current: &mut Vec<usize>, remaining: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if remaining.is_empty() {
            out.push(current.clone());
            return;
        }
        for i in 0..remaining.len() {
            let v = remaining.remove(i);
            current.push(v);
            helper(current, remaining, out);
            current.pop();
            remaining.insert(i, v);
        }
    }
    let mut out = Vec::new();
    helper(&mut Vec::new(), &mut (0..n).collect(), &mut out);
    out
}

fn brute_force_centralizer(pi: &Perm) -> usize {
    permutations(pi.size())
        .iter()
        .filter(|p| pi.apply_morphism(p) == *pi)
        .count()
}

#[test]
fn canonical_invariant_under_conjugation() {
    for n in 0..8 {
        let cycle: Vec<usize> = (0..n).collect();
        let pi = Perm::from_cycles(n, &[&cycle]);
        let conjugate = pi.apply_morphism(&random_permutation(n));
        assert_eq!(pi.canonical(), conjugate.canonical());
    }
}

#[test]
fn morphism_to_canonical() {
    let pi = Perm::from_cycles(6, &[&[0, 1, 2, 3, 4, 5]]);
    let phi = pi.morphism_to_canonical();
    assert_eq!(pi.apply_morphism(&phi), pi.canonical());
}

#[test]
fn identity_has_full_centralizer() {
    let id = Perm::from_cycles(4, &[]).canonical();
    assert_eq!(id.automorphisms().count(), 24); // S_4
}

#[test]
fn full_cycle_centralizer() {
    let c = Perm::from_cycles(4, &[&[0, 1, 2, 3]]).canonical();
    assert_eq!(c.automorphisms().count(), 4); // <c>
}

#[test]
fn centralizer_matches_brute_force() {
    let cases = vec![
        Perm::from_cycles(4, &[&[0, 1], &[2, 3]]),
        Perm::from_cycles(4, &[&[0, 1, 2]]),
        Perm::from_cycles(5, &[&[0, 1, 2], &[3, 4]]),
    ];
    for pi in cases {
        let expected = brute_force_centralizer(&pi);
        let g = pi.canonical();
        assert_eq!(g.automorphisms().count(), expected);
    }
}

#[test]
fn empty_permutation() {
    let empty = Perm::from_cycles(0, &[]);
    assert_eq!(empty, empty.canonical());
    assert_eq!(empty.automorphisms().count(), 1);
}

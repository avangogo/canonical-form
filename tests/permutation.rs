mod common;

use canonical_form::Canonize;
use common::check_contract;

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

    fn invariant_neighborhood(&self, u: usize) -> impl Iterator<Item = (usize, u64)> {
        std::iter::once((self.image[u], 0))
    }
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
fn empty_permutation() {
    let empty = Perm::from_cycles(0, &[]);
    assert_eq!(empty, empty.canonical());
    assert_eq!(empty.automorphisms().count(), 1);
}

#[test]
fn the_contract_holds_on_small_permutations() {
    for pi in [
        Perm::from_cycles(4, &[&[0, 1], &[2, 3]]),
        Perm::from_cycles(4, &[&[0, 1, 2]]),
        Perm::from_cycles(5, &[&[0, 1, 2], &[3, 4]]),
        Perm::from_cycles(5, &[&[0, 1, 2, 3, 4]]),
    ] {
        check_contract(&pi);
    }
}

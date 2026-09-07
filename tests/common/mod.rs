//! What every implementation of `Canonize` owes, checked against brute force.
//!
//! Each test in this directory implements the trait for a different kind of
//! structure. What they all have to satisfy is checked here rather than once
//! per file.
#![allow(dead_code)] // Each test binary uses its own share of this module.

use canonical_form::Canonize;
use rand::prelude::SliceRandom;
use rand::rngs::StdRng;
use std::fmt::Debug;

/// `g` with its elements relabelled at random.
pub fn relabel<F: Canonize>(g: &F, rng: &mut StdRng) -> F {
    let mut perm: Vec<_> = (0..g.size()).collect();
    perm.shuffle(rng);
    g.apply_morphism(&perm)
}

/// Check that `rounds` random relabellings of `g` all have its canonical form.
///
/// The sampled counterpart of [`check_contract`], for structures too large to
/// enumerate the permutations of.
pub fn check_sampled<F: Canonize + Debug>(name: &str, g: &F, rounds: usize, rng: &mut StdRng) {
    let expected = g.canonical();
    for i in 0..rounds {
        assert_eq!(
            expected,
            relabel(g, rng).canonical(),
            "canonical form is not isomorphism invariant on {name} (relabelling #{i})"
        );
    }
}

/// Every permutation of `0..n`.
pub fn permutations(n: usize) -> Vec<Vec<usize>> {
    fn helper(current: &mut Vec<usize>, rest: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if rest.is_empty() {
            out.push(current.clone());
            return;
        }
        for i in 0..rest.len() {
            let v = rest.remove(i);
            current.push(v);
            helper(current, rest, out);
            let _ = current.pop();
            rest.insert(i, v);
        }
    }
    let mut out = Vec::new();
    helper(&mut Vec::new(), &mut (0..n).collect(), &mut out);
    out
}

/// Check `g` against every one of its relabellings:
///
/// - all of them have the same canonical form;
/// - `morphism_to_canonical` reaches it;
/// - `automorphisms` is exactly the permutations that fix it.
///
/// This canonizes `g.size()!` structures, so `g` must be small.
pub fn check_contract<F: Canonize + Debug>(g: &F) {
    let perms = permutations(g.size());
    let expected = g.canonical();
    for p in &perms {
        let other = g.apply_morphism(p);
        assert_eq!(other.canonical(), expected, "relabelled by {p:?}");
        let phi = other.morphism_to_canonical();
        assert_eq!(other.apply_morphism(&phi), expected, "morphism for {p:?}");
    }
    let fixing = perms
        .iter()
        .filter(|p| expected.apply_morphism(p) == expected);
    assert_eq!(
        expected.automorphisms().count(),
        fixing.count(),
        "{expected:?}"
    );
    for phi in expected.automorphisms() {
        assert_eq!(
            expected.apply_morphism(&phi),
            expected,
            "{phi:?} fixes nothing"
        );
    }
}

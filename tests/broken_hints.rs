//! What the trait promises, and what happens when it is not kept.
//!
//! This would otherwise fail somewhere inside the refinement, on an index the
//! caller has no way to connect to the method that produced it.

use canonical_form::Canonize;

/// Hints pointing outside `0..size()`.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct HintOutOfRange;

impl Canonize for HintOutOfRange {
    fn size(&self) -> usize {
        3
    }
    fn apply_morphism(&self, _p: &[usize]) -> Self {
        Self
    }
    fn invariant_neighborhood(&self, u: usize) -> impl Iterator<Item = (usize, u64)> {
        std::iter::once((u + 7, 0))
    }
}

#[test]
#[should_panic(expected = "invariant_neighborhood(0) gave a hint to 7")]
fn a_hint_outside_the_elements_says_so() {
    let _ = HintOutOfRange.canonical();
}

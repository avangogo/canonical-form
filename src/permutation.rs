//! Reusable scratch for computing with permutations.
//!
//! A permutation of `0..n` is a slice `p` with `p[u]` the image of `u`.

/// Inverts permutations into a buffer it owns, reused from call to call.
#[derive(Clone, Debug, Default)]
pub struct Inverter {
    inverse: Vec<usize>,
}

impl Inverter {
    /// The inverse of `perm`, valid until the next call. Anything but a
    /// permutation of `0..perm.len()` panics on an index.
    pub fn invert(&mut self, perm: &[usize]) -> &[usize] {
        self.inverse.clear();
        self.inverse.resize(perm.len(), 0);
        for (u, &image) in perm.iter().enumerate() {
            self.inverse[image] = u;
        }
        &self.inverse
    }
}

/// A permutation split into its cycles, taken as a function of its points so
/// that a composition needs no materializing.
#[derive(Clone, Debug, Default)]
pub struct Cycles {
    /// The representative of each point's cycle: the point the walk entered
    /// that cycle from.
    cycle: Vec<usize>,
}

impl Cycles {
    /// Split the permutation `perm` of `0..n` into its cycles, in one pass.
    pub fn load(&mut self, n: usize, perm: impl Fn(usize) -> usize) {
        self.cycle.clear();
        self.cycle.resize(n, usize::MAX);
        for start in 0..n {
            if self.cycle[start] != usize::MAX {
                continue;
            }
            let mut u = start;
            while self.cycle[u] == usize::MAX {
                self.cycle[u] = start;
                u = perm(u);
            }
            debug_assert_eq!(u, start, "`perm` is not a permutation of `0..n`");
        }
    }

    /// The point naming the cycle of `u`, so that two points have the same
    /// representative exactly when a power of the permutation relates them.
    pub fn representative(&self, u: usize) -> usize {
        self.cycle[u]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invert_is_the_inverse() {
        let mut inverter = Inverter::default();
        let perm = [3, 0, 4, 1, 2];
        let inverse = inverter.invert(&perm).to_vec();
        assert_eq!(inverse, vec![1, 3, 4, 0, 2]);
        for (u, &image) in perm.iter().enumerate() {
            assert_eq!(inverse[image], u);
            assert_eq!(perm[inverse[u]], u);
        }
        assert_eq!(inverter.invert(&inverse), perm);
        assert_eq!(inverter.invert(&[]), []);
        assert_eq!(inverter.invert(&[0, 1, 2, 3]), [0, 1, 2, 3]);
    }

    #[test]
    fn reuse_does_not_keep_the_previous_result() {
        // The buffer is shared, so a shorter permutation must not read back
        // entries left by a longer one.
        let mut inverter = Inverter::default();
        let _ = inverter.invert(&[5, 4, 3, 2, 1, 0]);
        assert_eq!(inverter.invert(&[1, 2, 0]), [2, 0, 1]);
    }

    /// The cycles of `perm`, as sorted sets of points.
    fn split(perm: &[usize]) -> Vec<Vec<usize>> {
        let mut cycles = Cycles::default();
        cycles.load(perm.len(), |u| perm[u]);
        let mut by_representative: std::collections::BTreeMap<usize, Vec<usize>> =
            Default::default();
        for u in 0..perm.len() {
            by_representative
                .entry(cycles.representative(u))
                .or_default()
                .push(u);
        }
        by_representative.into_values().collect()
    }

    #[test]
    fn cycles_of_a_known_permutation() {
        // (0 2 4)(1 3)(5)
        assert_eq!(
            split(&[2, 3, 4, 1, 0, 5]),
            vec![vec![0, 2, 4], vec![1, 3], vec![5]]
        );
        assert_eq!(split(&[0, 1, 2]), vec![vec![0], vec![1], vec![2]]);
        assert!(split(&[]).is_empty());
    }

    #[test]
    fn reuse_does_not_keep_the_previous_cycles() {
        let mut cycles = Cycles::default();
        cycles.load(6, |u| (u + 1) % 6);
        cycles.load(3, |u| u);
        assert_eq!(
            (0..3).map(|u| cycles.representative(u)).collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
    }
}

use crate::hash::mix;
use crate::{Canonize, refinement::Partition};

/// The hints of every element of a structure, colors already hashed, in one
/// flat buffer: the search builds a refiner per call, and `n + 1` allocations
/// per canonization was the first cost to show.
#[derive(Debug, Clone)]
struct Hints {
    /// All hints, concatenated
    hints: Vec<(usize, u64)>,
    /// Where each element's hints start
    offsets: Vec<usize>,
}

impl Hints {
    /// Read the hints `g` gives to each of its elements.
    fn new<F>(g: &F) -> Self
    where
        F: Canonize,
    {
        let n = g.size();
        let mut hints = Vec::with_capacity(n);
        let mut offsets = Vec::with_capacity(n + 1);
        offsets.push(0);
        for u in 0..n {
            hints.extend(
                g.invariant_neighborhood(u)
                    .map(|(v, color)| (v, mix(color))),
            );
            offsets.push(hints.len());
        }
        if let Some(i) = hints.iter().position(|&(v, _)| v >= n) {
            let u = offsets.partition_point(|&start| start <= i) - 1;
            panic!(
                "invariant_neighborhood({u}) gave a hint to {}, \
                 but the elements of this structure are 0..{n}",
                hints[i].0
            );
        }
        Self { hints, offsets }
    }

    /// The hints of `u`.
    #[inline]
    fn of(&self, u: usize) -> &[(usize, u64)] {
        &self.hints[self.offsets[u]..self.offsets[u + 1]]
    }

    /// Sieve every element reached by a hint from `part`.
    fn sieve_from(&self, part: &[usize], partition: &mut Partition) {
        for &u in part {
            for &(v, color) in self.of(u) {
                partition.sieve(v, color);
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct WL1Refiner {
    /// What the structure says reaches what.
    hints: Hints,
    /// Parts whose hints are still to be propagated. Kept between calls:
    /// `refine` runs once per node of the search tree.
    stack: Vec<usize>,
    /// Copy of the part being propagated, so that the partition can be
    /// modified while its elements are read.
    part_buffer: Vec<usize>,
}

impl WL1Refiner {
    pub fn new<F>(g: &F) -> Self
    where
        F: Canonize,
    {
        let n = g.size();
        Self {
            hints: Hints::new(g),
            stack: Vec::with_capacity(n),
            part_buffer: Vec::with_capacity(n),
        }
    }

    /// The first refinement of `partition`, splitting by what the structure
    /// says of each element on its own before propagating anything.
    ///
    /// Propagating one part at a time reaches the same fixpoint, but this fold
    /// is one contiguous pass per element and usually leaves it nothing to do.
    /// Folding from `color` is what makes it free when there is none.
    #[inline]
    pub fn refine_root(&mut self, partition: &mut Partition, color: impl Fn(usize) -> u64) {
        partition.refine_by_value(|u| {
            // Not yet a sieve value, so nothing owed to the `0` convention.
            self.hints
                .of(u)
                .iter()
                .fold(color(u), |total, &(_, c)| total.wrapping_add(c))
        });
        self.stack.clear();
        self.stack.extend(0..partition.num_parts());
        self.propagate(partition);
    }

    /// Refine `partition` again, knowing that `new_part` is the only part
    /// created since it was last refined.
    #[inline]
    pub fn refine(&mut self, partition: &mut Partition, new_part: usize) {
        self.stack.clear();
        self.stack.push(new_part);
        self.propagate(partition);
    }

    /// Propagate the parts on the stack until the partition is equitable. Each
    /// is sieved on its own and split right away, so a sieve value is the sum of
    /// the hashed colors reaching an element from that one part.
    ///
    /// `#[inline]` is for the codegen units, not the call: without it, whether
    /// this folds into the tree walk depends on how `codegen-units = 16` cuts
    /// the crate, which an unrelated edit can flip.
    #[inline]
    fn propagate(&mut self, partition: &mut Partition) {
        while let Some(part) = self.stack.pop()
            && !partition.is_discrete()
        {
            self.part_buffer.clear();
            self.part_buffer.extend_from_slice(partition.part(part));
            self.hints.sieve_from(&self.part_buffer, partition);
            partition.split(|new| {
                self.stack.push(new);
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::example::Graph;
    use crate::hash::pointed_add;

    /// Panic unless `partition` is equitable: within one part, every element is
    /// reached by the same hints from every part. It says nothing of the order
    /// the parts were propagated in.
    fn assert_equitable(refiner: &WL1Refiner, partition: &Partition) {
        partition.check_consistent();

        // The premise of that seeding: every index in `0..num_parts()` names a
        // live part, and together they cover each element exactly once.
        let mut covered = vec![false; partition.num_elems()];
        for part in 0..partition.num_parts() {
            assert!(!partition.part(part).is_empty(), "part {part} is empty");
            for &e in partition.part(part) {
                assert!(!covered[e], "element {e} is in two parts");
                covered[e] = true;
            }
        }
        assert!(covered.iter().all(|&c| c), "some element is in no part");

        let mut reaching = vec![0u64; partition.num_elems()];
        for source in 0..partition.num_parts() {
            reaching.iter_mut().for_each(|r| *r = 0);
            for &u in partition.part(source) {
                for &(v, color) in refiner.hints.of(u) {
                    reaching[v] = pointed_add(reaching[v], color);
                }
            }
            for target in 0..partition.num_parts() {
                let part = partition.part(target);
                let first = reaching[part[0]];
                for &e in part {
                    assert_eq!(
                        reaching[e], first,
                        "part {target} is split by the hints coming from part {source}"
                    );
                }
            }
        }
    }

    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
    }

    fn random_graph(n: usize, density: u64, rng: &mut Rng) -> Graph {
        let mut edges = Vec::new();
        for u in 0..n {
            for v in 0..u {
                if rng.next() % 8 < density {
                    edges.push((u, v));
                }
            }
        }
        Graph::new(n, &edges)
    }

    #[test]
    fn the_refinement_is_equitable() {
        let mut rng = Rng(0x243f_6a88_85a3_08d3);
        for n in 1..12 {
            for density in 1..8 {
                for _ in 0..20 {
                    let g = random_graph(n, density, &mut rng);
                    for sigma in 0..=n {
                        let mut partition = Partition::with_singletons(n, sigma);
                        let mut refiner = WL1Refiner::new(&g);
                        refiner.refine_root(&mut partition, |u| g.invariant_color(u));
                        assert_equitable(&refiner, &partition);
                    }
                }
            }
        }
    }

    /// The same, after individualizing a vertex: the incremental refinement
    /// from a single new part has to reach the same fixpoint as a full one.
    #[test]
    fn the_incremental_refinement_is_equitable() {
        let mut rng = Rng(0x1319_8a2e_0370_7344);
        for n in 2..10 {
            for density in 1..8 {
                for _ in 0..20 {
                    let g = random_graph(n, density, &mut rng);
                    let mut partition = Partition::with_singletons(n, 0);
                    let mut refiner = WL1Refiner::new(&g);
                    refiner.refine_root(&mut partition, |u| g.invariant_color(u));
                    // The search only individualizes an element of a part it
                    // can split, so that is what is tried here.
                    for part in 0..partition.num_parts() {
                        if partition.part(part).len() < 2 {
                            continue;
                        }
                        for i in 0..partition.part(part).len() {
                            let e = partition.part(part)[i];
                            let mut branch = partition.clone();
                            let new_part = branch.individualize(e);
                            refiner.refine(&mut branch, new_part);
                            assert_equitable(&refiner, &branch);
                        }
                    }
                }
            }
        }
    }
}

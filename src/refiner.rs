use crate::{Canonize, refinement::Partition};

#[derive(Debug, Clone)]
pub struct WL1Refiner {
    invariants: Vec<Vec<Vec<usize>>>,
    invariants_size: usize,
}

fn precompute_invariant<F>(g: &F) -> Vec<Vec<Vec<usize>>>
where
    F: Canonize,
{
    let n = g.size();
    let mut res = Vec::with_capacity(n);
    for i in 0..n {
        res.push(g.invariant_neighborhood(i));
    }
    res
}

impl WL1Refiner {
    pub fn new<F>(g: &F) -> Self
    where
        F: Canonize,
    {
        let invariants = precompute_invariant(g);
        let invariants_size = invariants.first().map_or(0, |v| v.len());
        assert!(invariants.iter().all(|v| v.len() == invariants_size));
        Self {
            invariants_size,
            invariants,
        }
    }
    pub fn dummy() -> Self {
        Self {
            invariants: Vec::new(),
            invariants_size: 0,
        }
    }
    /// Compute the coarsest refinement of `partition` with part undistinguishable
    /// by the invarriants.
    /// If `new_part` is `Some(p)`, assumes that the partition is up-to-date up to the creation
    /// of the part `p`.
    pub fn refine(&self, partition: &mut Partition, new_part: Option<usize>) {
        if !partition.is_discrete() {
            let n = partition.num_elems();
            assert!(n >= 2);
            // Stack contains the new created partitions
            let mut stack: Vec<_> = match new_part {
                Some(p) => vec![p],
                None => partition.parts().collect(),
            };
            // base
            let max_step =
                ((n + 1 - partition.num_parts()) as u64).pow(self.invariants_size as u32);
            let threshold = u64::MAX / max_step; //
            let mut part_buffer = Vec::new();
            while !stack.is_empty() && !partition.is_discrete() {
                let mut weight = 1; // multiplicator to make the values in the sieve unique
                while let Some(part) = stack.pop() {
                    part_buffer.clear();
                    part_buffer.extend_from_slice(partition.part(part));
                    let factor = (part_buffer.len() + 1) as u64;
                    #[allow(clippy::needless_range_loop)]
                    for i in 0..self.invariants_size {
                        weight *= factor;
                        // Compute sieve
                        for &u in &part_buffer {
                            for &v in &self.invariants[u][i] {
                                partition.sieve(v, weight);
                            }
                        }
                    }
                    if weight > threshold {
                        break;
                    }
                }
                partition.split(|new| {
                    stack.push(new);
                });
            }
        }
    }
}

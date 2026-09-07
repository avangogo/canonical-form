//! The individualization-refinement search.
//!
//! [`IsoTree`] walks the tree of partitions obtained by individualizing vertices
//! and refining, [`canonical_constraint`] takes the largest leaf image, and the
//! rest is the pruning that keeps the walk off the whole symmetric group.

use crate::Canonize;
use crate::children::{ChildStack, Segment};
use crate::permutation::{Cycles, Inverter};
use crate::refinement::Partition;
use crate::refiner::WL1Refiner;
use std::collections::BTreeMap;
use std::collections::btree_map::Entry::{Occupied, Vacant};

#[derive(Clone, Copy, Debug)]
struct IsoTreeNode {
    /// Number of parts, which is how the partition backtracks to this node.
    nparts: usize,
    /// The vertices that can be individualized from this node.
    children: Segment,
}

impl IsoTreeNode {
    /// Describe the node the *already refined* `partition` is at.
    fn new(partition: &Partition, children: &mut ChildStack) -> Self {
        let children = children.push(match partition.smallest_non_singleton() {
            Some(set) => partition.part(set),
            None => &[],
        });
        Self {
            nparts: partition.num_parts(),
            children,
        }
    }
}

/// A node of the branch, with the vertex labelling the edge out of it.
#[derive(Clone, Copy, Debug)]
struct Step {
    node: IsoTreeNode,
    vertex: usize,
}

/// Depth-first walk of the tree of the normalization process. One partition
/// serves it all, undone on backtrack rather than copied, and `current` is
/// `None` exactly when the walk is over.
#[derive(Clone, Debug)]
pub struct IsoTree {
    refiner: WL1Refiner,
    partition: Partition,
    current: Option<IsoTreeNode>,
    /// The branch from the root to `current`, root first.
    path: Vec<Step>,
    children: ChildStack,
    /// The automorphism relating two leaves with the same image.
    sigma: AutomorphismCycles,
}

impl IsoTree {
    /// Walk of the tree of `g` modulo the permutations that fix `0..sigma`.
    pub fn new<F: Canonize>(g: &F, sigma: usize) -> Self {
        let mut partition = Partition::with_singletons(g.size(), sigma);
        let mut refiner = WL1Refiner::new(g);
        let mut children = ChildStack::default();
        refiner.refine_root(&mut partition, |u| g.invariant_color(u));
        let root = IsoTreeNode::new(&partition, &mut children);
        Self {
            refiner,
            partition,
            current: Some(root),
            path: Vec::new(),
            children,
            sigma: AutomorphismCycles::default(),
        }
    }

    /// `true` once every node has been visited.
    pub fn is_done(&self) -> bool {
        self.current.is_none()
    }

    /// The bijection of the current node, if it is a leaf.
    pub fn leaf(&self) -> Option<&[usize]> {
        self.partition.as_bijection()
    }

    fn path_edges(&self) -> Vec<usize> {
        self.path.iter().map(|step| step.vertex).collect()
    }

    /// The depth at which the current path leaves the one with edges `other`.
    fn fork_with(&self, other: &[usize]) -> usize {
        self.path
            .iter()
            .zip(other)
            .take_while(|(step, v)| step.vertex == **v)
            .count()
    }

    /// Move to the next node in depth-first order, or out of the walk.
    pub fn advance(&mut self) {
        while let Some(mut node) = self.current.take() {
            if let Some(v) = self.children.take(&mut node.children) {
                let new_part = self.partition.individualize(v);
                self.refiner.refine(&mut self.partition, new_part);
                let child = IsoTreeNode::new(&self.partition, &mut self.children);
                self.path.push(Step { node, vertex: v });
                self.current = Some(child);
                return;
            }
            self.backtrack_from(&node); // No child left.
        }
    }

    /// Make the parent of the abandoned `node` current again, its children
    /// going with it.
    fn backtrack_from(&mut self, node: &IsoTreeNode) {
        self.children.drop_top(&node.children);
        if let Some(step) = self.path.pop() {
            self.partition.undo(step.node.nparts);
            self.current = Some(step.node);
        }
    }

    /// Abandon the current branch below `depth`.
    fn prune_to(&mut self, depth: usize) {
        while self.path.len() > depth {
            let node = self.current.take().expect("the walk is not over");
            self.backtrack_from(&node);
        }
    }

    /// Backtrack to where this branch left the branch of `old`, and drop the
    /// children that a leaf equal to `old` makes redundant.
    ///
    /// The two leaves have the same image, so their bijections differ by an
    /// automorphism `sigma = old.phi^-1 . phi`, which fixes everything
    /// individualized above the fork and so permutes the children there. One
    /// cycle of `sigma` heads subtrees with the same leaf images, so keeping one
    /// child per cycle keeps every image, hence the largest.
    fn prune_isomorphic(&mut self, old: &Leaf) {
        let fork = self.fork_with(&old.path);
        // Nothing left at the fork for the automorphism to act on.
        if self.path[fork].node.children.is_exhausted() {
            self.prune_to(fork);
            return;
        }
        let phi = self
            .partition
            .as_bijection()
            .expect("only called on a leaf");
        self.sigma.load(&old.phi, phi);
        self.prune_to(fork);
        if let Some(node) = &mut self.current {
            // A taken child covers its cycle, explored or pruned; a dropped
            // one does not, hence reading only the taken.
            for &taken in self.children.taken(&node.children) {
                let _ = self.sigma.claim_cycle(taken);
            }
            self.children
                .retain_untaken(&mut node.children, |v| self.sigma.claim_cycle(v));
        }
    }
}

/// The cycles of the automorphism relating two leaves with the same image, and
/// which are accounted for. It is never materialized, only walked.
#[derive(Clone, Debug, Default)]
struct AutomorphismCycles {
    /// Scratch for building `sigma` during [`Self::load`].
    inverter: Inverter,
    cycles: Cycles,
    /// Whether a cycle, indexed by its representative, has been claimed.
    claimed: Vec<bool>,
}

impl AutomorphismCycles {
    /// Load `sigma = left^-1 . right` and split it into its cycles.
    fn load(&mut self, left: &[usize], right: &[usize]) {
        let n = right.len();
        let inverse = self.inverter.invert(left);
        self.cycles.load(n, |point| inverse[right[point]]);
        self.claimed.clear();
        self.claimed.resize(n, false);
    }

    /// Claim the cycle of `point`, returning whether it was still free.
    fn claim_cycle(&mut self, point: usize) -> bool {
        !std::mem::replace(&mut self.claimed[self.cycles.representative(point)], true)
    }
}

/// A leaf already visited, held in `zeta` under the image it produced.
struct Leaf {
    /// The edges of the path that reached the leaf.
    path: Vec<usize>,
    phi: Vec<usize>,
}

/// Normal form of `g` under the permutations fixing `0..sigma`, with a morphism
/// reaching it: the bijection of the first leaf that produced it.
pub fn canonical_constraint<F>(g: &F, sigma: usize) -> (F, Vec<usize>)
where
    F: Canonize,
{
    // The images of `g` already computed, under the leaf producing each.
    let mut zeta: BTreeMap<F, Leaf> = BTreeMap::new();
    let mut tree = IsoTree::new(g, sigma);
    while !tree.is_done() {
        if let Some(phi) = tree.leaf() {
            let image = g.apply_morphism(phi);
            match zeta.entry(image) {
                Occupied(entry) => {
                    // This branch is a copy of one already explored.
                    tree.prune_isomorphic(entry.get());
                }
                Vacant(entry) => {
                    let _ = entry.insert(Leaf {
                        path: tree.path_edges(),
                        phi: phi.to_vec(),
                    });
                }
            }
        }
        tree.advance();
    }
    let (g_max, leaf) = zeta.into_iter().next_back().unwrap();
    (g_max, leaf.phi)
}

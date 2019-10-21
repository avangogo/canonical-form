//! Algorithm to reduce combinatorial structures modulo isomorphism.
//!
//! This can typically be used to to test if two graphs are isomorphic.
//!
//! The algorithm manipulates its input as a black box by
//! the action of permutations
//! and by testing equallity with element of its orbit,
//! plus some user-defined functions
//! that help to break symmetries.
//!
//!```
//!use canonical_form::Canonize;
//!
//!// Simple Graph implementation as adjacency lists
//!#[derive(Ord, PartialOrd, PartialEq, Eq, Clone, Debug)]
//!struct Graph {
//!       adj: Vec<Vec<usize>>,
//!}
//!
//!
//!impl Graph {
//!    fn new(n: usize, edges: &[(usize, usize)]) -> Self {
//!        let mut adj = vec![Vec::new(); n];
//!        for &(u, v) in edges {
//!            adj[u].push(v);
//!            adj[v].push(u);
//!        }
//!        for list in &mut adj {
//!            list.sort() // Necessary to make the derived `==` correct
//!        }
//!        Graph { adj }
//!    }
//!}
//!
//!// The Canonize trait allows to use the canonial form algorithms
//!impl Canonize for Graph {
//!    fn size(&self) -> usize {
//!        self.adj.len()
//!    }
//!    fn apply_morphism(&self, perm: &[usize]) -> Self {
//!        let mut adj = vec![Vec::new(); self.size()];
//!        for (i, nbrs) in self.adj.iter().enumerate() {
//!            adj[perm[i]] = nbrs.iter().map(|&u| perm[u]).collect();
//!            adj[perm[i]].sort();
//!        }
//!        Graph { adj }
//!    }
//!    fn invariant_neighborhood(&self, u: usize) -> Vec<Vec<usize>> {
//!        vec![self.adj[u].clone()]
//!    }
//!}
//!
//! // Usage of library functions
//! // Two isomorphic graphs
//! let c5 = Graph::new(5, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]);
//! let other_c5 = Graph::new(5, &[(0, 2), (2, 1), (1, 4), (4, 3), (3, 0)]);
//! assert_eq!(c5.canonical(), other_c5.canonical());
//!
//! // Non-isomorphic graphs
//! let p5 = Graph::new(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
//! assert!(c5.canonical() != p5.canonical());
//!
//! // Recovering the permutation that gives the canonical form
//! let p = c5.morphism_to_canonical();
//! assert_eq!(c5.apply_morphism(&p), c5.canonical());
//!
//! // Enumerating automorphisms
//! assert_eq!(c5.canonical().automorphisms().count(), 10)
//!```

#![warn(
    missing_docs,
    missing_debug_implementations,
    missing_copy_implementations,
    trivial_casts,
    trivial_numeric_casts,
    unsafe_code,
    unstable_features,
    unused_import_braces,
    unused_qualifications,
    unused_labels,
    unused_results
)]

mod refine;

use crate::refine::Partition;
use std::collections::btree_map::Entry::{Occupied, Vacant};
use std::collections::BTreeMap;
use std::rc::Rc;
pub mod example;

/// Objects that can be reduced modulo the actions of a permutation group.
///
/// An object that implement this trait has a set of elements assimilated to
/// {0,...,n-1} on which the group of permutations can act.
/// The purpose of the trait is to compute a normal form of
/// the object modulo the permutation of its elements.
pub trait Canonize
where
    Self: Sized + Ord + Clone,
{
    /// Return the number of vertices.
    ///
    /// The elements of `self` are assimilated to the number of `0..self.size()`.
    fn size(&self) -> usize;

    /// Return the result of the action of a permuation `p` on the object.
    ///
    /// The input `p` is guarenteed to be a permutation represented
    /// as a slice of size `self.size()`
    /// where `p[u]` is the image of `u` by the permutation.
    ///
    /// The action of the permutation set on `Self` is assumed to be a group action.
    /// ```
    /// use canonical_form::Canonize;
    /// use canonical_form::example::Graph;
    ///
    /// let g = Graph::new(4, &[(1, 2), (2, 3)]);
    ///
    /// let p = &[1, 2, 0, 3];
    /// assert_eq!(g.apply_morphism(p), Graph::new(4, &[(2, 0), (0, 3)]));
    ///
    /// let identity = &[0, 1, 2, 3];
    /// assert_eq!(g.apply_morphism(identity), g);
    ///
    /// let q = &[1, 0, 3, 2];
    /// let p_circle_q = &[2, 1, 3, 0];
    /// assert_eq!(
    ///     g.apply_morphism(q).apply_morphism(p),
    ///     g.apply_morphism(p_circle_q)
    /// )
    /// ```
    fn apply_morphism(&self, p: &[usize]) -> Self;

    /// Optionally returns a value for each node that is invariant by isomorphism.
    ///
    /// If defined, the returned vector `c` must have size `self.len()`, where `c[u]`
    /// is the value associated to the element `u`. It must satisfy the property that if
    /// `c[u]` and `c[v]` are different then no automorphism of `self`
    /// maps `u` to `v`.
    fn invariant_coloring(&self) -> Option<Vec<u64>> {
        None
    }

    /// Return lists of vertices that are invariant isomorphism.
    ///
    /// This function helps the algorithm to be efficient.
    /// The output `invar` is such that each `invar[i]` is a vector
    /// of vertices `[v1, ..., vk]`
    /// (so the `vi`s are elements of `0..self.size()`) such that
    /// for every permutation `p`,
    /// `self.apply_morphism(p).invariant_neighborhood(p[u])[i]`
    /// is equal to `[p[v1], ..., p[vk]]` up to reordering.
    ///
    /// The length of the output (the number of lists) has to be independent from `u`.
    fn invariant_neighborhood(&self, _u: usize) -> Vec<Vec<usize>> {
        Vec::new()
    }

    /// Computes a canonical form of a combinatorial object.
    ///
    /// This is the main function provided by this trait.
    /// A canonical form is a function that assigns to an object `g` (e.g. a graph)
    /// another object of sane type `g.canonical()` that is isomorphic to `g`
    /// with the property that `g1` and `g2` are isomorphic if and only if
    /// `g1.canocial() == g2.canonical()`.
    /// ```
    /// use canonical_form::Canonize;
    /// use canonical_form::example::Graph;
    ///
    /// let p5 = Graph::new(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
    /// let other_p5 = Graph::new(5, &[(3, 4), (0, 4), (0, 2), (1, 2)]);
    /// assert_eq!(p5.canonical(), other_p5.canonical());
    ///
    /// let not_p5 = Graph::new(5, &[(0, 1), (1, 2), (2, 0), (3, 4)]);
    /// assert!(p5.canonical() != not_p5.canonical());
    /// ```
    fn canonical(&self) -> Self {
        self.canonical_typed(0)
    }

    /// The "typed" objects refers to the case where only
    /// the action of permutations that are constant
    /// on `0..sigma` are considered.
    ///
    /// So `g.canonical_typed(sigma)` returns a normal form of `g`
    /// modulo the permutations that stabilize the `sigma` first vertices.
    /// ```
    /// use canonical_form::Canonize;
    /// use canonical_form::example::Graph;
    ///
    /// // p5 with the edge (0, 1) at an end
    /// let p5 = Graph::new(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
    /// // p5 with the edge (0, 1) is in the middle
    /// let p5_bis = Graph::new(5, &[(3, 4), (4, 1), (1, 0), (0, 2)]);
    /// // If we fix the vertices 0 and 1, these two (rooted) graphs are different
    /// assert!(p5.canonical_typed(2) != p5_bis.canonical_typed(2));
    ///
    /// let p5_ter = Graph::new(5, &[(0, 1), (1, 3), (3, 4), (4, 2)]);
    /// assert_eq!(p5.canonical_typed(2), p5_ter.canonical_typed(2));
    /// ```
    fn canonical_typed(&self, sigma: usize) -> Self {
        let partition = Partition::with_singletons(self.size(), sigma);
        canonical_constraint(self, partition)
    }

    #[inline]
    /// Return a permutation `p` such that `self.apply_morphism(&p) = self.canonical()`.
    /// ```
    /// use canonical_form::Canonize;
    /// use canonical_form::example::Graph;
    ///
    /// let g = Graph::new(6, &[(0, 1), (1, 2), (2, 3), (3, 4), (3, 5)]);
    /// let p = g.morphism_to_canonical();
    /// assert_eq!(g.apply_morphism(&p), g.canonical());
    /// ```
    fn morphism_to_canonical(&self) -> Vec<usize> {
        self.morphism_to_canonical_typed(0)
    }

    /// Return a permutation `phi` such that
    /// `g.apply_morphism(&phi) = canonical_typed(&g, sigma)`.
    /// ```
    /// use canonical_form::Canonize;
    /// use canonical_form::example::Graph;
    ///
    /// let g = Graph::new(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
    /// let p = g.morphism_to_canonical_typed(2);
    /// assert_eq!(g.apply_morphism(&p), g.canonical_typed(2));
    /// ```
    fn morphism_to_canonical_typed(&self, sigma: usize) -> Vec<usize> {
        assert!(sigma <= self.size());
        let partition = Partition::with_singletons(self.size(), sigma);
        morphism_to_canonical_constraint(self, partition)
    }

    /// Iterator on the automorphism group of `g`.
    ///
    /// The input `g` must be in normal form.
    /// ```
    /// use canonical_form::Canonize;
    /// use canonical_form::example::Graph;
    ///
    /// let c6 = Graph::new(6, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]).canonical();
    ///
    /// let mut count = 0;
    /// for p in c6.automorphisms() {
    ///     assert_eq!(c6.apply_morphism(&p), c6);
    ///     count += 1;
    /// }
    /// assert_eq!(count, 12);
    /// ```
    #[inline]
    fn automorphisms(&self) -> AutomorphismIterator<Self> {
        self.stabilizer(0)
    }

    /// Iterator on the automorphisms of `g`
    /// that fix the `sigma` first vertices.
    ///
    /// The input `g` must be in normal form computed with `canonical_typed`.
    /// ```
    /// use canonical_form::Canonize;
    /// use canonical_form::example::Graph;
    ///
    /// // Cube graph with one fixed vertex
    /// let cube = Graph::new(8, &[(0, 1), (1, 2), (2, 3), (3, 0),
    ///                            (4, 5), (5, 6), (6, 7), (7, 4),
    ///                            (0, 4), (1, 5), (2, 6), (3, 7)]).canonical_typed(1);
    ///
    /// let mut count = 0;
    /// for p in cube.stabilizer(1) {
    ///     assert_eq!(cube.apply_morphism(&p), cube);
    ///     assert_eq!(p[0], 0);
    ///     count += 1;
    /// }
    /// assert_eq!(count, 6);
    /// ```
    #[inline]
    fn stabilizer(&self, sigma: usize) -> AutomorphismIterator<Self> {
        let mut partition = Partition::simple(self.size());
        for i in 0..sigma {
            let _ = partition.individualize(i);
        }
        AutomorphismIterator::new(self, partition)
    }
}

/// Return the next part to be refined.
/// This part is chosen as a smallest part with at least 2 elements.
/// Return None is the partition is discrete.
fn target_selector(part: &Partition) -> Option<usize> {
    let mut min = usize::max_value();
    let mut arg_min = None;
    for i in part.parts() {
        let length = part.part(i).len();
        if 2 <= length && (length < min) {
            min = length;
            arg_min = Some(i);
        }
    }
    arg_min
}

fn precompute_invariant<F>(g: &F) -> Vec<Vec<Vec<usize>>>
where
    F: Canonize,
{
    let n = g.size();
    let mut res = Vec::with_capacity(n);
    for i in 0..n {
        res.push(g.invariant_neighborhood(i))
    }
    res
}

/// Compute the coarsest refinement of `partition` with part undistinguishable
/// by the invarriants.
/// If `new_part` is `Some(p)`, assumes that the partition is up-to-date up to the creation
/// of the part `p`.
fn refine(partition: &mut Partition, invariants: &[Vec<Vec<usize>>], new_part: Option<usize>) {
    if !partition.is_discrete() {
        let n = partition.num_elems();
        assert!(n >= 2);
        let invariant_size = invariants[0].len();
        debug_assert!(invariants.iter().all(|v| v.len() == invariant_size));
        // Stack contains the new created partitions
        let mut stack: Vec<_> = match new_part {
            Some(p) => vec![p],
            None => partition.parts().collect(),
        };
        // base
        let max_step = ((n + 1 - partition.num_parts()) as u64).pow(invariant_size as u32);
        let threshold = u64::max_value() / max_step; //
        let mut part_buffer = Vec::new();
        while !stack.is_empty() && !partition.is_discrete() {
            let mut weight = 1; // multiplicator to make the values in the sieve unique
            while let Some(part) = stack.pop() {
                part_buffer.clear();
                part_buffer.extend_from_slice(partition.part(part));
                let factor = (part_buffer.len() + 1) as u64;
                for i in 0..invariant_size {
                    weight *= factor;
                    // Compute sieve
                    for &u in &part_buffer {
                        for &v in &invariants[u][i] {
                            partition.sieve(v, weight)
                        }
                    }
                }
                if weight > threshold {
                    break;
                };
            }
            partition.split(|new| {
                stack.push(new);
            })
        }
    }
}

/// Return the first index on which `u` and `v` differ.
fn fca(u: &[usize], v: &[usize]) -> usize {
    let mut i = 0;
    while i < u.len() && i < v.len() && u[i] == v[i] {
        i += 1;
    }
    i
}

/// Node of the tree of the normalization process
#[derive(Clone, Debug)]
struct IsoTreeNode {
    nparts: usize,
    children: Vec<usize>,
    inv: Rc<Vec<Vec<Vec<usize>>>>,
}

impl IsoTreeNode {
    fn root<F: Canonize>(partition: &mut Partition, g: &F) -> Self {
        let inv = Rc::new(precompute_invariant(g));
        if let Some(coloring) = g.invariant_coloring() {
            partition.refine_by_value(&coloring, |_| {})
        }
        Self::new(partition, inv, None)
    }
    fn new(
        partition: &mut Partition,
        inv: Rc<Vec<Vec<Vec<usize>>>>,
        new_part: Option<usize>,
    ) -> Self {
        refine(partition, &inv, new_part);
        Self {
            children: match target_selector(partition) {
                Some(set) => partition.part(set).to_vec(),
                None => Vec::new(),
            },
            nparts: partition.num_parts(),
            inv,
        }
    }
    fn explore(&self, v: usize, pi: &mut Partition) -> Self {
        debug_assert!(self.is_restored(pi));
        let new_part = pi.individualize(v);
        Self::new(pi, self.inv.clone(), new_part)
    }
    // Should never be used
    fn dummy() -> Self {
        Self {
            children: Vec::new(),
            nparts: 1,
            inv: Rc::new(Vec::new()),
        }
    }
    fn restore(&self, partition: &mut Partition) {
        partition.undo(self.nparts)
    }
    fn is_restored(&self, partition: &mut Partition) -> bool {
        partition.num_parts() == self.nparts
    }
}

/// Normal form of `g` under the action of isomorphisms that
/// stabilize the parts of `partition`.
fn canonical_constraint<F>(g: &F, mut partition: Partition) -> F
where
    F: Canonize,
{
    // contains the images of `g` already computed associated to the path to the corresponding leaf
    let mut zeta: BTreeMap<F, Vec<usize>> = BTreeMap::new();
    let mut tree = Vec::new(); // A stack of IsoTreeNode
    let mut path = Vec::new(); // Current path as a vector of chosen vertices
    let mut node = IsoTreeNode::root(&mut partition, g);
    loop {
        // If we have a leaf, treat it
        if let Some(phi) = partition.as_bijection() {
            match zeta.entry(g.apply_morphism(phi)) {
                Occupied(entry) =>
                // We are in a branch isomorphic to a branch we explored
                {
                    let k = fca(entry.get(), &path) + 1;
                    tree.truncate(k);
                    path.truncate(k);
                }
                Vacant(entry) => {
                    let _ = entry.insert(path.clone());
                }
            }
        };
        // If there is a child, explore it
        if let Some(u) = node.children.pop() {
            let new_node = node.explore(u, &mut partition);
            tree.push(node);
            path.push(u);
            node = new_node;
        } else {
            match tree.pop() {
                Some(n) => {
                    node = n;
                    let _ = path.pop();
                    node.restore(&mut partition); // backtrack the partition
                }
                None => break,
            }
        };
    }
    let (g_max, _) = zeta.into_iter().next_back().unwrap(); // return the largest image found
    g_max
}

/// Iterator on the automorphisms of a combinatorial structure.
#[derive(Clone, Debug)]
pub struct AutomorphismIterator<F> {
    tree: Vec<IsoTreeNode>,
    node: IsoTreeNode,
    partition: Partition,
    g: F,
}

impl<F: Canonize> AutomorphismIterator<F> {
    /// Iterator on the automorphisms of `g` that preserve `partition`.
    fn new(g: &F, mut partition: Partition) -> Self {
        debug_assert!(g == &canonical_constraint(g, partition.clone()));
        Self {
            tree: vec![IsoTreeNode::root(&mut partition, g)],
            partition,
            node: IsoTreeNode::dummy(), // Dummy node that will be unstacked at the first iteration
            g: g.clone(),
        }
    }
}

impl<F: Canonize> Iterator for AutomorphismIterator<F> {
    type Item = Vec<usize>;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(u) = self.node.children.pop() {
                let new_node = self.node.explore(u, &mut self.partition);
                let old_node = std::mem::replace(&mut self.node, new_node);
                self.tree.push(old_node)
            } else {
                match self.tree.pop() {
                    Some(n) => {
                        n.restore(&mut self.partition);
                        self.node = n
                    }
                    None => return None,
                }
            }
            if let Some(phi) = self.partition.as_bijection() {
                if self.g.apply_morphism(phi) == self.g {
                    return Some(phi.to_vec());
                }
            };
        }
    }
}

/// Return a morphism `phi`
/// such that `g.apply_morphism(phi) = canonical_constraint(g, partition)`.
fn morphism_to_canonical_constraint<F>(g: &F, mut partition: Partition) -> Vec<usize>
where
    F: Canonize,
{
    // initialisation
    let mut tree = Vec::new();
    let mut node = IsoTreeNode::root(&mut partition, g);
    let mut max = None;
    let mut phimax = Vec::new();
    loop {
        if let Some(phi) = partition.as_bijection() {
            // If node is a leaf
            let phi_g = Some(g.apply_morphism(phi));
            if phi_g > max {
                max = phi_g;
                phimax = phi.to_vec();
            }
        };
        if let Some(u) = node.children.pop() {
            let new_node = node.explore(u, &mut partition);
            tree.push(node);
            node = new_node;
        } else {
            match tree.pop() {
                Some(n) => {
                    n.restore(&mut partition);
                    node = n
                }
                None => break,
            }
        }
    }
    phimax
}

/// Tests
#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Ord, PartialOrd, PartialEq, Eq, Clone, Debug)]
    struct Graph {
        adj: Vec<Vec<usize>>,
    }

    impl Graph {
        fn new(n: usize, edges: &[(usize, usize)]) -> Self {
            let mut adj = vec![Vec::new(); n];
            for &(u, v) in edges {
                adj[u].push(v);
                adj[v].push(u);
            }
            Graph { adj }
        }
    }

    impl Canonize for Graph {
        fn size(&self) -> usize {
            self.adj.len()
        }
        fn apply_morphism(&self, perm: &[usize]) -> Self {
            let mut adj = vec![Vec::new(); self.size()];
            for (i, nbrs) in self.adj.iter().enumerate() {
                adj[perm[i]] = nbrs.iter().map(|&u| perm[u]).collect();
                adj[perm[i]].sort();
            }
            Graph { adj }
        }
        fn invariant_neighborhood(&self, u: usize) -> Vec<Vec<usize>> {
            vec![self.adj[u].clone()]
        }
    }

    #[test]
    fn graph() {
        let c5 = Graph::new(5, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]);
        let other_c5 = Graph::new(5, &[(0, 2), (2, 1), (1, 4), (4, 3), (3, 0)]);
        assert_eq!(c5.canonical(), other_c5.canonical());

        let p5 = Graph::new(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
        assert!(c5.canonical() != p5.canonical());

        let p = c5.morphism_to_canonical();
        assert_eq!(c5.apply_morphism(&p), c5.canonical());
    }
    #[test]
    fn automorphisms_iterator() {
        let c4 = Graph::new(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]).canonical();
        let mut count = 0;
        for phi in c4.automorphisms() {
            assert_eq!(c4.apply_morphism(&phi), c4);
            count += 1;
        }
        assert_eq!(count, 8)
    }
}

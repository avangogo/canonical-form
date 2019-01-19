//! Algorithms to reduce combinatorial structures modulo isomorphism.
//!
//! This can typically be use to to test if two graphs are isomorphic.

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
)]

mod refine;

use crate::refine::Partition;
use std::collections::BTreeMap;

/// Objects that can be reduced modulo the actions of a permutation group.
///
/// An object implement this trait if it has elements
/// and if the group of permutations on of this elements
/// acts on the object.
pub trait Canonize
where
    Self: Sized + Ord + Clone,
{
    /// Returns the number of vertices.
    ///
    /// The elements of `x` are assimilated to the number of `0..x.self()`.
    fn size(&self) -> usize;
    
    /// Returns the result of the action of a permuation `perm` on the object.
    ///
    /// The permutation `perm` is represented as a slice of size `self.size()`
    /// where `perm[u]` is the image of `u` by the permutation.
    fn apply_morphism(&self, perm: &[usize]) -> Self;

    /// Can return a coloring that is invariant by isomorphism.
    ///
    /// This coloring is expressed as a vector `c` such that `c[u]`
    /// is the color of `u`. It need to satisfy the property that if
    /// `c[u]` and `c[v]` are different then no automorphism of `self`
    /// maps `u` to `v`.
    fn invariant_coloring(&self) -> Option<Vec<u64>> {
        None
    }

    /// Return lists of vertices that are invariant isomorphism.
    ///
    /// The output `inv` is a vector such that each `inv[i]` is a vector
    /// of distinct vertices `[v1, ..., vk]`
    /// (so the `vi`  elements of `0..self.size()`) such that
    /// for every permutation `perm`,
    /// `self.invariant_neighborhood(perm[u])[i]`
    /// is equal to `[perm[v1], ..., perm[vk]]` up to reordering.
    ///
    /// The length of the output (the number of lists) has to be independent
    /// of `u`.
    fn invariant_neighborhood(&self, _u: usize) -> Vec<Vec<usize>> {
        Vec::new()
    }
}

/// Return the next part to be refined.
/// This part is chosen as a smallest part with at least 2 elements.
/// Return None is the partition is discrete.
fn target_selector(part: &Partition) -> Option<usize> {
    let mut min = std::usize::MAX;
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

/// Apply to `g` the permutation that correspond to the partitions of `part`.
fn apply<F>(part: &Partition, g: &F) -> F
where
    F: Canonize,
{
    let morphism = part.to_bijection().unwrap();
    g.apply_morphism(&morphism)
}

/// Compute the coarsest refinement of `partition` with part undistinguishable
/// by the invarriants.
fn refine<F>(partition: &mut Partition, g: &F)
where
    F: Canonize,
{
    if !partition.is_discrete() {
        let n = g.size();
        assert!(n >= 2);
        // Apply coloring if any
        if let Some(coloring) = g.invariant_coloring() {
            partition.refine_by_value(&coloring, |_| {})
        }
        // Precompute invariants
        let mut invariants = Vec::with_capacity(n);
        for i in 0..n {
            invariants.push(g.invariant_neighborhood(i))
        }
        let invariant_size = invariants[0].len();
        assert!(invariants.iter().all(|v| v.len() == invariant_size));
        let mut crible = vec![0; n];
        // Stack contains the new created partitions
        let mut stack: Vec<_> = partition.parts();
        // base
        let m = (n + 1) as u64;
        let step = m.pow(invariant_size as u32);
        let threshold = std::u64::MAX / step;
        while !stack.is_empty() && !partition.is_discrete() {
            // Re-initialize crible
            for e in &mut crible {
                *e = 0;
            }
            let mut weight = 1;
            while let Some(part) = stack.pop() {
                for i in 0..invariant_size {
                    weight *= m;
                    // Compute crible and new_part
                    for &u in partition.part(part) {
                        for &v in &invariants[u][i] {
                            crible[v] += weight;
                        }
                    }
                }
                if weight > threshold {
                    break;
                };
            }
            partition.refine_by_value(&crible, |new| {
                stack.push(new);
            });
        }
    }
}

// Return the first index on which `u` and `v` differ.
fn fca(u: &[usize], v: &[usize]) -> usize {
    let mut i = 0;
    while i < u.len() && i < v.len() && u[i] == v[i] {
        i += 1;
    }
    i
}

/// Node of the tree of the normalization process
#[derive(Clone,Debug)]
struct IsoTreeNode {
    pi: Partition,
    children: Vec<usize>,
}

impl IsoTreeNode {
    fn new<F: Canonize>(mut partition: Partition, g: &F) -> Self {
        refine(&mut partition, g);
        IsoTreeNode {
            children: get_children(&partition),
            pi: partition,
        }
    }
    fn explore<F: Canonize>(&self, v: usize, g: &F) -> Self {
        let mut new_pi = self.pi.clone();
        new_pi.individualize(v);
        refine(&mut new_pi, g);
        Self::new(new_pi, g)
    }
    fn is_leaf(&mut self) -> bool {
        self.pi.is_discrete()
    }
    fn empty() -> Self {
        IsoTreeNode {
            children: Vec::new(),
            pi: Partition::simple(0),
        }
    }
}

fn get_children(pi: &Partition) -> Vec<usize> {
    match target_selector(pi) {
        Some(p) => pi.part(p).to_vec(),
        None => Vec::new(),
    }
}

/// Normal form of `g` under the action of isomorphisms that
/// stabilize the parts of `partition`.
fn canonical_form_constraint<F>(g: &F, partition: Partition) -> F
    where
    F: Canonize,
{
    let mut zeta: BTreeMap<F, Vec<usize>> = BTreeMap::new();
    // initialisation
    let mut tree = Vec::new();
    let mut path = Vec::new();
    let mut node = IsoTreeNode::new(partition, g);
    loop {
        if node.is_leaf() {
            let gv = apply(&node.pi, g);
            if let Some(u) = zeta.get(&gv) {
                let k = fca(u, &path) + 1;
                tree.truncate(k);
                path.truncate(k);
            } else {
                zeta.insert(gv, path.clone());
            }
        };
        if let Some(u) = node.children.pop() {
            path.push(u);
            let new_node = node.explore(u, g);
            tree.push(node);
            node = new_node;
        } else {
            match tree.pop() {
                Some(n) => {
                    node = n;
                    path.pop();
                }
                None => break,
            }
        }
    }
    zeta.keys().next_back().unwrap().clone()
}

/// Iterator on the automorphisms of a combinatorial structure.
#[derive(Clone,Debug)]
pub struct AutomorphismIterator<F> {
    tree: Vec<IsoTreeNode>,
    node: IsoTreeNode,
    g: F,
}

impl<F: Canonize> AutomorphismIterator<F> {
    /// Iterator on the automorphisms of `g` that preserve `partition`.
    fn new(g: &F, partition: Partition) -> Self {
        assert!(*g == canonical_form_constraint(g, partition.clone()));
        AutomorphismIterator {
            tree: vec![IsoTreeNode::new(partition, g)],
            node: IsoTreeNode::empty(), // Dummy node that will be unstacked at the first iteration
            g: g.clone(),
        }
    }
    /// Create the iterator on the automorphisms of `g`
    /// rooted on the `sigma` first vertices.
    pub fn typed(g: &F, sigma: usize) -> Self {
        let mut partition = Partition::simple(g.size());
        for i in 0..sigma {
            partition.individualize(i)
        }
        Self::new(g, partition)
    }
}

impl<F: Canonize> Iterator for AutomorphismIterator<F> {
    type Item = Vec<usize>;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(u) = self.node.children.pop() {
                let new_node = self.node.explore(u, &self.g);
                let old_node = std::mem::replace(&mut self.node, new_node);
                self.tree.push(old_node)
            } else {
                match self.tree.pop() {
                    Some(n) => self.node = n,
                    None => return None,
                }
            }
            if self.node.is_leaf() {
                let morphism = self.node.pi.to_bijection().unwrap();
                let gv = self.g.apply_morphism(&morphism);
                if gv == self.g {
                    return Some(morphism);
                }
            };
        }
    }
}

/// Normal form of `g` modulo the isomorphisms that
/// stabilize the sigma first vertices.
pub fn canonical_form_typed<F>(g: &F, sigma: usize) -> F
where
    F: Canonize,
{
    let mut partition = Partition::simple(g.size());
    for i in 0..sigma {
        partition.individualize(i)
    }
    canonical_form_constraint(&g, partition)
}

/// Computes a normal form of a combinatorial object.
///
/// A normal form is a function that assigns to an object `g` (e.g. a graph)
/// an object `canonical_form(g)` that is isomorphic to `g`
/// such that if `g1` and `g2` are isomorphic then
/// `canonical_form(g)=canonical_form(g)`.
/// The exact specification of this function does not matter.
pub fn canonical_form<F>(g: &F) -> F
where
    F: Canonize,
{
    canonical_form_typed(g, 0)
}

/// Return a morphism `phi`
/// such that `g.apply_morphism(phi) = canonical_form_constraint(g, partition)`.
fn canonical_form_morphism_constraint<F>(g: &F, partition: Partition) -> Vec<usize>
where
    F: Canonize,
{
    // initialisation
    let mut tree = Vec::new();
    let mut node = IsoTreeNode::new(partition, g);
    let mut max = None;
    let mut phimax = Vec::new();
    loop {
        if node.is_leaf() {
            let phi = node.pi.to_bijection().unwrap();
            let gv = Some(g.apply_morphism(&phi));
            if gv > max {
                max = gv;
                phimax = phi;
            }
        };
        if let Some(u) = node.children.pop() {
            let new_node = node.explore(u, g);
            tree.push(node);
            node = new_node;
        } else {
            match tree.pop() {
                Some(n) => node = n,
                None => break,
            }
        }
    }
    phimax
}

/// Return a morphism `phi` such that
/// `g.apply_morphism(phi) = canonical_form_typed(g, sigma)`.
pub fn canonical_form_typed_morphism<F>(g: &F, sigma: usize) -> Vec<usize>
where
    F: Canonize,
{
    assert!(sigma <= g.size());
    let mut partition = Partition::simple(g.size());
    for v in 0..sigma {
        partition.individualize(v);
    }
    canonical_form_morphism_constraint(g, partition)
}

/// Return a morphism `phi` such that `g.apply_morphism(phi) = canonical_form(g)`.
pub fn canonical_form_morphism<F>(g: &F) -> Vec<usize>
where
    F: Canonize,
{
    canonical_form_typed_morphism(g, 0)
}

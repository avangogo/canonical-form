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
//!    fn invariant_neighborhood(&self, u: usize) -> impl Iterator<Item = (usize, u64)> {
//!        self.adj[u].iter().map(|&v| (v, 0))
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

mod children;
mod hash;
mod permutation;
mod refinement;
mod refiner;
mod search;

use crate::search::{IsoTree, canonical_constraint};
pub mod example;

/// Objects that can be reduced modulo the actions of a permutation group.
///
/// An object that implements this trait has a set of elements assimilated to
/// {0,...,n-1} on which the group of permutations can act.
/// The purpose of the trait is to compute a normal form of
/// the object modulo the permutation of its elements.
pub trait Canonize
where
    Self: Sized + Ord + Clone,
{
    /// Return the number of vertices.
    ///
    /// The elements of `self` are assimilated to the numbers `0..self.size()`.
    fn size(&self) -> usize;

    /// Return the result of the action of a permutation `p` on the object.
    ///
    /// The input `p` is guaranteed to be a permutation represented
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

    /// Returns a value for the element `u` that is invariant by isomorphism.
    ///
    /// The value must satisfy the property that if `self.invariant_color(u)`
    /// and `self.invariant_color(v)` are different then no automorphism of
    /// `self` maps `u` to `v`.
    ///
    /// The default gives the same value to every element, which costs nothing.
    fn invariant_color(&self, _u: usize) -> u64 {
        0
    }

    /// Return colored hints on the neighborhood of `u` that are invariant by
    /// isomorphism.
    ///
    /// This function helps the algorithm to be efficient.
    /// The output is a sequence of pairs `(v, color)` where `v` is an element of
    /// `0..self.size()`. For every permutation `p`,
    /// `self.apply_morphism(p).invariant_neighborhood(p[u])` must be equal to
    /// `self.invariant_neighborhood(u)` with every `v` replaced by `p[v]`,
    /// up to reordering.
    ///
    /// A color is an opaque label: only equality between colors is used, and
    /// their numeric values carry neither order nor magnitude. Two hints of
    /// color `21` are never confused with one hint of color `42`.
    fn invariant_neighborhood(&self, _u: usize) -> impl Iterator<Item = (usize, u64)> {
        std::iter::empty()
    }

    /// Computes a canonical form of a combinatorial object.
    ///
    /// This is the main function provided by this trait.
    /// A canonical form is a function that assigns to an object `g` (e.g. a graph)
    /// another object of same type `g.canonical()` that is isomorphic to `g`
    /// with the property that `g1` and `g2` are isomorphic if and only if
    /// `g1.canonical() == g2.canonical()`.
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

    /// The "typed" objects refer to the case where only
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
        canonical_constraint(self, sigma).0
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
        canonical_constraint(self, sigma).1
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
        AutomorphismIterator::new(self, sigma)
    }
}

/// Iterator on the automorphisms of a combinatorial structure.
#[derive(Clone, Debug)]
pub struct AutomorphismIterator<F> {
    tree: IsoTree,
    g: F,
}

impl<F: Canonize> AutomorphismIterator<F> {
    /// Iterator on the automorphisms of `g` that fix `0..sigma`.
    fn new(g: &F, sigma: usize) -> Self {
        debug_assert!(g == &canonical_constraint(g, sigma).0);
        Self {
            tree: IsoTree::new(g, sigma),
            g: g.clone(),
        }
    }
}

impl<F: Canonize> Iterator for AutomorphismIterator<F> {
    type Item = Vec<usize>;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        while !self.tree.is_done() {
            let automorphism = match self.tree.leaf() {
                Some(phi) if self.g.apply_morphism(phi) == self.g => Some(phi.to_vec()),
                _ => None,
            };
            self.tree.advance();
            if automorphism.is_some() {
                return automorphism;
            }
        }
        None
    }
}

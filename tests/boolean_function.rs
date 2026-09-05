extern crate canonical_form;

use canonical_form::Canonize;

/// A boolean function on `n` variables, given by its truth table.
/// `table[x]` is `f` evaluated at the input where bit `i` of `x` is the
/// value of variable `i`. Canonization is up to permutation of variables
/// (the "P" in NPN-classification).
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct BoolFn {
    n: usize,
    table: Vec<bool>,
}

impl BoolFn {
    fn new(n: usize, table: Vec<bool>) -> Self {
        assert_eq!(table.len(), 1 << n);
        BoolFn { n, table }
    }
}

impl Canonize for BoolFn {
    fn size(&self) -> usize {
        self.n
    }

    fn apply_morphism(&self, p: &[usize]) -> Self {
        let n = self.n;
        let mut table = vec![false; 1 << n];
        for y in 0..(1usize << n) {
            let mut z = 0usize;
            for i in 0..n {
                if (y >> p[i]) & 1 == 1 {
                    z |= 1 << i;
                }
            }
            table[y] = self.table[z];
        }
        BoolFn { n, table }
    }

    // Per-variable influence (number of inputs on which flipping that
    // variable changes the output). Invariant under relabeling, and lets
    // the refiner split variables instead of falling back to brute-forcing
    // all n! permutations.
    fn invariant_coloring(&self) -> Option<Vec<u64>> {
        let mut influence = vec![0u64; self.n];
        for u in 0..self.n {
            let bit = 1usize << u;
            let mut count = 0u64;
            for x in 0..(1usize << self.n) {
                if x & bit == 0 && self.table[x] != self.table[x | bit] {
                    count += 1;
                }
            }
            influence[u] = count;
        }
        Some(influence)
    }
}

fn majority(n: usize) -> BoolFn {
    let table = (0..(1usize << n)).map(|x| x.count_ones() as usize > n / 2).collect();
    BoolFn::new(n, table)
}

fn xor(n: usize) -> BoolFn {
    let table = (0..(1usize << n)).map(|x: usize| x.count_ones() % 2 == 1).collect();
    BoolFn::new(n, table)
}

fn permutations(n: usize) -> Vec<Vec<usize>> {
    fn helper(current: &mut Vec<usize>, remaining: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if remaining.is_empty() {
            out.push(current.clone());
            return;
        }
        for i in 0..remaining.len() {
            let v = remaining.remove(i);
            current.push(v);
            helper(current, remaining, out);
            current.pop();
            remaining.insert(i, v);
        }
    }
    let mut out = Vec::new();
    helper(&mut Vec::new(), &mut (0..n).collect(), &mut out);
    out
}

fn brute_force_automorphisms(f: &BoolFn) -> usize {
    permutations(f.n)
        .iter()
        .filter(|p| f.apply_morphism(p) == *f)
        .count()
}

#[test]
fn canonical_invariant_under_relabeling() {
    // x0 AND (x1 OR x2)
    let f = BoolFn::new(3, vec![false, false, false, true, false, true, true, true]);
    for p in permutations(3) {
        assert_eq!(f.canonical(), f.apply_morphism(&p).canonical());
    }
}

#[test]
fn morphism_to_canonical() {
    let f = majority(3);
    let phi = f.morphism_to_canonical();
    assert_eq!(f.apply_morphism(&phi), f.canonical());
}

#[test]
fn fully_symmetric_functions_have_full_automorphism_group() {
    for n in 1..5 {
        let factorial: u64 = (1..=n as u64).product();
        assert_eq!(majority(n).canonical().automorphisms().count() as u64, factorial);
        assert_eq!(xor(n).canonical().automorphisms().count() as u64, factorial);
    }
}

#[test]
fn automorphisms_match_brute_force() {
    let cases = vec![
        // x0 AND (x1 OR x2): symmetric in x1, x2 only
        BoolFn::new(3, vec![false, false, false, true, false, true, true, true]),
        // x0 XOR x1, ignoring x2
        BoolFn::new(3, (0..8).map(|x: usize| (x & 1) ^ ((x >> 1) & 1) == 1).collect()),
        majority(4),
    ];
    for f in cases {
        let expected = brute_force_automorphisms(&f);
        let g = f.canonical();
        assert_eq!(g.automorphisms().count(), expected);
    }
}

#[test]
fn constant_function() {
    let f = BoolFn::new(2, vec![true, true, true, true]);
    let g = f.canonical();
    assert_eq!(g.automorphisms().count(), 2); // both variables irrelevant, S_2
}

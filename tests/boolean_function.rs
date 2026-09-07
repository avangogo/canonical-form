mod common;

use canonical_form::Canonize;
use common::check_contract;

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
        for (y, slot) in table.iter_mut().enumerate() {
            let mut z = 0usize;
            for (i, &pi) in p.iter().enumerate() {
                if (y >> pi) & 1 == 1 {
                    z |= 1 << i;
                }
            }
            *slot = self.table[z];
        }
        BoolFn { n, table }
    }

    // Per-variable influence (number of inputs on which flipping that
    // variable changes the output). Invariant under relabeling, and lets
    // the refiner split variables instead of falling back to brute-forcing
    // all n! permutations.
    fn invariant_color(&self, u: usize) -> u64 {
        let bit = 1usize << u;
        let mut count = 0u64;
        for x in 0..(1usize << self.n) {
            if x & bit == 0 && self.table[x] != self.table[x | bit] {
                count += 1;
            }
        }
        count
    }
}

fn majority(n: usize) -> BoolFn {
    let table = (0..(1usize << n))
        .map(|x| x.count_ones() as usize > n / 2)
        .collect();
    BoolFn::new(n, table)
}

fn xor(n: usize) -> BoolFn {
    let table = (0..(1usize << n))
        .map(|x: usize| x.count_ones() % 2 == 1)
        .collect();
    BoolFn::new(n, table)
}

#[test]
fn fully_symmetric_functions_have_full_automorphism_group() {
    for n in 1..5 {
        let factorial: u64 = (1..=n as u64).product();
        assert_eq!(
            majority(n).canonical().automorphisms().count() as u64,
            factorial
        );
        assert_eq!(xor(n).canonical().automorphisms().count() as u64, factorial);
    }
}

#[test]
fn constant_function() {
    let f = BoolFn::new(2, vec![true, true, true, true]);
    let g = f.canonical();
    assert_eq!(g.automorphisms().count(), 2); // both variables irrelevant, S_2
}

#[test]
fn the_contract_holds_on_every_function_of_three_variables() {
    for n in 1..=3 {
        for bits in 0..(1u32 << (1 << n)) {
            let table = (0..(1 << n)).map(|x| (bits >> x) & 1 == 1).collect();
            check_contract(&BoolFn::new(n, table));
        }
    }
}

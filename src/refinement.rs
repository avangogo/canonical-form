//! The ordered partition the refinement and the search work on.

use crate::hash::pointed_add;

/// An ordered partition of `0..n`, its parts named by an index in `0..k`.
#[derive(Clone, Debug)]
pub struct Partition {
    /// The elements, grouped by part.
    elems: Vec<usize>,
    /// `elems[rev_elems[i]] = i`
    rev_elems: Vec<usize>,
    /// `i` is in the part `set_id[i]`
    set_id: Vec<usize>,
    sets: Vec<Set>,
    /// What [`Self::split`] cuts by and resets, `0` for an untouched element.
    sieve: Vec<u64>,
    /// The parts holding an `e` with `sieve[e] != 0`.
    touched: Vec<usize>,
}

/// The part `elems[begin..end]`, of which `elems[begin..mid]` is touched.
#[derive(Clone, Debug, Copy)]
struct Set {
    begin: usize,
    end: usize,
    mid: usize,
}

impl Set {
    const fn len(&self) -> usize {
        self.end - self.begin
    }
}

impl Partition {
    /// The partition `{0}...{k-1}{k..size-1}`.
    pub fn with_singletons(size: usize, k: usize) -> Self {
        assert!(k <= size);
        let mut sets = Vec::with_capacity(size);
        let num_parts = if k == size { k } else { k + 1 };
        let mut set_id = if size > 0 {
            vec![num_parts - 1; size]
        } else {
            Vec::new()
        };
        for (i, set_i) in set_id.iter_mut().enumerate().take(k) {
            sets.push(Set {
                begin: i,
                mid: i,
                end: i + 1,
            });
            *set_i = i;
        }
        if size > k {
            sets.push(Set {
                begin: k,
                mid: k,
                end: size,
            });
        }
        Self {
            elems: (0..size).collect(),
            rev_elems: (0..size).collect(),
            set_id,
            sets,
            sieve: vec![0; size],
            touched: Vec::with_capacity(size),
        }
    }

    #[inline]
    fn swap(&mut self, i1: usize, i2: usize) {
        if i1 != i2 {
            let e1 = self.elems[i1];
            let e2 = self.elems[i2];
            self.elems[i1] = e2;
            self.elems[i2] = e1;
            self.rev_elems[e1] = i2;
            self.rev_elems[e2] = i1;
        }
    }

    /// Add `x` to the sieve value of `e`, and move `e` into the touched prefix
    /// of its part. `0` means untouched, hence [`pointed_add`] and not a sum.
    #[inline]
    pub fn sieve(&mut self, e: usize, x: u64) {
        if self.sieve[e] == 0 {
            let set_id = self.set_id[e];
            let set = &mut self.sets[set_id];
            if set.len() == 1 {
                return; // Nothing to split.
            }
            if set.mid == set.begin {
                self.touched.push(set_id);
            }
            // Move `e` into `elems[set.begin..set.mid]`.
            let new_pos = set.mid;
            set.mid += 1;
            self.swap(new_pos, self.rev_elems[e]);
        }
        self.sieve[e] = pointed_add(self.sieve[e], x);
    }

    /// Split every touched part by its sieve values, resetting them, and call
    /// `callback` on each part created.
    pub fn split<F>(&mut self, mut callback: F)
    where
        F: FnMut(usize),
    {
        // Parts are touched in a labelling-dependent order, and cut in one
        // that is not.
        self.touched.sort_unstable();
        // By index, so that `self` is free for `cut` to borrow.
        for i in 0..self.touched.len() {
            self.cut(self.touched[i], &mut callback);
        }
        self.touched.clear();
    }

    /// Cut the touched prefix of `s` into one part per sieve value, reporting
    /// them to `callback`. `s` keeps the last value.
    fn cut(&mut self, s: usize, callback: &mut impl FnMut(usize)) {
        let set = self.sets[s];
        let (begin, end) = (set.begin, set.mid);
        self.sets[s].mid = begin;
        self.sort_by_sieve(begin, end);

        // The key starts past the touched prefix, from an element untouched
        // whenever that prefix is short. Its `0` is a value no touched element
        // can carry, so the untouched tail keeps `s`: that is what `0` is
        // reserved for, and why sieve values go through `pointed_add`.
        let mut current_set = s;
        let mut current_key = self.sieve[self.elems[set.end - 1]];
        for i in (begin..end).rev() {
            let e = self.elems[i];
            if self.sieve[e] != current_key {
                current_key = self.sieve[e];
                self.sets[current_set].begin = i + 1;
                self.sets[current_set].mid = i + 1;
                self.sets.push(Set {
                    begin,
                    mid: begin,
                    end: i + 1,
                });
                current_set = self.num_parts() - 1;
                callback(current_set);
            }
            self.set_id[e] = current_set;
            self.sieve[e] = 0;
        }
    }

    /// Order `elems[begin..end]` by sieve value, keeping `rev_elems` in step.
    fn sort_by_sieve(&mut self, begin: usize, end: usize) {
        let sieve = &self.sieve;
        let first_key = sieve[self.elems[begin]];
        // The common case is a part that does not split at all.
        if self.elems[begin..end]
            .iter()
            .all(|&e| sieve[e] == first_key)
        {
            return;
        }
        self.elems[begin..end].sort_unstable_by_key(|&e| sieve[e]);
        for i in begin..end {
            self.rev_elems[self.elems[i]] = i;
        }
    }

    /// Separate the elements that `value` gives different values.
    ///
    /// Every element of a part that can still be split gets one, so the part is
    /// touched in full and the bookkeeping [`Self::sieve`] does per element has
    /// nothing to do. Singletons are never written to, nor cleared.
    pub fn refine_by_value(&mut self, value: impl Fn(usize) -> u64) {
        for part in 0..self.num_parts() {
            let set = self.sets[part];
            if set.len() >= 2 {
                for i in set.begin..set.end {
                    let e = self.elems[i];
                    // Not a plain write: an element of a touched prefix must
                    // not read as untouched, and `0` is a value `value` gives.
                    self.sieve[e] = pointed_add(0, value(e));
                }
                self.sets[part].mid = set.end;
                self.touched.push(part);
            }
        }
        self.split(|_| ());
    }

    /// The part that `s` was cut off from. Both [`Self::cut`] and
    /// [`Self::individualize`] take a prefix of the range of the part they
    /// leave behind, so the element past the end of `s` is still in its
    /// parent.
    fn parent_set(&self, s: usize) -> usize {
        self.set_id[self.elems[self.sets[s].end]]
    }

    /// Undo every part created after the first `n_parts`.
    pub fn undo(&mut self, n_parts: usize) {
        for s in (n_parts..self.num_parts()).rev() {
            let set = self.sets[s];
            let parent = self.parent_set(s);
            for e in &self.elems[set.begin..set.end] {
                self.set_id[*e] = parent;
            }
            self.sets[parent].begin = set.begin;
            self.sets[parent].mid = set.begin;
        }
        self.sets.truncate(n_parts);
    }

    #[inline]
    pub fn part(&self, part: usize) -> &[usize] {
        let set = self.sets[part];
        &self.elems[set.begin..set.end]
    }

    #[inline]
    pub const fn num_parts(&self) -> usize {
        self.sets.len()
    }

    #[inline]
    pub const fn num_elems(&self) -> usize {
        self.elems.len()
    }

    #[inline]
    /// Whether every part is a singleton.
    pub const fn is_discrete(&self) -> bool {
        self.num_elems() == self.num_parts()
    }

    /// Put `e`, which must not be alone in its part, in a part of its own.
    pub fn individualize(&mut self, e: usize) -> usize {
        let s = self.set_id[e];
        debug_assert!(
            self.sets[s].len() >= 2,
            "individualizing {e}, which is alone in its part already"
        );
        let i = self.rev_elems[e];
        self.swap(i, self.sets[s].begin);
        let delimiter = self.sets[s].begin + 1;
        let new_set = self.num_parts();
        self.set_id[e] = new_set;
        let new_begin = self.sets[s].begin;
        self.sets.push(Set {
            begin: new_begin,
            mid: new_begin,
            end: delimiter,
        });
        self.sets[s].begin = delimiter;
        self.sets[s].mid = delimiter;
        new_set
    }

    /// The bijection sending each element to the position of its part, when
    /// every part is a singleton.
    pub fn as_bijection(&self) -> Option<&[usize]> {
        if self.is_discrete() {
            Some(&self.rev_elems)
        } else {
            None
        }
    }

    /// A smallest part with at least two elements, ties broken by position.
    pub fn smallest_non_singleton(&self) -> Option<usize> {
        self.sets
            .iter()
            .enumerate()
            .filter(|(_, set)| set.len() >= 2)
            .min_by_key(|(_, set)| (set.len(), set.begin))
            .map(|(i, _)| i)
    }

    /// Panic unless every field agrees with every other. A full scan, for the
    /// refinement tests to run on a partition at rest.
    #[cfg(test)]
    pub(crate) fn check_consistent(&self) {
        let n = self.elems.len();
        assert_eq!(self.rev_elems.len(), n);
        for i in 0..n {
            assert_eq!(self.rev_elems[self.elems[i]], i);
            assert_eq!(self.elems[self.rev_elems[i]], i);
        }
        for (i, set) in self.sets.iter().enumerate() {
            assert!(set.begin < set.end);
            assert!(set.begin <= set.mid);
            assert!(set.mid <= set.end);
            for j in set.begin..set.end {
                assert_eq!(self.set_id[self.elems[j]], i);
            }
            for j in set.begin..set.mid {
                assert_ne!(self.sieve[self.elems[j]], 0);
            }
            for j in set.mid..set.end {
                assert_eq!(self.sieve[self.elems[j]], 0);
            }
        }
    }
}

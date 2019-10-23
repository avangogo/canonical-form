//! Data structure for refining partitions with extra properties for canonization of combinatorical
//! structures.

use std::fmt::{Display, Formatter, Result};

/// A partition of a set `{0..n-1}`.
/// Each part of the partition is refered by an index in `{0..k-1}`,
/// where `k` is the number of parts.
#[derive(Clone, Debug)]
pub struct Partition {
    /// Vector of contained elements
    elems: Vec<usize>,
    /// `elems[rev_elem[i]] = i`
    rev_elems: Vec<usize>,
    /// i is in part set_id[i]
    set_id: Vec<usize>,
    /// List of parts indexed by their id
    sets: Vec<Set>,
    /// mechanism to split sets:
    /// `sieve` contains a value for element (0 by default)
    /// when split is called, the partition is refined according to the values in `sieve`
    /// and `sieve` is reset.
    sieve: Vec<u64>,
    /// Contains the sets containing a `e` with `sieve[e] != 0`
    touched: Vec<usize>,
}

/// A part of a partition.
/// The part `s` is `elems[s.begin..s.end]`
/// The subset of `s` with non-zero `sieve` is `elems[s.begin..s.mid]`
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
    /// Return the partition of `{0..n-1}` with one part.
    pub fn simple(size: usize) -> Self {
        Self::with_singletons(size, 0)
    }

    /// Return the partition of `{0}...{k}{k+1..n-1}` with k+1 parts.
    pub fn with_singletons(size: usize, k: usize) -> Self {
        assert!(k <= size);
        let mut sets = Vec::with_capacity(size);
        let num_parts = if k == size { k } else { k + 1 };
        let mut set_id = if size > 0 {
            vec![num_parts - 1; size]
        } else {
            Vec::new()
        };
        // Create the singletons
        for (i, set_i) in set_id.iter_mut().enumerate().take(k) {
            sets.push(Set {
                begin: i,
                mid: i,
                end: i + 1,
            });
            *set_i = i
        }
        // Create the set with the rest, if it is non-empty
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

    /// Exchange the elements of indices `i1` and `i2` of a partition.
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

    /// Add `x` to the sieve value of `e` and update the related administration
    pub fn sieve(&mut self, e: usize, x: u64) {
        if self.sieve[e] == 0 {
            let set = &mut self.sets[self.set_id[e]];
            if set.len() == 1 {
                // A part of size one cannot be split further: we ignore the call
                return;
            };
            if set.mid == set.begin {
                self.touched.push(self.set_id[e]);
            };
            // update the partition so that `e` is in `elems[s.begin..s.mid]`
            let new_pos = set.mid;
            set.mid += 1;
            self.swap(new_pos, self.rev_elems[e]);
        }
        self.sieve[e] += x;
    }

    /// Split the partitions according to the values in sieve
    /// Call callback on each partition created
    pub fn split<F>(&mut self, mut callback: F)
    where
        F: FnMut(usize),
    {
        self.touched.sort();
        for &s in &self.touched {
            let set = self.sets[s];
            let begin = set.begin;
            let end = set.mid;
            self.sets[s].mid = begin;
            let sieve = &self.sieve;

            self.elems[begin..end].sort_by_key(|e| sieve[*e]);

            let mut current_set = s;
            let mut current_key = self.sieve[self.elems[set.end - 1]];

            for i in (begin..end).rev() {
                let elem_i = self.elems[i];
                if self.sieve[elem_i] != current_key {
                    current_key = self.sieve[elem_i];
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
                self.set_id[elem_i] = current_set;
                self.rev_elems[elem_i] = i;
                self.sieve[elem_i] = 0;
            }
        }
        self.touched.clear();
    }

    /// Separate elements with different keys.
    pub fn refine_by_value<F>(&mut self, key: &[u64], callback: F)
    where
        F: FnMut(usize),
    {
        for (i, &key) in key.iter().enumerate() {
            self.sieve(i, key);
        }
        self.split(callback)
    }

    fn parent_set(&self, s: usize) -> usize {
        self.set_id[self.elems[self.sets[s].end]]
    }

    /// Delete the last sets created until there are only nsets left
    pub fn undo(&mut self, nparts: usize) {
        for s in (nparts..self.num_parts()).rev() {
            let set = self.sets[s];
            let parent = self.parent_set(s);
            for e in &mut self.elems[set.begin..set.end] {
                self.set_id[*e] = parent
            }
            self.sets[parent].begin = set.begin;
            self.sets[parent].mid = set.begin;
        }
        self.sets.truncate(nparts);
    }

    #[inline]
    /// Return the slice of the elements of the part `part`.
    pub fn part(&self, part: usize) -> &[usize] {
        &self.elems[self.sets[part].begin..self.sets[part].end]
    }

    #[inline]
    /// Number of parts in the partition.
    pub fn num_parts(&self) -> usize {
        self.sets.len()
    }

    /// Number of elememts in the partition.
    pub fn num_elems(&self) -> usize {
        self.elems.len()
    }

    #[inline]
    /// Return `true` if the partition contains only cells of size 1.
    pub fn is_discrete(&self) -> bool {
        self.elems.len() == self.sets.len()
    }

    /// Refine the partition such that `e` is in a cell of size 1.
    pub fn individualize(&mut self, e: usize) -> Option<usize> {
        let s = self.set_id[e];
        if self.sets[s].end - self.sets[s].begin >= 2 {
            let i = self.rev_elems[e];
            self.swap(i, self.sets[s].begin);
            let delimiter = self.sets[s].begin + 1;
            //
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
            Some(new_set)
        } else {
            None
        }
    }

    /// If the partition conatains only cells of size 1, returns the bijection
    /// that map each element to the position of the corresponding cell.
    pub fn as_bijection(&self) -> Option<&[usize]> {
        if self.is_discrete() {
            Some(&self.rev_elems)
        } else {
            None
        }
    }

    /// Return the list of the cell in the order of the partition.
    pub const fn parts(&self) -> PartsIterator<'_> {
        PartsIterator {
            partition: self,
            pos: 0,
        }
    }

    /// Panic if the data structure does not correspond to a partition.
    fn check_consistent(&self) {
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
                assert_eq!(self.set_id[self.elems[j]], i)
            }
            for j in set.begin..set.mid {
                assert!(self.sieve[j] != 0)
            }
            for j in set.mid..set.end {
                assert_eq!(self.sieve[j], 0)
            }
        }
    }
}

pub struct PartsIterator<'a> {
    partition: &'a Partition,
    pos: usize,
}

impl<'a> Iterator for PartsIterator<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(&e) = self.partition.elems.get(self.pos) {
            let s = self.partition.set_id[e];
            self.pos = self.partition.sets[s].end;
            Some(s)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let n_parts = self.partition.num_parts();
        let lower = if self.pos < n_parts {
            n_parts - self.pos
        } else {
            0
        };
        (lower, Some(n_parts))
    }
}

impl Display for Partition {
    fn fmt(&self, f: &mut Formatter) -> Result {
        self.check_consistent();
        write!(f, "(")?;
        for i in 0..self.elems.len() {
            if i > 0 {
                if self.sets[self.set_id[self.elems[i]]].begin == i {
                    write!(f, ")(")?;
                } else {
                    write!(f, ",")?;
                }
            }
            write!(f, "{}", self.elems[i])?;
        }
        write!(f, ")")?;
        Ok(())
    }
}

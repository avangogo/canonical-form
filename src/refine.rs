//! Data structure for refining partitions.

use std::fmt::Display;
use std::fmt::*;

/// A refining partition of a set `{0..n-1}`.
#[derive(Clone, Debug)]
pub struct Partition {
    elems: Vec<usize>,     // vector of elements
    rev_elems: Vec<usize>, // elems[rev_elem[i]] = i
    set_id: Vec<usize>,    // i is in part set_id[i]
    begin: Vec<usize>,
    end: Vec<usize>,
}

impl Partition {
    /// Return the partition of `{0..n-1}` with one part.
    pub fn simple(size: usize) -> Self {
        Partition {
            // allow zero size partition ?
            elems: (0..size).collect(),
            rev_elems: (0..size).collect(),
            set_id: vec![0; size],
            begin: if size > 0 { vec![0] } else { vec![] },
            end: if size > 0 { vec![size] } else { vec![] },
        }
    }

    /// Exchange the elements of indices `i1` and `i2` of a partition.
    #[inline]
    fn swap(&mut self, i1: usize, i2: usize) {
        if i1 != i2 {
            debug_assert_eq!(self.set_id[self.elems[i1]], self.set_id[self.elems[i2]]);
            let e1 = self.elems[i1];
            let e2 = self.elems[i2];
            self.elems[i1] = e2;
            self.elems[i2] = e1;
            self.rev_elems[e1] = i2;
            self.rev_elems[e2] = i1;
        }
    }

    /// Separate elements with different keys.
    pub fn refine_by_value<F, T>(&mut self, key: &[T], mut callback: F)
    where
        F: FnMut(usize),
        T: Ord,
    {
        for s in 0..self.num_parts() {
            let begin = self.begin[s];
            let end = self.end[s];
            if end > begin + 1 {
                self.elems[begin..end].sort_by_key(|&e| &key[e]);
                for i in begin..end {
                    self.rev_elems[self.elems[i]] = i;
                }
                let mut current_set = s;
                let mut current_key = &key[self.elems[begin]];
                for i in (begin + 1)..end {
                    if key[self.elems[i]] != *current_key {
                        current_key = &key[self.elems[i]];
                        self.end[current_set] = i;
                        self.begin.push(i);
                        self.end.push(end);
                        current_set = self.num_parts() - 1;
                        callback(current_set);
                    }
                    self.set_id[self.elems[i]] = current_set;
                }
            }
        }
    }

    #[inline]
    /// Return the slice of the elements of the part `part`.
    pub fn part(&self, part: usize) -> &[usize] {
        &self.elems[self.begin[part]..self.end[part]]
    }

    #[inline]
    /// Number of parts in the partition.
    pub fn num_parts(&self) -> usize {
        self.begin.len()
    }

    #[inline]
    /// Return `true` if the partition contains only cells of size 1.
    pub fn is_discrete(&self) -> bool {
        self.elems.len() == self.begin.len()
    }

    /// Refine the partition such that `e` is in a cell of size 1.
    pub fn individualize(&mut self, e: usize) {
        let s = self.set_id[e];
        if self.end[s] - self.begin[s] >= 2 {
            let i = self.rev_elems[e];
            self.swap(i, self.begin[s]);
            let delimiter = self.begin[s] + 1;
            //
            self.set_id[e] = self.num_parts();
            self.begin.push(self.begin[s]);
            self.end.push(delimiter);
            self.begin[s] = delimiter;
            self.check_consistent()
        }
    }

    /// If the partition conatains only cells of size 1, returns the bijection
    /// that map each element to the position of the corresponding cell.
    pub fn to_bijection(&self) -> Option<Vec<usize>> {
        if !self.is_discrete() {
            None
        } else {
            Some(self.rev_elems.clone())
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
        let p = self.begin.len();
        assert_eq!(self.end.len(), p);
        for i in 0..p {
            assert!(self.begin[i] < self.end[i]);
            for j in self.begin[i]..self.end[i] {
                assert_eq!(self.set_id[self.elems[j]], i)
            }
        }
    }

    /// Return the list of the cell in the order of the partition.
    pub fn parts(&self) -> Vec<usize> {
        let mut res = Vec::new();
        let mut i = 0;
        while i < self.elems.len() {
            let s = self.set_id[self.elems[i]];
            res.push(s);
            i = self.end[s]
        }
        let mut test: Vec<_> = self.elems.iter().map(|&i| self.set_id[i]).collect();
        test.dedup();
        assert_eq!(res, test);
        res
    }
}

impl Display for Partition {
    fn fmt(&self, f: &mut Formatter) -> Result {
        self.check_consistent();
        write!(f, "(").unwrap();
        for i in 0..self.elems.len() {
            if i > 0 {
                if self.begin[self.set_id[self.elems[i]]] == i {
                    write!(f, ")(").unwrap();
                } else {
                    write!(f, ",").unwrap();
                }
            }
            write!(f, "{}", self.elems[i]).unwrap();
        }
        write!(f, ")").unwrap();
        Ok(())
    }
}

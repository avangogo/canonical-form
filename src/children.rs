//! One buffer for the children of every node on a branch.

/// The children of the nodes of a branch, as a stack of segments.
///
/// Nodes are created and dropped in stack order, so one buffer serves the whole
/// branch. **Only the top segment reaches the end of the buffer**, which is why
/// every method but [`Self::push`] is meant for the node being visited.
#[derive(Clone, Debug, Default)]
pub struct ChildStack {
    buffer: Vec<usize>,
}

/// One node's children inside a [`ChildStack`]: `begin..next` still to
/// individualize, `next..` already taken.
#[derive(Clone, Copy, Debug)]
pub struct Segment {
    begin: usize,
    next: usize,
}

impl Segment {
    pub const fn is_exhausted(&self) -> bool {
        self.next == self.begin
    }
}

impl ChildStack {
    pub fn push(&mut self, children: &[usize]) -> Segment {
        let begin = self.buffer.len();
        self.buffer.extend_from_slice(children);
        Segment {
            begin,
            next: self.buffer.len(),
        }
    }

    /// Take the next child of `segment`, `None` once it has none left.
    pub fn take(&self, segment: &mut Segment) -> Option<usize> {
        if segment.is_exhausted() {
            return None;
        }
        segment.next -= 1;
        Some(self.buffer[segment.next])
    }

    pub fn taken(&self, segment: &Segment) -> &[usize] {
        &self.buffer[segment.next..]
    }

    pub fn drop_top(&mut self, segment: &Segment) {
        self.buffer.truncate(segment.begin);
    }

    /// Keep the children of the top segment that `keep` accepts. The ones
    /// already taken move down to close the gap, so the segment still reaches
    /// the end of the buffer.
    pub fn retain_untaken(&mut self, segment: &mut Segment, mut keep: impl FnMut(usize) -> bool) {
        let mut kept = segment.begin;
        for i in segment.begin..segment.next {
            let child = self.buffer[i];
            if keep(child) {
                self.buffer[kept] = child;
                kept += 1;
            }
        }
        let taken = self.buffer.len() - segment.next;
        self.buffer.copy_within(segment.next.., kept);
        segment.next = kept;
        self.buffer.truncate(kept + taken);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn children_come_back_in_reverse() {
        let mut stack = ChildStack::default();
        let mut segment = stack.push(&[4, 7, 9]);
        assert_eq!(stack.take(&mut segment), Some(9));
        assert_eq!(stack.take(&mut segment), Some(7));
        assert_eq!(stack.take(&mut segment), Some(4));
        assert_eq!(stack.take(&mut segment), None);
        assert!(segment.is_exhausted());
    }

    #[test]
    fn taken_children_stay_available() {
        let mut stack = ChildStack::default();
        let mut segment = stack.push(&[1, 2, 3]);
        assert!(stack.taken(&segment).is_empty());
        let _ = stack.take(&mut segment);
        let _ = stack.take(&mut segment);
        assert_eq!(stack.taken(&segment), [2, 3]);
    }

    #[test]
    fn a_segment_only_sees_its_own_children() {
        let mut stack = ChildStack::default();
        let mut lower = stack.push(&[1, 2]);
        let mut upper = stack.push(&[8, 9]);
        assert_eq!(stack.take(&mut upper), Some(9));
        stack.drop_top(&upper);
        // The lower segment is untouched and is the top again.
        assert_eq!(stack.take(&mut lower), Some(2));
        assert_eq!(stack.taken(&lower), [2]);
    }

    #[test]
    fn retain_keeps_the_taken_children_at_the_top() {
        let mut stack = ChildStack::default();
        let _lower = stack.push(&[100]);
        let mut segment = stack.push(&[1, 2, 3, 4, 5]);
        let _ = stack.take(&mut segment); // 5
        let _ = stack.take(&mut segment); // 4
        stack.retain_untaken(&mut segment, |v| v % 2 == 1);
        assert_eq!(stack.taken(&segment), [4, 5]);
        assert_eq!(stack.take(&mut segment), Some(3));
        assert_eq!(stack.take(&mut segment), Some(1));
        assert_eq!(stack.take(&mut segment), None);
        assert_eq!(stack.taken(&segment), [1, 3, 4, 5]);
    }

    #[test]
    fn retain_everything_and_nothing() {
        let mut stack = ChildStack::default();
        let mut segment = stack.push(&[1, 2, 3]);
        stack.retain_untaken(&mut segment, |_| true);
        assert_eq!(stack.take(&mut segment), Some(3));
        stack.retain_untaken(&mut segment, |_| false);
        assert!(segment.is_exhausted());
        assert_eq!(stack.taken(&segment), [3]);
    }
}

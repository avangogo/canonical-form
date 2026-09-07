//! Fingerprinting a multiset of colors: hash each with [`mix`], sum them with
//! [`pointed_add`]. The sum wraps, so two multisets collide with probability
//! about `2^-63`, which makes them indistinguishable and nothing worse. The
//! constants are fixed, so fingerprints are reproducible across processes.

/// The `splitmix64` finalizer, preceded by an increment.
///
/// Every step is invertible, so this is a permutation of `u64` and no collision
/// comes from here. It makes colors opaque: raw, two colors `21` would sum like
/// one color `42`. The increment is there because the bare finalizer fixes `0`,
/// which would be invisible in a sum.
pub const fn mix(z: u64) -> u64 {
    let mut z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Bit set on every result of [`pointed_add`], keeping it clear of `0`.
const HIGH_BIT: u64 = 1 << 63;

/// Add `a` to `b`, never returning `0`: the low 63 bits add modulo `2^63`, and
/// the high bit is set, which leaves `0` to mean "nothing added".
#[inline]
pub const fn pointed_add(a: u64, b: u64) -> u64 {
    HIGH_BIT | a.wrapping_add(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pointed_add_is_never_zero() {
        assert_ne!(pointed_add(0, 0), 0);
        // A pair of increments that cancel each other must not read as `0`.
        let x = mix(7);
        assert_ne!(pointed_add(pointed_add(0, x), x.wrapping_neg()), 0);
    }

    #[test]
    fn pointed_add_does_not_depend_on_order() {
        // The increments of a multiset arrive in an arbitrary order, so the
        // fold must depend only on the multiset itself.
        let mut xs: Vec<u64> = (0..12u64).map(mix).collect();
        xs.extend([1 << 63, u64::MAX, 0, 1]);
        let fold = |order: &[u64]| order.iter().fold(0u64, |acc, &x| pointed_add(acc, x));
        let expected = fold(&xs);
        for k in 0..xs.len() {
            let mut p = xs.clone();
            p.rotate_left(k);
            assert_eq!(fold(&p), expected);
            p.reverse();
            assert_eq!(fold(&p), expected);
        }
    }

    #[test]
    fn pointed_add_adds_modulo_two_to_the_63() {
        let mut acc = 0;
        let mut expected = 0u64;
        for x in [1, u64::MAX, 1 << 62, mix(3), 1 << 62, mix(99)] {
            acc = pointed_add(acc, x);
            expected = expected.wrapping_add(x) & !HIGH_BIT;
            assert_eq!(acc & !HIGH_BIT, expected);
            assert_ne!(acc, 0);
        }
    }

    #[test]
    fn no_small_color_hashes_to_zero() {
        // A color hashing to `0` would be neutral for the sum, so its
        // multiplicity would be invisible.
        assert!((0..100_000u64).all(|c| mix(c) != 0));
    }
}

//! A small seedable pseudo-random generator, standing in for the `rand` crate.
//!
//! Behemoth's conformer search uses `rand`. Beavyr keeps five dependencies on purpose, and what
//! the search actually needs is narrow -- uniform draws from a range, and a weighted coin -- so it
//! is provided here instead of adding a sixth.
//!
//! There is a second reason to prefer this. `from_entropy` in the original means two runs of the
//! same search give different answers, which is awkward when a user reports that a conformer
//! search missed something. [`StdRng::from_entropy`] here still varies between runs, but every
//! search can be made exactly reproducible by passing a seed, and the algorithm is fixed in this
//! file rather than tracking a dependency's version.
//!
//! The generator is xoshiro256**, seeded through SplitMix64. Both are public-domain designs by
//! Blackman and Vigna. This is not cryptographic and is not meant to be.
//!
//! The API deliberately mirrors the slice of `rand` it replaces -- `seed_from_u64`,
//! `from_entropy`, `gen_range`, `gen_bool` -- so the imported call sites need no edit.

use std::ops::{Range, RangeInclusive};

/// xoshiro256** state.
#[derive(Debug, Clone)]
pub struct StdRng {
    s: [u64; 4],
}

impl StdRng {
    /// A generator fixed by `seed`. Two runs with the same seed draw the same numbers.
    pub fn seed_from_u64(seed: u64) -> Self {
        // SplitMix64 expands one word into the four xoshiro needs, which is the standard way to
        // avoid a poor initial state from a small seed.
        let mut z = seed;
        let mut next = || {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            x ^ (x >> 31)
        };
        Self {
            s: [next(), next(), next(), next()],
        }
    }

    /// A generator seeded from the clock, for when reproducibility is not wanted.
    pub fn from_entropy() -> Self {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos() as u64)
            .unwrap_or(0x2545_F491_4F6C_DD1D);
        // Mixing in a stack address separates generators created in the same nanosecond on
        // different threads, which the conformer search does.
        let here = &nanos as *const u64 as u64;
        Self::seed_from_u64(nanos ^ here.rotate_left(17))
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// A uniform double in `[0, 1)`.
    #[inline]
    fn next_f64(&mut self) -> f64 {
        // The top 53 bits are the ones with full quality in xoshiro**, and 53 is exactly the
        // mantissa width, so this wastes nothing and never returns 1.0.
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// A uniform draw below `bound`, without the modulo bias a plain `%` would introduce.
    #[inline]
    fn below(&mut self, bound: u64) -> u64 {
        debug_assert!(bound > 0);
        // Lemire's method: reject only the short final window.
        let zone = u64::MAX - u64::MAX % bound;
        loop {
            let value = self.next_u64();
            if value < zone {
                return value % bound;
            }
        }
    }

    /// A uniform draw from `range`.
    pub fn gen_range<T, R: SampleRange<T>>(&mut self, range: R) -> T {
        range.sample(self)
    }

    /// `true` with probability `p`.
    pub fn gen_bool(&mut self, p: f64) -> bool {
        if p <= 0.0 {
            false
        } else if p >= 1.0 {
            true
        } else {
            self.next_f64() < p
        }
    }
}

/// Ranges that [`StdRng::gen_range`] can draw from.
///
/// Only the four shapes the imported search uses are implemented; adding another is a few lines,
/// and leaving them out means an unsupported draw is a compile error rather than a surprise.
pub trait SampleRange<T> {
    fn sample(self, rng: &mut StdRng) -> T;
}

impl SampleRange<f64> for Range<f64> {
    fn sample(self, rng: &mut StdRng) -> f64 {
        assert!(self.start < self.end, "empty f64 range");
        self.start + rng.next_f64() * (self.end - self.start)
    }
}

impl SampleRange<usize> for Range<usize> {
    fn sample(self, rng: &mut StdRng) -> usize {
        assert!(self.start < self.end, "empty usize range");
        self.start + rng.below((self.end - self.start) as u64) as usize
    }
}

impl SampleRange<usize> for RangeInclusive<usize> {
    fn sample(self, rng: &mut StdRng) -> usize {
        let (start, end) = (*self.start(), *self.end());
        assert!(start <= end, "empty usize range");
        start + rng.below((end - start) as u64 + 1) as usize
    }
}

impl SampleRange<i32> for RangeInclusive<i32> {
    fn sample(self, rng: &mut StdRng) -> i32 {
        let (start, end) = (*self.start(), *self.end());
        assert!(start <= end, "empty i32 range");
        start + rng.below((end as i64 - start as i64) as u64 + 1) as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_same_seed_gives_the_same_sequence() {
        let mut a = StdRng::seed_from_u64(12345);
        let mut b = StdRng::seed_from_u64(12345);
        for _ in 0..64 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn different_seeds_diverge() {
        let mut a = StdRng::seed_from_u64(1);
        let mut b = StdRng::seed_from_u64(2);
        assert!((0..8).any(|_| a.next_u64() != b.next_u64()));
    }

    #[test]
    fn ranges_stay_inside_their_bounds() {
        let mut rng = StdRng::seed_from_u64(7);
        for _ in 0..10_000 {
            let f: f64 = rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI);
            assert!(f >= -std::f64::consts::PI && f < std::f64::consts::PI);

            let u: usize = rng.gen_range(3..9);
            assert!((3..9).contains(&u));

            let ui: usize = rng.gen_range(1..=4);
            assert!((1..=4).contains(&ui));

            let i: i32 = rng.gen_range(-179..=180);
            assert!((-179..=180).contains(&i));
        }
    }

    /// An inclusive range of one value must be drawable, since the search does that when only one
    /// torsion is available to mutate.
    #[test]
    fn a_single_value_range_is_allowed() {
        let mut rng = StdRng::seed_from_u64(9);
        assert_eq!(rng.gen_range(5usize..=5), 5usize);
        assert_eq!(rng.gen_range(2usize..3), 2usize);
    }

    /// Every value in a small range must actually appear, or a biased draw would quietly shrink
    /// the search space.
    #[test]
    fn small_ranges_cover_every_value() {
        let mut rng = StdRng::seed_from_u64(11);
        let mut seen = [false; 6];
        for _ in 0..2_000 {
            seen[rng.gen_range(0..6usize)] = true;
        }
        assert!(seen.iter().all(|s| *s), "some values were never drawn");
    }

    #[test]
    fn gen_bool_respects_its_probability() {
        let mut rng = StdRng::seed_from_u64(13);
        assert!(!rng.gen_bool(0.0));
        assert!(rng.gen_bool(1.0));

        let trials = 20_000;
        let hits = (0..trials).filter(|_| rng.gen_bool(0.25)).count();
        let fraction = hits as f64 / trials as f64;
        assert!((fraction - 0.25).abs() < 0.02, "got {fraction}, want ~0.25");
    }

    /// A uniform double must fill [0, 1) without ever reaching 1.0, which would put `gen_range`
    /// one past the end of its range.
    #[test]
    fn uniform_doubles_stay_below_one() {
        let mut rng = StdRng::seed_from_u64(17);
        let mut lowest = f64::MAX;
        let mut highest = f64::MIN;
        for _ in 0..100_000 {
            let v = rng.next_f64();
            assert!((0.0..1.0).contains(&v));
            lowest = lowest.min(v);
            highest = highest.max(v);
        }
        assert!(lowest < 0.01 && highest > 0.99, "poor spread: {lowest}..{highest}");
    }
}

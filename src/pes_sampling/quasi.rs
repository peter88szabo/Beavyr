//! Random-number streams for PES sampling: pseudo-random (xoshiro) or a
//! scrambled Halton low-discrepancy sequence.
//!
//! Every accepted-or-rejected sample attempt is one point of the sequence.
//! Its random choices (mode energies/phases, TS displacement, rotor angles,
//! orientations) consume the point's coordinates in a fixed order, so each
//! choice always lands on the same sequence dimension. Dimension k uses the
//! k-th prime as base with a seeded random digit permutation (π(0) = 0), which
//! removes the correlations plain Halton shows between large bases.
use super::state::RandomNumbers;
use crate::optimizer::prng::StdRng;

/// Coordinates beyond this many per sample fall back to pseudo-random numbers.
pub const MAX_DIMENSIONS: usize = 512;

pub trait Unit {
    /// A number in [0, 1).
    fn unit(&mut self) -> f64;
}

impl Unit for StdRng {
    fn unit(&mut self) -> f64 {
        self.gen_range(0.0..1.0)
    }
}

pub struct Draws {
    rng: StdRng,
    halton: Option<Halton>,
}

struct Halton {
    index: u64,
    dimension: usize,
    primes: Vec<u64>,
    permutations: Vec<Vec<u32>>,
    scramble: StdRng,
}

pub fn primes(count: usize) -> Vec<u64> {
    let mut found: Vec<u64> = Vec::with_capacity(count);
    let mut n = 2u64;
    while found.len() < count {
        if found.iter().take_while(|p| *p * *p <= n).all(|p| n % p != 0) {
            found.push(n);
        }
        n += 1;
    }
    found
}

/// Radical inverse of `n` in base `p`, with digit permutation `perm`.
pub fn scrambled_radical_inverse(mut n: u64, p: u64, perm: &[u32]) -> f64 {
    let mut value = 0.;
    let mut factor = 1. / p as f64;
    while n > 0 {
        value += perm[(n % p) as usize] as f64 * factor;
        n /= p;
        factor /= p as f64;
    }
    value
}

impl Draws {
    pub fn new(kind: RandomNumbers, seed: u64) -> Self {
        let halton = (kind == RandomNumbers::QuasiRandom).then(|| Halton {
            index: 0,
            dimension: 0,
            primes: primes(MAX_DIMENSIONS),
            permutations: vec![],
            // Independent of the pseudo-random stream, still fixed by the seed.
            scramble: StdRng::seed_from_u64(seed ^ 0x9E37_79B9_7F4A_7C15),
        });
        Self {
            rng: StdRng::seed_from_u64(seed),
            halton,
        }
    }

    /// Start the next sample attempt: the next point of the sequence.
    pub fn next_point(&mut self) {
        if let Some(h) = &mut self.halton {
            h.index += 1;
            h.dimension = 0;
        }
    }
}

impl Unit for Draws {
    fn unit(&mut self) -> f64 {
        let Some(h) = &mut self.halton else {
            return self.rng.unit();
        };
        let k = h.dimension;
        h.dimension += 1;
        if k >= MAX_DIMENSIONS {
            return self.rng.unit();
        }
        while h.permutations.len() <= k {
            let p = h.primes[h.permutations.len()] as usize;
            // Fisher-Yates on 1..p, keeping 0 fixed so finite expansions stay finite.
            let mut perm: Vec<u32> = (0..p as u32).collect();
            for i in (2..p).rev() {
                let j = 1 + (h.scramble.unit() * i as f64) as usize;
                perm.swap(i, j.min(i));
            }
            h.permutations.push(perm);
        }
        scrambled_radical_inverse(h.index, h.primes[k], &h.permutations[k]).min(1. - f64::EPSILON)
    }
}

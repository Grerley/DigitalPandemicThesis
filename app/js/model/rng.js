// rng.js — Deterministic, seedable pseudo-random number generators.
// -----------------------------------------------------------------------------
// All stochastic components of the lab (network Monte-Carlo, synthetic cohort)
// draw from an explicit, seeded RNG so that "same seed -> same result" holds
// exactly, as required for reproducibility. We use `mulberry32`: a small, fast,
// well-distributed 32-bit generator. It is NOT cryptographic — it is chosen for
// speed and perfect reproducibility across engines.

/**
 * Create a seeded uniform generator on [0,1).
 * @param {number} seed integer seed
 * @returns {() => number}
 */
export function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0;
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** A small RNG object bundling common draws around a seeded uniform stream. */
export class RNG {
  constructor(seed = 1) {
    this.seed = seed >>> 0;
    this._u = mulberry32(this.seed);
    this._spare = null;
  }
  /** Uniform in [0,1). */
  next() {
    return this._u();
  }
  /** Uniform integer in [0, n). */
  int(n) {
    return Math.floor(this._u() * n);
  }
  /** Standard normal via Box–Muller (cached spare for efficiency). */
  normal(mean = 0, sd = 1) {
    if (this._spare !== null) {
      const s = this._spare;
      this._spare = null;
      return mean + sd * s;
    }
    let u = 0, v = 0, s = 0;
    do {
      u = this._u() * 2 - 1;
      v = this._u() * 2 - 1;
      s = u * u + v * v;
    } while (s >= 1 || s === 0);
    const mul = Math.sqrt((-2 * Math.log(s)) / s);
    this._spare = v * mul;
    return mean + sd * (u * mul);
  }
  /** Bernoulli(p). */
  bernoulli(p) {
    return this._u() < p ? 1 : 0;
  }
  /** Weighted categorical sample: returns index i with prob weights[i]/sum. */
  categorical(weights) {
    let sum = 0;
    for (const w of weights) sum += w;
    let r = this._u() * sum;
    for (let i = 0; i < weights.length; i++) {
      r -= weights[i];
      if (r < 0) return i;
    }
    return weights.length - 1;
  }
  /** In-place Fisher–Yates shuffle. */
  shuffle(arr) {
    for (let i = arr.length - 1; i > 0; i--) {
      const j = this.int(i + 1);
      [arr[i], arr[j]] = [arr[j], arr[i]];
    }
    return arr;
  }
  /** Sample k distinct integers from [0, n) without replacement. */
  sampleIndices(n, k) {
    const idx = Array.from({ length: n }, (_, i) => i);
    this.shuffle(idx);
    return idx.slice(0, Math.min(k, n));
  }
}

// Numerics for the cancellation demo: three ways to compute the population variance in
// plain doubles, and an accurate reference computed in double-double arithmetic.

/** Population variance, textbook one-pass formula: sum(x^2)/n - mean(x)^2. */
export function onePass(x: Float64Array): number {
  let s = 0, s2 = 0;
  for (let i = 0; i < x.length; i++) { s += x[i]; s2 += x[i] * x[i]; }
  const n = x.length, m = s / n;
  return s2 / n - m * m;
}

/** Population variance, two-pass: centre on the mean first, then sum the squares. */
export function twoPass(x: Float64Array): number {
  let s = 0;
  for (let i = 0; i < x.length; i++) s += x[i];
  const m = s / x.length;
  let q = 0;
  for (let i = 0; i < x.length; i++) { const d = x[i] - m; q += d * d; }
  return q / x.length;
}

/** Population variance, Welford's one-pass running update (also centres). */
export function welford(x: Float64Array): number {
  let m = 0, m2 = 0;
  for (let i = 0; i < x.length; i++) {
    const d = x[i] - m;
    m += d / (i + 1);
    m2 += d * (x[i] - m);
  }
  return m2 / x.length;
}

// ---- double-double arithmetic (about 32 significant digits) ----
type DD = [number, number];
function twoSum(a: number, b: number): DD {
  const s = a + b, bb = s - a;
  return [s, (a - (s - bb)) + (b - bb)];
}
function quickTwoSum(a: number, b: number): DD { const s = a + b; return [s, b - (s - a)]; }
function split(a: number): DD { const c = 134217729 * a, hi = c - (c - a); return [hi, a - hi]; }
function twoProd(a: number, b: number): DD {
  const p = a * b, [ah, al] = split(a), [bh, bl] = split(b);
  return [p, ((ah * bh - p) + ah * bl + al * bh) + al * bl];
}
function ddAdd(a: DD, b: DD): DD {
  const [s, e] = twoSum(a[0], b[0]);
  return quickTwoSum(s, e + a[1] + b[1]);
}
function ddDiv(a: DD, n: number): DD {
  const q1 = a[0] / n, p = twoProd(q1, n);
  const r = ddAdd(a, [-p[0], -p[1]]);
  return quickTwoSum(q1, (r[0] + r[1]) / n);
}

/**
 * Exact-to-double population variance of the stored numbers x (the reference "true"
 * value): mean and squared deviations in double-double, so no digits are lost even for
 * an offset of 1e12 with sd 0.01.
 */
export function refVariance(x: Float64Array): number {
  const n = x.length;
  let s: DD = [0, 0];
  for (let i = 0; i < n; i++) s = ddAdd(s, [x[i], 0]);
  const m = ddDiv(s, n);
  let q: DD = [0, 0];
  for (let i = 0; i < n; i++) {
    const [dh, e] = twoSum(x[i], -m[0]);
    const d = quickTwoSum(dh, e - m[1]);
    const p = twoProd(d[0], d[0]);
    q = ddAdd(q, quickTwoSum(p[0], p[1] + 2 * d[0] * d[1]));
  }
  return ddDiv(q, n)[0];
}

/** Correct significant digits, -log10(relative error), capped to [0, 16]. */
export function digits(est: number, ref: number): number {
  const rel = Math.abs(est - ref) / Math.abs(ref);
  if (!Number.isFinite(rel)) return 0;
  return Math.max(0, Math.min(16, -Math.log10(Math.max(rel, 1e-16))));
}

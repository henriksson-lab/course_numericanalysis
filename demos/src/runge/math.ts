// Numerics for the Runge explorer: the logistic OD curve, node sets, barycentric
// interpolation, the Lebesgue function and a least-squares polynomial fit.
import { rng } from "../lib/ui";

export const K = 1.2, N0 = 0.05, R = 0.8;
export const od = (t: number): number => K / (1 + (K / N0 - 1) * Math.exp(-R * t));

export type NodeSet = "equi" | "cheb";

/** n+1 sampling times on [0, T]: equally spaced, or Chebyshev points of the first kind. */
export function nodes(set: NodeSet, n: number, T: number): number[] {
  return Array.from({ length: n + 1 }, (_, k) => set === "equi"
    ? T * k / n
    : T / 2 - T / 2 * Math.cos((2 * k + 1) * Math.PI / (2 * n + 2)));
}

// The lecture's noisy hourly readings (logistic + N(0, 0.02^2), set.seed(1) in R).
export const Y12 = [0.0375, 0.1095, 0.1959, 0.4207, 0.6259, 0.8279, 1.0188, 1.1207, 1.1673, 1.1736, 1.2210, 1.2037, 1.1857];
export const SD = 0.02;

/** Reading errors for a design, from a fixed seed; the lecture's hourly 0-12 h design gets the lecture's values. */
export function noiseFor(set: NodeSet, n: number, T: number): number[] {
  if (set === "equi" && n === 12 && T === 12) return Y12.map((y, i) => y - od(i));
  const g = rng(1000 * T + 100 * n + (set === "cheb" ? 7 : 0) + 1);
  return Array.from({ length: n + 1 }, () => SD * g.norm());
}

/** Barycentric weights, computed on s = 2t/T - 1 so that they neither overflow nor underflow. */
export function baryWeights(t: number[], T: number): number[] {
  const s = t.map(v => 2 * v / T - 1);
  return s.map((sj, j) => {
    let p = 1;
    for (let k = 0; k < s.length; k++) if (k !== j) p *= sj - s[k];
    return 1 / p;
  });
}

/** Evaluate the interpolant at x (barycentric formula of the second kind). */
export function baryEval(t: number[], y: number[], w: number[], x: number): number {
  let num = 0, den = 0;
  for (let j = 0; j < t.length; j++) {
    const d = x - t[j];
    if (d === 0) return y[j];
    const c = w[j] / d;
    num += c * y[j]; den += c;
  }
  return num / den;
}

/** Lebesgue function sum_i |L_i(x)|. */
export function lebesgue(t: number[], w: number[], x: number): number {
  let num = 0, den = 0;
  for (let j = 0; j < t.length; j++) {
    const d = x - t[j];
    if (d === 0) return 1;
    const c = w[j] / d;
    num += Math.abs(c); den += c;
  }
  return num / Math.abs(den);
}

/** Least-squares polynomial of degree d (Chebyshev basis on [0, T], modified Gram-Schmidt QR). */
export function lsFit(t: number[], y: number[], d: number, T: number): { f: (x: number) => number; sd: number } {
  const m = t.length, p = Math.min(d, m - 1) + 1;
  const basis = (x: number): number[] => {
    const s = 2 * x / T - 1, b = [1, s];
    for (let k = 2; k < p; k++) b.push(2 * s * b[k - 1] - b[k - 2]);
    return b.slice(0, p);
  };
  const A = t.map(basis);                                   // m x p
  const Q = A.map(r => r.slice());
  const Rm: number[][] = Array.from({ length: p }, () => new Array(p).fill(0));
  for (let j = 0; j < p; j++) {
    for (let pass = 0; pass < 2; pass++)                    // re-orthogonalise once
      for (let i = 0; i < j; i++) {
        let r = 0;
        for (let k = 0; k < m; k++) r += Q[k][i] * Q[k][j];
        Rm[i][j] += r;
        for (let k = 0; k < m; k++) Q[k][j] -= r * Q[k][i];
      }
    let nr = 0;
    for (let k = 0; k < m; k++) nr += Q[k][j] ** 2;
    nr = Math.sqrt(nr);
    Rm[j][j] = nr;
    for (let k = 0; k < m; k++) Q[k][j] /= nr;
  }
  const qy = Array.from({ length: p }, (_, j) => { let s = 0; for (let k = 0; k < m; k++) s += Q[k][j] * y[k]; return s; });
  const c = new Array(p).fill(0);
  for (let j = p - 1; j >= 0; j--) {
    let s = qy[j];
    for (let i = j + 1; i < p; i++) s -= Rm[j][i] * c[i];
    c[j] = s / Rm[j][j];
  }
  const f = (x: number) => basis(x).reduce((a, b, i) => a + b * c[i], 0);
  let rss = 0;
  for (let k = 0; k < m; k++) rss += (y[k] - f(t[k])) ** 2;
  return { f, sd: m > p ? Math.sqrt(rss / (m - p)) : 0 };
}

/** Evaluation grid on [0, T]: step T/2400 (0.02 h for a 48-h run). */
export const GRID = 2400;
export const grid = (T: number): number[] => Array.from({ length: GRID + 1 }, (_, i) => T * i / GRID);

/** Largest |f - p| on the grid, with its location. */
export function maxErr(g: number[], p: (x: number) => number): { e: number; at: number } {
  let e = 0, at = 0;
  for (const x of g) { const v = Math.abs(p(x) - od(x)); if (v > e || !Number.isFinite(v)) { e = v; at = x; } }
  return { e, at };
}

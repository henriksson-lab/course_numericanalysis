// Numerics for the design demo: the oral one-compartment PK model of the
// identifiability unit, a Levenberg-Marquardt fit in log-parameters, the SVD of the
// log-parameter Jacobian (one-sided Jacobi), the covariance s^2 (J^T J)^{-1}, the
// profile of ka (ke and V re-fitted for each fixed ka) and the SSE over (ka, ke) with
// V re-fitted (closed form, since C is linear in 1/V).

export const DOSE = 100;
export const SIGMA = 0.1;                        // noise SD of the simulated data (mg/L)
export const TRUTH = [1, 0.2, 20];               // ka (1/h), ke (1/h), V (L)
export const T_DATA = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12];
export const Y_DATA = [1.02, 1.88, 2.34, 2.98, 3.27, 3.26, 3.17, 2.77, 1.92, 1.23, 1.00, 0.61];

/** ka/(ka-ke) (e^{-ke t} - e^{-ka t}), stable when ka is close to ke. */
export function shape(t: number, ka: number, ke: number): number {
  const d = ka - ke;
  if (Math.abs(d * t) < 1e-12) return ka * t * Math.exp(-ke * t);
  return ka * Math.exp(-ke * t) * (-Math.expm1(-d * t)) / d;
}
export const pk = (t: number, ka: number, ke: number, V: number): number => DOSE / V * shape(t, ka, ke);

/** The noise of the course data, in units of SIGMA: each sample keeps its own value. */
export const EPS_DATA = T_DATA.map((t, i) => (Y_DATA[i] - pk(t, TRUTH[0], TRUTH[1], TRUTH[2])) / SIGMA);

const sse = (r: number[]): number => r.reduce((s, v) => s + v * v, 0);

/** Solve a small dense system A x = b (Gaussian elimination, partial pivoting). */
function solve(A: number[][], b: number[]): number[] | null {
  const n = b.length, M = A.map((row, i) => [...row, b[i]]);
  for (let k = 0; k < n; k++) {
    let p = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[p][k])) p = i;
    if (Math.abs(M[p][k]) < 1e-300) return null;
    [M[k], M[p]] = [M[p], M[k]];
    for (let i = k + 1; i < n; i++) {
      const f = M[i][k] / M[k][k];
      for (let j = k; j <= n; j++) M[i][j] -= f * M[k][j];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = M[i][n];
    for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j];
    x[i] = s / M[i][i];
  }
  return x;
}

/** Central-difference Jacobian of a vector function (rows: outputs, columns: inputs). */
function jacobian(f: (p: number[]) => number[], p: number[], h = 1e-6): number[][] {
  const cols = p.map((_, j) => {
    const a = p.slice(), b = p.slice();
    a[j] += h; b[j] -= h;
    const fa = f(a), fb = f(b);
    return fa.map((v, i) => (v - fb[i]) / (2 * h));
  });
  return cols[0].map((_, i) => cols.map(c => c[i]));
}

/** Levenberg-Marquardt on residuals r(p) = y - model(p); returns p and the SSE. */
export function lm(model: (p: number[]) => number[], y: number[], p0: number[], maxIt = 100): { p: number[]; S: number } {
  const res = (p: number[]) => { const m = model(p); return y.map((v, i) => v - m[i]); };
  let p = p0.slice(), r = res(p), S = sse(r), lam = 1e-3;
  for (let it = 0; it < maxIt; it++) {
    const J = jacobian(model, p);                         // d model / d p
    const k = p.length;
    const A = Array.from({ length: k }, (_, a) => Array.from({ length: k }, (_, b) => J.reduce((s, row) => s + row[a] * row[b], 0)));
    const g = Array.from({ length: k }, (_, a) => J.reduce((s, row, i) => s + row[a] * r[i], 0));
    let improved = false;
    for (let tries = 0; tries < 12; tries++) {
      const B = A.map((row, a) => row.map((v, b) => a === b ? v * (1 + lam) + 1e-12 : v));
      const d = solve(B, g);
      if (!d || d.some(v => !Number.isFinite(v))) { lam *= 10; continue; }
      const step = d.map(v => Math.max(-2, Math.min(2, v)));      // log-parameters: at most e^2
      const pn = p.map((v, j) => v + step[j]);
      const rn = res(pn), Sn = sse(rn);
      if (Number.isFinite(Sn) && Sn <= S) {
        const done = S - Sn <= 1e-14 * (1 + S) && Math.max(...step.map(Math.abs)) < 1e-9;
        p = pn; r = rn; S = Sn; lam = Math.max(lam / 10, 1e-9); improved = true;
        if (done) return { p, S };
        break;
      }
      lam *= 10;
    }
    if (!improved) break;
  }
  return { p, S };
}

/** Model at sample times ts for log-parameters (log ka, log ke, log V). */
export const modelLog = (ts: number[]) => (lp: number[]): number[] =>
  ts.map(t => pk(t, Math.exp(lp[0]), Math.exp(lp[1]), Math.exp(lp[2])));

/** One-sided Jacobi SVD of an n x p matrix (n >= p): singular values (descending), V. */
export function svd(J: number[][]): { s: number[]; V: number[][] } {
  const n = J.length, p = J[0].length;
  const U = J.map(r => r.slice());
  const V = Array.from({ length: p }, (_, i) => Array.from({ length: p }, (_, j) => i === j ? 1 : 0));
  for (let sweep = 0; sweep < 60; sweep++) {
    let off = 0;
    for (let a = 0; a < p - 1; a++) for (let b = a + 1; b < p; b++) {
      let al = 0, be = 0, ga = 0;
      for (let i = 0; i < n; i++) { al += U[i][a] ** 2; be += U[i][b] ** 2; ga += U[i][a] * U[i][b]; }
      if (Math.abs(ga) <= 1e-15 * Math.sqrt(al * be) || ga === 0) continue;
      off = Math.max(off, Math.abs(ga) / Math.sqrt(al * be));
      const z = (be - al) / (2 * ga);
      const tn = Math.sign(z || 1) / (Math.abs(z) + Math.sqrt(1 + z * z));
      const c = 1 / Math.sqrt(1 + tn * tn), s = c * tn;
      for (let i = 0; i < n; i++) { const x = U[i][a], y = U[i][b]; U[i][a] = c * x - s * y; U[i][b] = s * x + c * y; }
      for (let i = 0; i < p; i++) { const x = V[i][a], y = V[i][b]; V[i][a] = c * x - s * y; V[i][b] = s * x + c * y; }
    }
    if (off < 1e-15) break;
  }
  const s = Array.from({ length: p }, (_, j) => Math.sqrt(U.reduce((acc, r) => acc + r[j] * r[j], 0)));
  const order = s.map((_, j) => j).sort((a, b) => s[b] - s[a]);
  return { s: order.map(j => s[j]), V: V.map(row => order.map(j => row[j])) };
}

export interface Analysis {
  th: number[];          // estimate (ka, ke, V)
  S: number;             // SSE at the estimate
  s2: number;            // variance used for SEs and Delta: s^2 = S/(n-3), or SIGMA^2 if noise-free
  sEstimated: boolean;
  sv: number[];          // singular values of the log-parameter Jacobian
  V: number[][];         // right singular vectors (columns)
  cov: number[][];       // covariance of (ka, ke, V), natural units
  se: number[];
  corr: number[][];
  singular: boolean;
}

/** Fit the model to (ts, ys) and compute the linearised uncertainty. */
export function analyse(ts: number[], ys: number[], noisy: boolean): Analysis {
  const f = lm(modelLog(ts), ys, TRUTH.map(Math.log));
  const th = f.p.map(Math.exp);
  const n = ts.length;
  const sEstimated = noisy && n > 3;
  const s2 = sEstimated ? f.S / (n - 3) : SIGMA * SIGMA;
  const J = jacobian(modelLog(ts), f.p);                       // n x 3, log-parameters
  const { s: sv, V } = svd(J);
  const singular = !(sv[2] > 1e-9 * sv[0]);
  const covLog = V.map((_, a) => V.map((_, b) => s2 * sv.reduce((acc, si, k) => acc + V[a][k] * V[b][k] / (si * si), 0)));
  const cov = covLog.map((row, a) => row.map((v, b) => v * th[a] * th[b]));
  const se = cov.map((row, a) => singular ? Infinity : Math.sqrt(row[a]));
  const corr = cov.map((row, a) => row.map((v, b) => singular ? NaN : v / (se[a] * se[b])));
  return { th, S: f.S, s2, sEstimated, sv, V, cov, se, corr, singular };
}

/** min over V of the SSE at fixed (ka, ke) (closed form, A = D/V). */
export function sseKaKe(ts: number[], ys: number[], ka: number, ke: number): number {
  let sgy = 0, sgg = 0;
  const g = ts.map(t => shape(t, ka, ke));
  for (let i = 0; i < ts.length; i++) { sgy += g[i] * ys[i]; sgg += g[i] * g[i]; }
  const a = sgg > 0 ? Math.max(0, sgy / sgg) : 0;
  let s = 0;
  for (let i = 0; i < ts.length; i++) { const r = ys[i] - a * g[i]; s += r * r; }
  return s;
}

/** Profile of ka: min over (ke, V) of the SSE with ka fixed. A coarse scan of ke (V in
 *  closed form) picks the start, then LM on (log ke, log V). */
export function profileKa(ts: number[], ys: number[], ka: number): number {
  let best = Infinity, lkeBest = Math.log(0.2);
  for (let i = 0; i <= 50; i++) {
    const lke = Math.log(0.005) + i * (Math.log(20) - Math.log(0.005)) / 50;
    const v = sseKaKe(ts, ys, ka, Math.exp(lke));
    if (v < best) { best = v; lkeBest = lke; }
  }
  // V for the start: closed form at the best ke
  const g = ts.map(t => shape(t, ka, Math.exp(lkeBest)));
  const a = g.reduce((s, v, i) => s + v * ys[i], 0) / Math.max(1e-300, g.reduce((s, v) => s + v * v, 0));
  const lV0 = a > 0 ? Math.log(DOSE / a) : Math.log(20);
  const model = (q: number[]) => ts.map(t => pk(t, ka, Math.exp(q[0]), Math.exp(q[1])));
  const f = lm(model, ys, [lkeBest, lV0], 40);
  return Math.min(f.S, best);
}

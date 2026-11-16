// Numerics for the gnlm demo: the oral one-compartment PK model (Bateman function),
// a forward-difference Jacobian, Householder QR least squares, and single passes of
// Gauss-Newton and Levenberg-Marquardt exactly as in the lecture's R code
// (nonlinear-least-squares/lec.qmd): GN solves J delta = r by QR; LM solves the stacked
// problem [J; sqrt(lambda) Dg] delta = [r; 0] with Dg = diag(column lengths of J),
// accepts if the RSS goes down (lambda / 10), otherwise rejects (lambda * 10).

export const DOSE = 100;
export const TT = [0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12];
export const Y = [1.02, 1.88, 2.34, 2.98, 3.27, 3.26, 3.17, 2.77, 1.92, 1.23, 1.00, 0.61];
export const P_HAT = [1.0150721, 0.1906157, 20.2162104];
export const P_FLIP = [0.1906157, 1.0150721, 3.7963083];
export const RSS_HAT = 0.06646744;

/** Bateman function as written on the slides (0/0 at ka = ke, as in R). */
export function bateman(t: number, ka: number, ke: number, V: number): number {
  return DOSE * ka / (V * (ka - ke)) * (Math.exp(-ke * t) - Math.exp(-ka * t));
}

/** Shape g(t) = ka/(ka-ke) (e^{-ke t} - e^{-ka t}), stable near ka = ke (for the landscape). */
export function shape(t: number, ka: number, ke: number): number {
  const d = ka - ke;
  if (Math.abs(d * t) < 1e-12) return ka * t * Math.exp(-ke * t);
  return ka * Math.exp(-ke * t) * (-Math.expm1(-d * t)) / d;
}

/** RSS with V chosen optimally for this (ka, ke) ("profiled"). Returns [rss, V]. */
export function rssProfile(ka: number, ke: number): [number, number] {
  let sgy = 0, sgg = 0;
  const g = TT.map(t => shape(t, ka, ke));
  for (let i = 0; i < TT.length; i++) { sgy += g[i] * Y[i]; sgg += g[i] * g[i]; }
  const a = sgy / sgg;
  let s = 0;
  for (let i = 0; i < TT.length; i++) { const r = Y[i] - a * g[i]; s += r * r; }
  return [s, DOSE / a];
}

export type Param = "lin" | "log";

/** Natural parameters (ka, ke, V) from the working parameters. */
export const toNat = (x: number[], par: Param): number[] => par === "log" ? x.map(Math.exp) : x.slice();
export const fromNat = (p: number[], par: Param): number[] => par === "log" ? p.map(Math.log) : p.slice();

export function model(x: number[], par: Param): number[] {
  const [ka, ke, V] = toNat(x, par);
  return TT.map(t => bateman(t, ka, ke, V));
}

export function rss(x: number[], par: Param): number {
  const m = model(x, par);
  let s = 0;
  for (let i = 0; i < m.length; i++) { const r = Y[i] - m[i]; s += r * r; }
  return s;
}

/** Forward-difference Jacobian; step 1e-7 p_j (lecture's jac) or 1e-7 max(1, |q_j|) in logs. */
export function jacobian(x: number[], par: Param): number[][] {
  const m0 = model(x, par);
  const J = m0.map(() => [0, 0, 0]);
  for (let j = 0; j < 3; j++) {
    const h = par === "lin" ? 1e-7 * x[j] : 1e-7 * Math.max(1, Math.abs(x[j]));
    const xp = x.slice(); xp[j] += h;
    const hj = xp[j] - x[j];
    const m1 = model(xp, par);
    for (let i = 0; i < m0.length; i++) J[i][j] = (m1[i] - m0[i]) / hj;
  }
  return J;
}

/**
 * Least squares A x ~ b by Householder QR (A is m x n, m >= n). Returns null if A is
 * numerically rank deficient (|R_jj| < 1e-7 times the largest |R_ii|, like qr.solve's
 * rank check) or contains non-finite entries.
 */
export function qrSolve(A: number[][], b: number[]): number[] | null {
  const m = A.length, n = A[0].length;
  const R = A.map(r => r.slice()), y = b.slice();
  for (const r of R) for (const v of r) if (!Number.isFinite(v)) return null;
  for (const v of y) if (!Number.isFinite(v)) return null;
  for (let k = 0; k < n; k++) {
    let nrm = 0;
    for (let i = k; i < m; i++) nrm += R[i][k] * R[i][k];
    nrm = Math.sqrt(nrm);
    if (nrm === 0) return null;
    const alpha = R[k][k] > 0 ? -nrm : nrm;
    const v = new Array(m).fill(0);
    for (let i = k; i < m; i++) v[i] = R[i][k];
    v[k] -= alpha;
    let vv = 0;
    for (let i = k; i < m; i++) vv += v[i] * v[i];
    if (vv === 0) continue;
    for (let j = k; j < n; j++) {
      let s = 0;
      for (let i = k; i < m; i++) s += v[i] * R[i][j];
      s = 2 * s / vv;
      for (let i = k; i < m; i++) R[i][j] -= s * v[i];
    }
    let s = 0;
    for (let i = k; i < m; i++) s += v[i] * y[i];
    s = 2 * s / vv;
    for (let i = k; i < m; i++) y[i] -= s * v[i];
  }
  let rmax = 0;
  for (let k = 0; k < n; k++) rmax = Math.max(rmax, Math.abs(R[k][k]));
  for (let k = 0; k < n; k++) if (!(Math.abs(R[k][k]) > 1e-7 * rmax)) return null;
  const x = new Array(n).fill(0);
  for (let k = n - 1; k >= 0; k--) {
    let s = y[k];
    for (let j = k + 1; j < n; j++) s -= R[k][j] * x[j];
    x[k] = s / R[k][k];
  }
  return x;
}

/** 2-norm condition number of a 12 x 3 matrix via the eigenvalues of J^T J (Jacobi). */
export function cond(J: number[][]): number {
  const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const r of J) for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) A[a][b] += r[a] * r[b];
  for (let sweep = 0; sweep < 50; sweep++) {
    let off = 0;
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 3; q++) off += A[p][q] * A[p][q];
    if (off < 1e-30 * (A[0][0] ** 2 + A[1][1] ** 2 + A[2][2] ** 2)) break;
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 3; q++) {
      if (A[p][q] === 0) continue;
      const th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
      const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < 3; k++) {
        const akp = A[k][p], akq = A[k][q];
        A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq;
      }
      for (let k = 0; k < 3; k++) {
        const apk = A[p][k], aqk = A[q][k];
        A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk;
      }
    }
  }
  const ev = [A[0][0], A[1][1], A[2][2]].map(v => Math.max(v, 0));
  return Math.sqrt(Math.max(...ev) / Math.min(...ev));
}

export type Status = "running" | "converged" | "diverged" | "singular" | "maxit";

export interface Run {
  method: "gn" | "lm";
  par: Param;
  x: number[];          // current working parameters
  rss: number;
  lambda: number;       // LM only
  path: number[][];     // accepted points (working parameters), path[0] = start
  rejected: number[][]; // LM rejected trial points
  passes: number;       // loop passes (LM: trials)
  accepted: number;
  nRejected: number;
  last: "" | "accepted" | "rejected" | "taken";
  lastTrialRss: number;
  status: Status;
}

export function newRun(method: "gn" | "lm", par: Param, start: number[], lambda0 = 0.01): Run {
  const x = fromNat(start, par);
  return {
    method, par, x, rss: rss(x, par), lambda: lambda0, path: [x.slice()], rejected: [],
    passes: 0, accepted: 0, nRejected: 0, last: "", lastTrialRss: NaN, status: "running",
  };
}

const norm = (v: number[]): number => Math.sqrt(v.reduce((s, a) => s + a * a, 0));
export const MAXIT = { gn: 30, lm: 200 };
const TOL = 1e-8;

/** One loop pass of the lecture's gauss_newton / lev_mar. Mutates and returns the run. */
export function step(run: Run): Run {
  if (run.status !== "running") return run;
  const { par } = run;
  run.passes++;
  const r = model(run.x, par).map((m, i) => Y[i] - m);
  const J = jacobian(run.x, par);
  if (run.method === "gn") {
    const d = qrSolve(J, r);
    if (!d) { run.status = "singular"; run.passes--; return run; }
    run.x = run.x.map((v, j) => v + d[j]);
    run.rss = rss(run.x, par);
    run.path.push(run.x.slice());
    run.accepted++;
    run.last = "taken";
    run.lastTrialRss = run.rss;
    const nat = toNat(run.x, par);
    if (!Number.isFinite(run.rss) || run.rss > 1e12 || nat.some(v => !Number.isFinite(v))) run.status = "diverged";
    else if (norm(d) < TOL * norm(run.x)) run.status = "converged";
    else if (run.passes >= MAXIT.gn) run.status = "maxit";
    return run;
  }
  // Levenberg-Marquardt
  const colLen = [0, 1, 2].map(j => Math.sqrt(J.reduce((s, row) => s + row[j] * row[j], 0)));
  const sl = Math.sqrt(run.lambda);
  const A = J.concat([[sl * colLen[0], 0, 0], [0, sl * colLen[1], 0], [0, 0, sl * colLen[2]]]);
  const d = qrSolve(A, r.concat([0, 0, 0]));
  if (!d) { run.status = "singular"; run.passes--; return run; }
  const xt = run.x.map((v, j) => v + d[j]);
  const rt = rss(xt, par);
  run.lastTrialRss = rt;
  if (Number.isFinite(rt) && rt < run.rss) {
    run.x = xt; run.rss = rt; run.lambda /= 10;
    run.path.push(xt.slice()); run.accepted++; run.last = "accepted";
    if (norm(d) < TOL * norm(run.x)) run.status = "converged";
  } else {
    run.lambda *= 10; run.rejected.push(xt); run.nRejected++; run.last = "rejected";
    // the rejected step is already below the tolerance: nothing left to gain
    if (norm(d) < TOL * norm(run.x) || run.lambda > 1e16) run.status = "converged";
  }
  if (run.status === "running" && run.passes >= MAXIT.lm) run.status = "maxit";
  return run;
}

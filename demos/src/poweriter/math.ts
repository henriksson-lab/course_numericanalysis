// Small dense linear algebra for the power-iteration demo (n <= 4): eigenvalues from the
// characteristic polynomial, eigenvectors by inverse iteration, and the iteration itself.
export type Mat = number[][];
export interface C { re: number; im: number }

const cadd = (a: C, b: C): C => ({ re: a.re + b.re, im: a.im + b.im });
const csub = (a: C, b: C): C => ({ re: a.re - b.re, im: a.im - b.im });
const cmul = (a: C, b: C): C => ({ re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re });
const cdiv = (a: C, b: C): C => {
  const d = b.re * b.re + b.im * b.im;
  return { re: (a.re * b.re + a.im * b.im) / d, im: (a.im * b.re - a.re * b.im) / d };
};
export const cabs = (a: C): number => Math.hypot(a.re, a.im);

export const matvec = (A: Mat, x: number[]): number[] => A.map(r => r.reduce((s, a, j) => s + a * x[j], 0));
export const dot = (x: number[], y: number[]): number => x.reduce((s, v, i) => s + v * y[i], 0);
export const norm2 = (x: number[]): number => Math.sqrt(dot(x, x));
export const sum = (x: number[]): number => x.reduce((s, v) => s + v, 0);

/** Characteristic polynomial coefficients c[0..n] (c[n] = 1), Faddeev–LeVerrier. */
function charPoly(A: Mat): number[] {
  const n = A.length;
  const c = new Array(n + 1).fill(0);
  c[n] = 1;
  let M: Mat = A.map(r => r.map(() => 0));
  for (let k = 1; k <= n; k++) {
    // M <- A M + c[n-k+1] I
    const AM = A.map((r, i) => r.map((_, j) => r.reduce((s, a, l) => s + a * M[l][j], 0)));
    M = AM.map((r, i) => r.map((v, j) => v + (i === j ? c[n - k + 1] : 0)));
    const AMk = A.map((r, i) => r.reduce((s, a, l) => s + a * M[l][i], 0));
    c[n - k] = -sum(AMk) / k;
  }
  return c;
}

function polyEval(c: number[], z: C): C {
  let p: C = { re: c[c.length - 1], im: 0 };
  for (let i = c.length - 2; i >= 0; i--) p = cadd(cmul(p, z), { re: c[i], im: 0 });
  return p;
}

/** All eigenvalues of a small real matrix, largest modulus first. */
export function eigenvalues(A: Mat): C[] {
  const n = A.length;
  let z: C[];
  if (n === 1) z = [{ re: A[0][0], im: 0 }];
  else if (n === 2) {
    const tr = A[0][0] + A[1][1], det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    const disc = tr * tr / 4 - det;
    if (disc >= 0) {
      const s = Math.sqrt(disc);
      // stable: the larger-magnitude root first, the other from det / root
      const r1 = tr / 2 + (tr >= 0 ? s : -s);
      const r2 = r1 !== 0 ? det / r1 : tr / 2 - s;
      z = [{ re: r1, im: 0 }, { re: r2, im: 0 }];
    } else {
      const s = Math.sqrt(-disc);
      z = [{ re: tr / 2, im: s }, { re: tr / 2, im: -s }];
    }
  } else {
    const c = charPoly(A);
    const R = 1 + Math.max(...c.slice(0, n).map(Math.abs));
    z = Array.from({ length: n }, (_, i) => {
      const a = 2 * Math.PI * i / n + 0.4;
      return { re: 0.5 * R * Math.cos(a), im: 0.5 * R * Math.sin(a) };
    });
    for (let it = 0; it < 2000; it++) {   // Durand–Kerner
      let move = 0;
      for (let i = 0; i < n; i++) {
        let den: C = { re: 1, im: 0 };
        for (let j = 0; j < n; j++) if (j !== i) den = cmul(den, csub(z[i], z[j]));
        const d = cdiv(polyEval(c, z[i]), den);
        z[i] = csub(z[i], d);
        move = Math.max(move, cabs(d));
      }
      if (move < 1e-15 * R) break;
    }
    z = z.map(v => Math.abs(v.im) < 1e-9 * (1 + Math.abs(v.re)) ? { re: v.re, im: 0 } : v);
  }
  return z.sort((a, b) => cabs(b) - cabs(a) || b.re - a.re || b.im - a.im);
}

/** Solve M y = b by Gaussian elimination with partial pivoting. Exact zero pivots are
 *  replaced by a tiny number, as a shift exactly at an eigenvalue still gives a direction. */
export function solve(M: Mat, b: number[]): { y: number[]; tiny: boolean } {
  const n = M.length;
  const A = M.map((r, i) => [...r, b[i]]);
  const scl = Math.max(1e-300, ...M.flat().map(Math.abs));
  let tiny = false;
  for (let k = 0; k < n; k++) {
    let p = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(A[i][k]) > Math.abs(A[p][k])) p = i;
    [A[k], A[p]] = [A[p], A[k]];
    if (Math.abs(A[k][k]) < 1e-14 * scl) tiny = true;
    if (A[k][k] === 0) A[k][k] = 1e-17 * scl;
    for (let i = k + 1; i < n; i++) {
      const f = A[i][k] / A[k][k];
      for (let j = k; j <= n; j++) A[i][j] -= f * A[k][j];
    }
  }
  const y = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = A[i][n];
    for (let j = i + 1; j < n; j++) s -= A[i][j] * y[j];
    y[i] = s / A[i][i];
  }
  return { y, tiny };
}

const shifted = (A: Mat, s: number): Mat => A.map((r, i) => r.map((v, j) => v - (i === j ? s : 0)));

/** Eigenvector for a real eigenvalue lam (approximate), refined by inverse iteration;
 *  returns the vector (unit 2-norm) and the refined eigenvalue (Rayleigh quotient). */
export function eigvec(A: Mat, lam: number): { v: number[]; lam: number } {
  const n = A.length;
  const M = shifted(A, lam + 1e-10 * Math.max(1, Math.abs(lam)));
  let x = Array.from({ length: n }, (_, i) => 1 - 0.13 * i);
  for (let it = 0; it < 4; it++) {
    const { y } = solve(M, x);
    const s = norm2(y);
    x = y.map(v => v / s);
  }
  // sign: largest entry positive
  const big = x.reduce((m, v, i) => Math.abs(v) > Math.abs(x[m]) ? i : m, 0);
  if (x[big] < 0) x = x.map(v => -v);
  return { v: x, lam: dot(x, matvec(A, x)) };
}

export interface Target { lam: number; v: number[] }
export interface Analysis {
  eig: C[];
  target: Target | null;   // the eigenpair the iteration should converge to
  rate: number;            // predicted error ratio per step
  why: string;             // reason when there is no target
}

/** What power iteration (sigma === null) or inverse iteration with shift sigma converges to. */
export function analyse(A: Mat, sigma: number | null): Analysis {
  const eig = eigenvalues(A);
  const n = eig.length;
  const tol = 1e-9;
  if (sigma === null) {
    const m1 = cabs(eig[0]);
    const ties = eig.filter(e => cabs(e) > m1 * (1 - tol) - 1e-300);
    const rate = ties.length > 1 ? 1 : n > 1 ? cabs(eig[1]) / m1 : 0;
    const pick = ties.find(e => e.im === 0 && e.re > 0) ?? ties.find(e => e.im === 0);
    if (!pick || m1 === 0) return { eig, target: null, rate, why: "the largest eigenvalues are a complex pair: x keeps rotating" };
    const t = eigvec(A, pick.re);
    return { eig, target: t, rate, why: ties.length > 1 ? `${ties.length} eigenvalues share the largest modulus ${m1.toFixed(4)}: no dominant one` : "" };
  }
  const d = eig.map(e => cabs(csub(e, { re: sigma, im: 0 })));
  const order = d.map((_, i) => i).sort((a, b) => d[a] - d[b]);
  const d0 = d[order[0]], d1 = n > 1 ? d[order[1]] : Infinity;
  const rate = d1 === Infinity ? 0 : d0 / d1;
  const e0 = eig[order[0]];
  if (n > 1 && d1 - d0 <= tol * Math.max(1e-12, d1))
    return { eig, target: null, rate: 1, why: "two eigenvalues are equally near the shift: no single nearest one" };
  if (e0.im !== 0) return { eig, target: null, rate, why: "the eigenvalue nearest the shift is complex: x keeps rotating" };
  return { eig, target: eigvec(A, e0.re), rate, why: "" };
}

export interface Run {
  xs: number[][];      // x_0 .. x_K (normalised)
  lam: number[];       // estimate at step k (NaN at k = 0)
  res: number[];       // residual ||A x - lam x|| at step k
  useSum: boolean;
  tiny: boolean;       // A - sigma I numerically singular
}

/** Power iteration (sigma null) or inverse iteration, K steps. Sum normalisation when A
 *  and x0 are nonnegative and this is power iteration; otherwise the 2-norm. */
export function iterate(A: Mat, x0: number[], sigma: number | null, K: number): Run {
  const nonneg = A.every(r => r.every(v => v >= 0)) && x0.every(v => v >= 0) && sum(x0) > 0;
  const useSum = sigma === null && nonneg;
  const nrm = (y: number[]) => useSum ? sum(y) : norm2(y);
  const s0 = nrm(x0);
  const xs = [x0.map(v => v / s0)], lam = [NaN], res = [NaN];
  let tiny = false;
  const M = sigma === null ? null : shifted(A, sigma);
  for (let k = 1; k <= K; k++) {
    const x = xs[k - 1];
    if (!x.every(Number.isFinite)) { xs.push(x); lam.push(NaN); res.push(NaN); continue; }
    if (M === null) {
      const y = matvec(A, x);
      const l = useSum ? sum(y) : dot(x, y) / dot(x, x);
      res.push(norm2(y.map((v, i) => v - l * x[i])));
      lam.push(l);
      const s = nrm(y);
      xs.push(s !== 0 ? y.map(v => v / s) : y.map(() => NaN));
    } else {
      const r = solve(M, x);
      tiny = tiny || r.tiny;
      const s = norm2(r.y);
      const xn = r.y.map(v => v / s);
      const Ax = matvec(A, xn);
      const l = dot(xn, Ax);
      xs.push(xn); lam.push(l);
      res.push(norm2(Ax.map((v, i) => v - l * xn[i])));
    }
  }
  return { xs, lam, res, useSum, tiny };
}

/** Error of x against the target eigenvector: in proportions when the eigenvector is
 *  positive, otherwise as unit vectors with the better sign. Max-abs difference. */
export function vecErr(x: number[], v: number[]): number {
  const pos = v.every(c => c >= 0) || v.every(c => c <= 0);
  const sx = sum(x);
  if (pos && Math.abs(sx) > 1e-300) {
    const sv = sum(v);
    return Math.max(...x.map((c, i) => Math.abs(c / sx - v[i] / sv)));
  }
  const nx = norm2(x), nv = norm2(v);
  const s = dot(x, v) >= 0 ? 1 : -1;
  return Math.max(...x.map((c, i) => Math.abs(c / nx - s * v[i] / nv)));
}

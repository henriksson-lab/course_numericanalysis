// Numerics for the qralgo demo: the (shifted) QR algorithm with deflation, as in exercise
// C2 of eigen-algorithms, and a general eigenvalue solver (Hessenberg + Francis double
// shift, after Numerical Recipes' elmhes/hqr) for the sensitivity tab.

export type Mat = number[][];
export interface Cx { re: number; im: number }
export type Shift = "none" | "rayleigh" | "wilkinson";

export const copy = (A: Mat): Mat => A.map(r => [...r]);
export const eye = (n: number): Mat => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
export const cabs = (z: Cx): number => Math.hypot(z.re, z.im);

/** Householder QR of a square matrix: M = Q R. */
export function qr(M: Mat): { Q: Mat; R: Mat } {
  const m = M.length;
  const R = copy(M), Q = eye(m);
  for (let j = 0; j < m - 1; j++) {
    let nx = 0;
    for (let i = j; i < m; i++) nx = Math.hypot(nx, R[i][j]);
    if (nx === 0) continue;
    const alpha = R[j][j] >= 0 ? -nx : nx;
    const v = new Array(m).fill(0);
    for (let i = j; i < m; i++) v[i] = R[i][j];
    v[j] -= alpha;
    let vv = 0;
    for (let i = j; i < m; i++) vv += v[i] * v[i];
    if (vv === 0) continue;
    for (let c = 0; c < m; c++) {                  // R <- H R
      let s = 0;
      for (let i = j; i < m; i++) s += v[i] * R[i][c];
      s = 2 * s / vv;
      for (let i = j; i < m; i++) R[i][c] -= s * v[i];
    }
    for (let r = 0; r < m; r++) {                  // Q <- Q H
      let s = 0;
      for (let i = j; i < m; i++) s += Q[r][i] * v[i];
      s = 2 * s / vv;
      for (let i = j; i < m; i++) Q[r][i] -= s * v[i];
    }
    for (let i = j + 1; i < m; i++) R[i][j] = 0;
  }
  return { Q, R };
}

/** Shift from the trailing 2 x 2 of the active block (size m). `complex` if that block has a complex pair. */
export function shiftOf(A: Mat, m: number, type: Shift): { mu: number; complex: boolean } {
  if (type === "none" || m < 2) return { mu: 0, complex: false };
  const c = A[m - 1][m - 1];
  if (type === "rayleigh") return { mu: c, complex: false };
  const a = A[m - 2][m - 2], p = A[m - 2][m - 1] * A[m - 1][m - 2];
  const d = (a - c) / 2, disc = d * d + p;
  if (disc < 0) return { mu: (a + c) / 2, complex: true };   // no real eigenvalue: use the real part
  const den = Math.abs(d) + Math.sqrt(disc);
  return { mu: den === 0 ? c : c - (d >= 0 ? 1 : -1) * p / den, complex: false };
}

export interface Run {
  mats: Mat[];          // A_0 .. A_K
  active: number[];     // size of the active block after step k
  mus: number[];        // shift used in step k (index k-1)
  defK: number[];       // step at which subdiagonal i (a_{i+2,i+1}) was deflated, or -1
  doneK: number;        // step at which the active block reached size 1, or -1
}

/**
 * The QR algorithm with deflation, exactly as qr_eig in exercise C2: the active leading block
 * A_k - mu I = Q R, A_{k+1} = R Q + mu I; once |a_{m,m-1}| < tol the last row/column is
 * removed. The rest of the matrix is transformed with the same Q, so A_k stays similar to A.
 */
export function runQR(A0: Mat, type: Shift, K: number, tol = 1e-12): Run {
  const n = A0.length;
  let A = copy(A0), m = n;
  const run: Run = { mats: [copy(A)], active: [n], mus: [], defK: new Array(Math.max(0, n - 1)).fill(-1), doneK: n <= 1 ? 0 : -1 };
  for (let k = 1; k <= K; k++) {
    if (m <= 1) break;
    const { mu } = shiftOf(A, m, type);
    const M: Mat = [];
    for (let i = 0; i < m; i++) M.push(A[i].slice(0, m).map((v, j) => v - (i === j ? mu : 0)));
    const { Q, R } = qr(M);
    const B = copy(A);
    for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) {
      let s = 0;
      for (let l = i; l < m; l++) s += R[i][l] * Q[l][j];
      B[i][j] = s + (i === j ? mu : 0);
    }
    for (let i = 0; i < m; i++) for (let j = m; j < n; j++) {    // Q^T A[active, rest]
      let s = 0;
      for (let l = 0; l < m; l++) s += Q[l][i] * A[l][j];
      B[i][j] = s;
    }
    for (let i = m; i < n; i++) for (let j = 0; j < m; j++) {    // A[rest, active] Q
      let s = 0;
      for (let l = 0; l < m; l++) s += A[i][l] * Q[l][j];
      B[i][j] = s;
    }
    A = B;
    run.mus.push(mu);
    if (Math.abs(A[m - 1][m - 2]) < tol) {
      run.defK[m - 2] = k;
      m--;
      if (m <= 1) run.doneK = k;
    }
    run.mats.push(copy(A));
    run.active.push(m);
  }
  return run;
}

/** All eigenvalues of a real square matrix: Hessenberg reduction + Francis double-shift QR. */
export function eigvals(A0: Mat): Cx[] {
  const n = A0.length;
  if (n === 0) return [];
  // 1-based copy, as in Numerical Recipes
  const a: number[][] = Array.from({ length: n + 1 }, () => new Array(n + 1).fill(0));
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) a[i + 1][j + 1] = A0[i][j];
  // elmhes: reduction to Hessenberg form by elimination with pivoting
  for (let m = 2; m < n; m++) {
    let x = 0, i = m;
    for (let j = m; j <= n; j++) if (Math.abs(a[j][m - 1]) > Math.abs(x)) { x = a[j][m - 1]; i = j; }
    if (i !== m) {
      for (let j = m - 1; j <= n; j++) { const t = a[i][j]; a[i][j] = a[m][j]; a[m][j] = t; }
      for (let j = 1; j <= n; j++) { const t = a[j][i]; a[j][i] = a[j][m]; a[j][m] = t; }
    }
    if (x !== 0) {
      for (i = m + 1; i <= n; i++) {
        let y = a[i][m - 1];
        if (y !== 0) {
          y /= x;
          a[i][m - 1] = y;
          for (let j = m; j <= n; j++) a[i][j] -= y * a[m][j];
          for (let j = 1; j <= n; j++) a[j][m] += y * a[j][i];
        }
      }
    }
  }
  for (let i = 3; i <= n; i++) for (let j = 1; j <= i - 2; j++) a[i][j] = 0;
  // hqr
  const wr = new Array(n + 1).fill(NaN), wi = new Array(n + 1).fill(0);
  let anorm = 0;
  for (let i = 1; i <= n; i++) for (let j = Math.max(i - 1, 1); j <= n; j++) anorm += Math.abs(a[i][j]);
  let nn = n, t = 0;
  const sign = (u: number, v: number) => (v >= 0 ? Math.abs(u) : -Math.abs(u));
  outer: while (nn >= 1) {
    let its = 0, l = 1;
    do {
      for (l = nn; l >= 2; l--) {
        let s = Math.abs(a[l - 1][l - 1]) + Math.abs(a[l][l]);
        if (s === 0) s = anorm;
        if (Math.abs(a[l][l - 1]) + s === s) { a[l][l - 1] = 0; break; }
      }
      let x = a[nn][nn];
      if (l === nn) {
        wr[nn] = x + t; wi[nn--] = 0;
      } else {
        let y = a[nn - 1][nn - 1];
        let w = a[nn][nn - 1] * a[nn - 1][nn];
        if (l === nn - 1) {
          const p = 0.5 * (y - x), q = p * p + w;
          let z = Math.sqrt(Math.abs(q));
          x += t;
          if (q >= 0) {
            z = p + sign(z, p);
            wr[nn - 1] = wr[nn] = x + z;
            if (z) wr[nn] = x - w / z;
            wi[nn - 1] = wi[nn] = 0;
          } else {
            wr[nn - 1] = wr[nn] = x + p;
            wi[nn] = z; wi[nn - 1] = -z;
          }
          nn -= 2;
        } else {
          if (its === 60) break outer;
          if (its === 10 || its === 20 || its === 40) {
            t += x;
            for (let i = 1; i <= nn; i++) a[i][i] -= x;
            const s = Math.abs(a[nn][nn - 1]) + Math.abs(a[nn - 1][nn - 2]);
            y = x = 0.75 * s;
            w = -0.4375 * s * s;
          }
          ++its;
          let m = nn - 2, p = 0, q = 0, r = 0, z = 0;
          for (m = nn - 2; m >= l; m--) {
            z = a[m][m];
            r = x - z;
            let s = y - z;
            p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
            q = a[m + 1][m + 1] - z - r - s;
            r = a[m + 2][m + 1];
            s = Math.abs(p) + Math.abs(q) + Math.abs(r);
            p /= s; q /= s; r /= s;
            if (m === l) break;
            const u = Math.abs(a[m][m - 1]) * (Math.abs(q) + Math.abs(r));
            const v = Math.abs(p) * (Math.abs(a[m - 1][m - 1]) + Math.abs(z) + Math.abs(a[m + 1][m + 1]));
            if (u + v === v) break;
          }
          for (let i = m + 2; i <= nn; i++) {
            a[i][i - 2] = 0;
            if (i !== m + 2) a[i][i - 3] = 0;
          }
          for (let k = m; k <= nn - 1; k++) {
            if (k !== m) {
              p = a[k][k - 1];
              q = a[k + 1][k - 1];
              r = 0;
              if (k !== nn - 1) r = a[k + 2][k - 1];
              if ((x = Math.abs(p) + Math.abs(q) + Math.abs(r)) !== 0) { p /= x; q /= x; r /= x; }
            }
            const s = sign(Math.sqrt(p * p + q * q + r * r), p);
            if (s !== 0) {
              if (k === m) {
                if (l !== m) a[k][k - 1] = -a[k][k - 1];
              } else a[k][k - 1] = -s * x;
              p += s;
              x = p / s; y = q / s; z = r / s;
              q /= p; r /= p;
              for (let j = k; j <= nn; j++) {
                p = a[k][j] + q * a[k + 1][j];
                if (k !== nn - 1) { p += r * a[k + 2][j]; a[k + 2][j] -= p * z; }
                a[k + 1][j] -= p * y;
                a[k][j] -= p * x;
              }
              const mmin = nn < k + 3 ? nn : k + 3;
              for (let i = l; i <= mmin; i++) {
                p = x * a[i][k] + y * a[i][k + 1];
                if (k !== nn - 1) { p += z * a[i][k + 2]; a[i][k + 2] -= p * r; }
                a[i][k + 1] -= p * q;
                a[i][k] -= p;
              }
            }
          }
        }
      }
    } while (l < nn - 1);
  }
  const out: Cx[] = [];
  for (let i = 1; i <= n; i++) out.push({ re: wr[i], im: wi[i] });
  return out;
}

/** Eigenvalues sorted by decreasing modulus (a conjugate pair: + imaginary part first). */
export function sortByModulus(ev: Cx[]): Cx[] {
  return [...ev].sort((u, v) => {
    const d = cabs(v) - cabs(u);
    if (Math.abs(d) > 1e-9 * Math.max(1, cabs(u))) return d;
    if (Math.abs(u.re - v.re) > 1e-9) return v.re - u.re;
    return v.im - u.im;
  });
}

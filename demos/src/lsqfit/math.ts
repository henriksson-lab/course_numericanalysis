// Numerics for the lsqfit demo: polynomial least squares in doubles by the normal
// equations (A^T A formed, then Gaussian elimination with partial pivoting, as R's solve
// and numpy.linalg.solve do) and by Householder QR (R x = Q^T b), the singular values of
// A by one-sided Jacobi (for kappa), and the exact least-squares solution in rational
// arithmetic (BigInt), used as the reference for "digits correct".

export const EPS = 2.220446049250313e-16;

/** Design matrix with columns t^0 ... t^d, row-major array of rows. */
export function vander(t: number[], d: number): number[][] {
  return t.map(v => { const r = [1]; for (let k = 1; k <= d; k++) r.push(r[k - 1] * v); return r; });
}

export interface Solve { x: number[] | null; note: string }

/** Normal equations: form A^T A and A^T b in doubles, solve by LU with partial pivoting. */
export function solveNormal(A: number[][], b: number[]): Solve {
  const m = A.length, n = A[0].length;
  const M: number[][] = Array.from({ length: n }, () => new Array(n + 1).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) { let s = 0; for (let k = 0; k < m; k++) s += A[k][i] * A[k][j]; M[i][j] = s; }
    let s = 0; for (let k = 0; k < m; k++) s += A[k][i] * b[k]; M[i][n] = s;
  }
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    if (M[p][c] === 0 || !Number.isFinite(M[p][c])) return { x: null, note: "exactly singular" };
    [M[c], M[p]] = [M[p], M[c]];
    for (let r = c + 1; r < n; r++) {
      const f = M[r][c] / M[c][c];
      for (let j = c; j <= n; j++) M[r][j] -= f * M[c][j];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = M[i][n];
    for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j];
    x[i] = s / M[i][i];
  }
  return { x: x.every(Number.isFinite) ? x : null, note: x.every(Number.isFinite) ? "" : "overflow" };
}

/** Householder QR: reduce [A | b] by reflections, then back substitution R x = (Q^T b)_1..n. */
export function solveQR(A: number[][], b: number[]): Solve {
  const m = A.length, n = A[0].length;
  const W = A.map((r, i) => [...r, b[i]]);
  for (let k = 0; k < n; k++) {
    let nrm = 0;
    for (let i = k; i < m; i++) nrm = Math.hypot(nrm, W[i][k]);
    if (nrm === 0) return { x: null, note: "zero column" };
    const alpha = W[k][k] >= 0 ? -nrm : nrm;          // sign choice avoids cancellation
    const v = new Array(m).fill(0);
    for (let i = k; i < m; i++) v[i] = W[i][k];
    v[k] -= alpha;
    let vv = 0; for (let i = k; i < m; i++) vv += v[i] * v[i];
    if (vv === 0) continue;
    for (let j = k; j <= n; j++) {
      let s = 0; for (let i = k; i < m; i++) s += v[i] * W[i][j];
      const f = 2 * s / vv;
      for (let i = k; i < m; i++) W[i][j] -= f * v[i];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = W[i][n];
    for (let j = i + 1; j < n; j++) s -= W[i][j] * x[j];
    x[i] = s / W[i][i];
  }
  return { x: x.every(Number.isFinite) ? x : null, note: "" };
}

/** Singular values of A (m >= n) by one-sided Jacobi, descending. */
export function singularValues(A: number[][]): number[] {
  const m = A.length, n = A[0].length;
  const U = A.map(r => [...r]);
  for (let sweep = 0; sweep < 60; sweep++) {
    let rot = 0;
    for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
      let a = 0, bb = 0, g = 0;
      for (let i = 0; i < m; i++) { a += U[i][p] * U[i][p]; bb += U[i][q] * U[i][q]; g += U[i][p] * U[i][q]; }
      if (g === 0 || Math.abs(g) <= 1e-15 * Math.sqrt(a * bb)) continue;
      rot++;
      const zeta = (bb - a) / (2 * g);
      const t = Math.sign(zeta || 1) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta));
      const c = 1 / Math.sqrt(1 + t * t), s = c * t;
      for (let i = 0; i < m; i++) {
        const up = U[i][p], uq = U[i][q];
        U[i][p] = c * up - s * uq; U[i][q] = s * up + c * uq;
      }
    }
    if (!rot) break;
  }
  const sv: number[] = [];
  for (let j = 0; j < n; j++) { let s = 0; for (let i = 0; i < m; i++) s = Math.hypot(s, U[i][j]); sv.push(s); }
  return sv.sort((a, b) => b - a);
}

// ---------------------------------------------------------------- exact rationals
export interface Rat { p: bigint; q: bigint }
const babs = (a: bigint) => (a < 0n ? -a : a);
function gcd(a: bigint, b: bigint): bigint { a = babs(a); b = babs(b); while (b) [a, b] = [b, a % b]; return a; }
function norm(r: Rat): Rat {
  if (r.q < 0n) r = { p: -r.p, q: -r.q };
  const g = gcd(r.p, r.q) || 1n;
  return { p: r.p / g, q: r.q / g };
}
function bitlen(a: bigint): number { return a === 0n ? 0 : a.toString(2).length; }

/** Rational to the nearest double (to well below the last bit). */
export function ratToNum(r: Rat): number {
  if (r.p === 0n) return 0;
  const neg = (r.p < 0n) !== (r.q < 0n);
  const p = babs(r.p), q = babs(r.q);
  const sh = 70 - (bitlen(p) - bitlen(q));        // quotient gets about 70 bits
  const quo = sh >= 0 ? (p << BigInt(sh)) / q : p / (q << BigInt(-sh));
  const v = Number(quo) * Math.pow(2, -sh);
  return neg ? -v : v;
}

/**
 * Exact least-squares coefficients for columns u^0..u^d with integer u_i and integer
 * readings y_i (any common scale factors are applied by the caller). Fraction-free
 * (Bareiss) elimination on the integer normal equations: exact, so no conditioning issue.
 */
export function exactLS(u: bigint[], y: bigint[], d: number): Rat[] {
  const n = d + 1;
  const pw = u.map(v => { const r = [1n]; for (let k = 1; k <= 2 * d; k++) r.push(r[k - 1] * v); return r; });
  const M: bigint[][] = Array.from({ length: n }, (_, i) => {
    const row: bigint[] = [];
    for (let j = 0; j < n; j++) { let s = 0n; for (const r of pw) s += r[i + j]; row.push(s); }
    let s = 0n; pw.forEach((r, k) => { s += r[i] * y[k]; }); row.push(s);
    return row;
  });
  let prev = 1n;
  for (let c = 0; c < n; c++) {
    let p = c;
    while (p < n && M[p][c] === 0n) p++;
    if (p === n) throw new Error("singular");
    [M[c], M[p]] = [M[p], M[c]];
    for (let r = c + 1; r < n; r++) {
      for (let j = c + 1; j <= n; j++) M[r][j] = (M[c][c] * M[r][j] - M[r][c] * M[c][j]) / prev;
      M[r][c] = 0n;
    }
    prev = M[c][c];
  }
  const x: Rat[] = new Array(n);
  for (let i = n - 1; i >= 0; i--) {
    // x_i = (M[i][n] - sum M[i][j] x_j) / M[i][i]
    let num: Rat = { p: M[i][n], q: 1n };
    for (let j = i + 1; j < n; j++) num = norm({ p: num.p * x[j].q - M[i][j] * x[j].p * num.q, q: num.q * x[j].q });
    x[i] = norm({ p: num.p, q: num.q * M[i][i] });
  }
  return x;
}

// Numerics for the deconv demo: the svd unit's simulated signature matrix (40 genes x 5
// cell types), its SVD by one-sided Jacobi (Hestenes: rotate pairs of columns of A until
// they are orthogonal; the column norms are then the singular values) and the filtered
// SVD solutions x = sum_i phi_i (u_i^T b) v_i for least squares, truncated SVD and ridge.
// One-sided Jacobi works on A itself, not on A^T A, so sigma_5 = 2.87 is accurate to
// working precision even though kappa^2 = 5e4 (the lecture's warning about A^T A).

export const EPS = 2.220446049250313e-16;
export const M = 40, N = 5;
export const TYPES = ["CD4T", "CD8T", "B", "NK", "Mono"];
export const P_TRUE = [0.30, 0.15, 0.15, 0.10, 0.30];

/** Signature with CD8T = CD4T except gene 1 x (1 - d) and gene 2 x (1 + d); d = 0.02 is the lecture's. */
export function signature(d: number): number[][] {
  const S: number[][] = [];
  for (let g = 1; g <= M; g++) {
    const row = new Array(N).fill(10);
    if (g <= 16) { row[0] = 100; row[1] = 100; }
    else row[2 + Math.floor((g - 17) / 8)] = 100;
    for (let j = 1; j <= N; j++) row[j - 1] *= 1 + 0.5 * Math.sin(g * j);
    row[1] = row[0];
    if (g === 1) row[1] = (1 - d) * row[0];
    if (g === 2) row[1] = (1 + d) * row[0];
    S.push(row);
  }
  return S;
}

export const matVec = (A: number[][], x: number[]): number[] => A.map(r => r.reduce((s, a, j) => s + a * x[j], 0));
export const norm = (x: number[]): number => Math.sqrt(x.reduce((s, v) => s + v * v, 0));

export interface SVD { d: number[]; U: number[][]; V: number[][] }   // U: n columns of length m; V: n columns of length n

/** Thin SVD by one-sided Jacobi. Returns columns u_i, v_i, sorted by decreasing sigma; v_i signed so its largest entry is positive. */
export function svd(A: number[][]): SVD {
  const m = A.length, n = A[0].length;
  const C = Array.from({ length: n }, (_, j) => A.map(r => r[j]));            // columns of A
  const V = Array.from({ length: n }, (_, j) => Array.from({ length: n }, (_, i): number => (i === j ? 1 : 0)));
  for (let sweep = 0; sweep < 60; sweep++) {
    let rotated = false;
    for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
      let al = 0, be = 0, ga = 0;
      for (let k = 0; k < m; k++) { al += C[p][k] ** 2; be += C[q][k] ** 2; ga += C[p][k] * C[q][k]; }
      if (ga === 0 || Math.abs(ga) <= EPS * Math.sqrt(al * be)) continue;
      rotated = true;
      const zeta = (be - al) / (2 * ga);
      const t = Math.sign(zeta || 1) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta));
      const c = 1 / Math.sqrt(1 + t * t), s = c * t;
      for (let k = 0; k < m; k++) { const a = C[p][k], b = C[q][k]; C[p][k] = c * a - s * b; C[q][k] = s * a + c * b; }
      for (let k = 0; k < n; k++) { const a = V[p][k], b = V[q][k]; V[p][k] = c * a - s * b; V[q][k] = s * a + c * b; }
    }
    if (!rotated) break;
  }
  const d = C.map(norm);
  const order = d.map((_, i) => i).sort((a, b) => d[b] - d[a]);
  const out: SVD = { d: [], U: [], V: [] };
  for (const i of order) {
    let v = V[i], u = C[i].map(x => (d[i] > 0 ? x / d[i] : 0));
    const big = v.reduce((bi, x, k) => (Math.abs(x) > Math.abs(v[bi]) ? k : bi), 0);
    if (v[big] < 0) { v = v.map(x => -x); u = u.map(x => -x); }
    out.d.push(d[i]); out.U.push(u); out.V.push(v);
  }
  return out;
}

export type Method = "ls" | "tsvd" | "ridge";

/** Multipliers phi_i of (u_i^T b) in x = sum phi_i (u_i^T b) v_i; filter factor f_i = phi_i sigma_i. */
export function phis(sv: SVD, method: Method, k: number, lambda: number): number[] {
  const tol = Math.max(M, N) * EPS * sv.d[0];                    // pseudoinverse tolerance (rounding level)
  return sv.d.map((s, i) => {
    if (method === "ridge") return s / (s * s + lambda * lambda);
    if (s <= tol) return 0;
    if (method === "tsvd" && i >= k) return 0;
    return 1 / s;
  });
}

/** The filtered pseudoinverse M = sum_i phi_i v_i u_i^T (n x m), so that x = M b. */
export function filterMatrix(sv: SVD, phi: number[]): number[][] {
  const m = sv.U[0].length, n = sv.V.length;
  const Mx = Array.from({ length: n }, () => new Array(m).fill(0));
  for (let i = 0; i < n; i++) {
    if (phi[i] === 0) continue;
    for (let r = 0; r < n; r++) {
      const c = phi[i] * sv.V[i][r];
      for (let g = 0; g < m; g++) Mx[r][g] += c * sv.U[i][g];
    }
  }
  return Mx;
}

/**
 * Expected error of x = M b over noise b = b0 * (1 + s z), z ~ N(0, I):
 * bias = ||M b0 - p|| (no noise), noise = sqrt(E ||M delta b||^2) = sqrt(sum_g (s b0_g)^2 ||M[:, g]||^2),
 * total = sqrt(bias^2 + noise^2) (root mean square error).
 */
export function expectedError(Mx: number[][], b0: number[], p: number[], s: number) {
  const x0 = matVec(Mx, b0);
  const bias = norm(x0.map((v, i) => v - p[i]));
  let v = 0;
  for (let g = 0; g < b0.length; g++) {
    let c = 0;
    for (let r = 0; r < Mx.length; r++) c += Mx[r][g] ** 2;
    v += (s * b0[g]) ** 2 * c;
  }
  const noise = Math.sqrt(v);
  return { bias, noise, total: Math.sqrt(bias * bias + noise * noise) };
}

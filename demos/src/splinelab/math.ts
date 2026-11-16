// Numerics for the spline lab: piecewise linear, cubic splines in second-derivative form
// (natural and not-a-knot, both solved as tridiagonal systems by the Thomas algorithm),
// a monotone cubic (PCHIP, Fritsch–Butland slopes as in scipy), the interpolating
// polynomial (barycentric form) and the smoothing spline (Reinsch: (R + λQᵀQ)γ = Qᵀy,
// pentadiagonal, banded Cholesky), with df = tr A and GCV.

export interface Curve {
  f(t: number): number;
  d(t: number): number;
}

/** Thomas algorithm: a sub, b diag, c super (a[0], c[n-1] unused). */
export function thomas(a: number[], b: number[], c: number[], d: number[]): number[] {
  const n = b.length;
  const cp = new Array<number>(n), dp = new Array<number>(n), x = new Array<number>(n);
  cp[0] = c[0] / b[0]; dp[0] = d[0] / b[0];
  for (let i = 1; i < n; i++) {
    const m = b[i] - a[i] * cp[i - 1];
    cp[i] = c[i] / m;
    dp[i] = (d[i] - a[i] * dp[i - 1]) / m;
  }
  x[n - 1] = dp[n - 1];
  for (let i = n - 2; i >= 0; i--) x[i] = dp[i] - cp[i] * x[i + 1];
  return x;
}

/** Index i with x[i] <= t < x[i+1], clamped to [0, n-2]. */
function seg(x: number[], t: number): number {
  let lo = 0, hi = x.length - 2;
  if (t <= x[0]) return 0;
  if (t >= x[hi]) return hi;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (x[m] <= t) lo = m; else hi = m; }
  return lo;
}

export function linear(x: number[], y: number[]): Curve {
  return {
    f: t => { const i = seg(x, t), h = x[i + 1] - x[i]; return y[i] + (y[i + 1] - y[i]) * (t - x[i]) / h; },
    d: t => { const i = seg(x, t); return (y[i + 1] - y[i]) / (x[i + 1] - x[i]); },
  };
}

/** Cubic spline from values y and second derivatives M at the knots (Fact 3). */
export function fromM(x: number[], y: number[], M: number[]): Curve {
  return {
    f: t => {
      const i = seg(x, t), h = x[i + 1] - x[i], A = (x[i + 1] - t) / h, B = (t - x[i]) / h;
      return A * y[i] + B * y[i + 1] + h * h / 6 * ((A ** 3 - A) * M[i] + (B ** 3 - B) * M[i + 1]);
    },
    d: t => {
      const i = seg(x, t), h = x[i + 1] - x[i], A = (x[i + 1] - t) / h, B = (t - x[i]) / h;
      return (y[i + 1] - y[i]) / h + h / 6 * (-(3 * A * A - 1) * M[i] + (3 * B * B - 1) * M[i + 1]);
    },
  };
}

/** Interior equations of Fact 4, rows i = 1..n-2 (0-based). */
function interior(x: number[], y: number[]) {
  const n = x.length, h = x.slice(1).map((v, i) => v - x[i]);
  const s = h.map((hi, i) => (y[i + 1] - y[i]) / hi);
  const a: number[] = [], b: number[] = [], c: number[] = [], r: number[] = [];
  for (let i = 1; i < n - 1; i++) {
    a.push(h[i - 1]); b.push(2 * (h[i - 1] + h[i])); c.push(h[i]); r.push(6 * (s[i] - s[i - 1]));
  }
  return { h, a, b, c, r };
}

export function naturalM(x: number[], y: number[]): number[] {
  const n = x.length;
  if (n < 3) return new Array(n).fill(0);
  const { a, b, c, r } = interior(x, y);
  return [0, ...thomas(a, b, c, r), 0];
}

/** Not-a-knot: S''' continuous at x[1] and x[n-2]; M[0], M[n-1] eliminated, still tridiagonal. */
export function notAKnotM(x: number[], y: number[]): number[] {
  const n = x.length;
  if (n < 4) return naturalM(x, y);
  const { h, a, b, c, r } = interior(x, y);
  const m = n - 2;
  // left: M0 = ((h0+h1) M1 - h0 M2) / h1
  const h0 = h[0], h1 = h[1];
  b[0] += h0 * (h0 + h1) / h1;
  c[0] -= h0 * h0 / h1;
  // right: M_{n-1} = ((hl+hk) M_{n-2} - hl M_{n-3}) / hk, hk = h[n-3], hl = h[n-2]
  const hk = h[n - 3], hl = h[n - 2];
  b[m - 1] += hl * (hl + hk) / hk;
  a[m - 1] -= hl * hl / hk;
  let Mi: number[];
  if (m === 2 && n === 4) {
    // 2x2 system after elimination (a single cubic through 4 points)
    const det = b[0] * b[1] - c[0] * a[1];
    Mi = [(r[0] * b[1] - c[0] * r[1]) / det, (b[0] * r[1] - a[1] * r[0]) / det];
  } else Mi = thomas(a, b, c, r);
  const M0 = ((h0 + h1) * Mi[0] - h0 * Mi[1]) / h1;
  const Mn = ((hl + hk) * Mi[m - 1] - hl * Mi[m - 2]) / hk;
  return [M0, ...Mi, Mn];
}

/** Cubic Hermite curve from values and slopes. */
function hermite(x: number[], y: number[], dk: number[]): Curve {
  return {
    f: t => {
      const i = seg(x, t), h = x[i + 1] - x[i], u = (t - x[i]) / h;
      const h00 = (1 + 2 * u) * (1 - u) ** 2, h10 = u * (1 - u) ** 2, h01 = u * u * (3 - 2 * u), h11 = u * u * (u - 1);
      return h00 * y[i] + h10 * h * dk[i] + h01 * y[i + 1] + h11 * h * dk[i + 1];
    },
    d: t => {
      const i = seg(x, t), h = x[i + 1] - x[i], u = (t - x[i]) / h;
      const g00 = 6 * u * u - 6 * u, g10 = 3 * u * u - 4 * u + 1, g01 = -g00, g11 = 3 * u * u - 2 * u;
      return (g00 * y[i] + g01 * y[i + 1]) / h + g10 * dk[i] + g11 * dk[i + 1];
    },
  };
}

/** PCHIP (Fritsch–Butland), as scipy.interpolate.PchipInterpolator. */
export function pchip(x: number[], y: number[]): Curve {
  const n = x.length;
  if (n < 3) return linear(x, y);
  const h = x.slice(1).map((v, i) => v - x[i]);
  const s = h.map((hi, i) => (y[i + 1] - y[i]) / hi);
  const d = new Array<number>(n).fill(0);
  for (let k = 1; k < n - 1; k++) {
    if (s[k - 1] * s[k] <= 0) { d[k] = 0; continue; }
    const w1 = 2 * h[k] + h[k - 1], w2 = h[k] + 2 * h[k - 1];
    d[k] = (w1 + w2) / (w1 / s[k - 1] + w2 / s[k]);
  }
  const edge = (h0: number, h1: number, s0: number, s1: number): number => {
    let e = ((2 * h0 + h1) * s0 - h0 * s1) / (h0 + h1);
    if (Math.sign(e) !== Math.sign(s0)) e = 0;
    else if (Math.sign(s0) !== Math.sign(s1) && Math.abs(e) > 3 * Math.abs(s0)) e = 3 * s0;
    return e;
  };
  d[0] = edge(h[0], h[1], s[0], s[1]);
  d[n - 1] = edge(h[n - 2], h[n - 3], s[n - 2], s[n - 3]);
  return hermite(x, y, d);
}

/** Interpolating polynomial of degree n-1, barycentric form. */
export function polynomial(x: number[], y: number[]): Curve {
  const n = x.length;
  const w = x.map((xj, j) => { let p = 1; for (let k = 0; k < n; k++) if (k !== j) p *= xj - x[k]; return 1 / p; });
  const f = (t: number): number => {
    let N = 0, D = 0;
    for (let j = 0; j < n; j++) {
      const dt = t - x[j];
      if (dt === 0) return y[j];
      const aj = w[j] / dt; N += aj * y[j]; D += aj;
    }
    return N / D;
  };
  const d = (t: number): number => {
    for (let j = 0; j < n; j++) if (Math.abs(t - x[j]) < 1e-9) t += 1e-7;
    const p = f(t);
    let N = 0, D = 0;
    for (let j = 0; j < n; j++) { const dt = t - x[j]; N += w[j] * (y[j] - p) / (dt * dt); D += w[j] / dt; }
    return -N / D;
  };
  return { f, d };
}

// ---------------------------------------------------------------- smoothing spline

/** Banded Cholesky of a symmetric pentadiagonal matrix: diagonals d0, d1, d2 (length m, m-1, m-2). */
function cholBand(d0: number[], d1: number[], d2: number[]) {
  const m = d0.length;
  const L0 = new Array<number>(m).fill(0), L1 = new Array<number>(m).fill(0), L2 = new Array<number>(m).fill(0);
  // L[i][i] = L0[i], L[i][i-1] = L1[i], L[i][i-2] = L2[i]
  for (let i = 0; i < m; i++) {
    if (i >= 2) L2[i] = d2[i - 2] / L0[i - 2];
    if (i >= 1) L1[i] = (d1[i - 1] - (i >= 2 ? L2[i] * L1[i - 1] : 0)) / L0[i - 1];
    L0[i] = Math.sqrt(d0[i] - L1[i] * L1[i] - L2[i] * L2[i]);
  }
  const solve = (r: number[]): number[] => {
    const z = new Array<number>(m);
    for (let i = 0; i < m; i++) z[i] = (r[i] - (i >= 1 ? L1[i] * z[i - 1] : 0) - (i >= 2 ? L2[i] * z[i - 2] : 0)) / L0[i];
    const g = new Array<number>(m);
    for (let i = m - 1; i >= 0; i--) g[i] = (z[i] - (i + 1 < m ? L1[i + 1] * g[i + 1] : 0) - (i + 2 < m ? L2[i + 2] * g[i + 2] : 0)) / L0[i];
    return g;
  };
  return solve;
}

export interface SmoothFit extends Curve {
  fitted: number[];
  df: number;
  rss: number;
  gcv: number;
}

/**
 * Smoothing spline minimising Σ(y - f)² + lam ∫ f''², lam in the units of x (as Python's
 * make_smoothing_spline). Reinsch form, Green & Silverman: (R + lam QᵀQ) γ = Qᵀy, ŷ = y - lam Q γ.
 */
export function smoothing(x: number[], y: number[], lam: number): SmoothFit {
  const n = x.length;
  const h = x.slice(1).map((v, i) => v - x[i]);
  const m = n - 2;
  // Q is n x m with column j (interior knot j+1): Q[j][j] = 1/h_j, Q[j+1][j] = -1/h_j - 1/h_{j+1}, Q[j+2][j] = 1/h_{j+1}
  const q0 = (j: number) => 1 / h[j], q1 = (j: number) => -1 / h[j] - 1 / h[j + 1], q2 = (j: number) => 1 / h[j + 1];
  const qcol = (j: number): [number, number][] => [[j, q0(j)], [j + 1, q1(j)], [j + 2, q2(j)]];
  const d0 = new Array<number>(m), d1 = new Array<number>(Math.max(0, m - 1)), d2 = new Array<number>(Math.max(0, m - 2));
  const dotQ = (j: number, k: number): number => {
    let s = 0;
    for (const [r1, v1] of qcol(j)) for (const [r2, v2] of qcol(k)) if (r1 === r2) s += v1 * v2;
    return s;
  };
  for (let j = 0; j < m; j++) {
    d0[j] = (h[j] + h[j + 1]) / 3 + lam * dotQ(j, j);
    if (j + 1 < m) d1[j] = h[j + 1] / 6 + lam * dotQ(j, j + 1);
    if (j + 2 < m) d2[j] = lam * dotQ(j, j + 2);
  }
  const solve = cholBand(d0, d1, d2);
  const Qt = (v: number[]): number[] => Array.from({ length: m }, (_, j) => qcol(j).reduce((s, [r, q]) => s + q * v[r], 0));
  const Qv = (g: number[]): number[] => {
    const out = new Array<number>(n).fill(0);
    for (let j = 0; j < m; j++) for (const [r, q] of qcol(j)) out[r] += q * g[j];
    return out;
  };
  const gam = solve(Qt(y));
  const Qg = Qv(gam);
  const fitted = y.map((v, i) => v - lam * Qg[i]);
  // df = tr A, A e_k = e_k - lam Q (R + lam QᵀQ)⁻¹ Qᵀ e_k; Qᵀ e_k has at most three non-zeros
  let df = 0;
  for (let k = 0; k < n; k++) {
    const r = new Array<number>(m).fill(0);
    for (let j = Math.max(0, k - 2); j <= Math.min(m - 1, k); j++) for (const [row, q] of qcol(j)) if (row === k) r[j] += q;
    const g = solve(r);
    let qk = 0;
    for (let j = Math.max(0, k - 2); j <= Math.min(m - 1, k); j++) for (const [row, q] of qcol(j)) if (row === k) qk += q * g[j];
    df += 1 - lam * qk;
  }
  const rss = y.reduce((s, v, i) => s + (v - fitted[i]) ** 2, 0);
  const gcv = n * rss / (n - df) ** 2;
  const curve = fromM(x, fitted, [0, ...gam, 0]);
  return { ...curve, fitted, df, rss, gcv };
}

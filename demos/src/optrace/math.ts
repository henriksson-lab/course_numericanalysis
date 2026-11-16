// Numerical core of the optrace demo: the 2-D objectives and four optimisers, written as
// on the gradient-methods slides (Armijo backtracking with c = 1e-4 and halving, Newton
// with the same line search, BFGS on the inverse Hessian with B0 = I, Nelder-Mead with
// R optim's coefficients, initial simplex and stopping rule).
export type Vec = [number, number];
export type Mat = [number, number, number, number];   // row major [a b; c d]

export interface Problem {
  id: string;
  label: string;        // preset button
  formula: string;      // one-line description
  xlab: string;
  ylab: string;
  dom: [number, number, number, number];   // x0, x1, y0, y1
  start: Vec;
  f: (x: Vec) => number;
  g: (x: Vec) => Vec;
  H: (x: Vec) => Mat;
  hessFD: boolean;      // Hessian by central differences of the analytic gradient
  xstar: Vec;
  gtol: number;
}

// ---- the 4PL dose-response data of the lecture (dose in uM) ----
const DOSE: number[] = [];
for (let i = 0; i <= 8; i++) { const c = 10 ** (-2 + 0.5 * i); DOSE.push(c, c); }
const Y = [12.9, 9.8, 2.9, 0.1, 4.2, 9.7, 17.7, 22.3, 49.7, 42.1, 67.3, 72.6, 97.4, 91.8, 94.6, 93.3, 92.9, 95.1];
export const BOT = 5.393180177899748, TOP = 96.173705843036529;
export const EC50 = 1.270980009586612, HILL = 1.239553930417529;

/** SSE and its gradient in (EC50, hill), EC50 in uM; bottom and top fixed at the best fit. */
function sse4(E: number, h: number): number {
  if (!(E > 0)) return NaN;
  let s = 0;
  for (let i = 0; i < DOSE.length; i++) {
    const m = BOT + (TOP - BOT) / (1 + (E / DOSE[i]) ** h);
    s += (Y[i] - m) ** 2;
  }
  return s;
}
function grad4(E: number, h: number): Vec {
  if (!(E > 0)) return [NaN, NaN];
  let gE = 0, gh = 0;
  for (let i = 0; i < DOSE.length; i++) {
    const u = (E / DOSE[i]) ** h, D = 1 + u;
    const r = Y[i] - (BOT + (TOP - BOT) / D);
    const dE = -(TOP - BOT) * h * u / (E * D * D);
    const dh = -(TOP - BOT) * u * Math.log(E / DOSE[i]) / (D * D);
    gE += -2 * r * dE; gh += -2 * r * dh;
  }
  return [gE, gh];
}

/** Hessian as the central-difference Jacobian of the gradient, symmetrised. */
function fdHess(g: (x: Vec) => Vec, x: Vec): Mat {
  const h0 = 1e-5 * Math.max(1, Math.abs(x[0])), h1 = 1e-5 * Math.max(1, Math.abs(x[1]));
  const a = g([x[0] + h0, x[1]]), b = g([x[0] - h0, x[1]]);
  const c = g([x[0], x[1] + h1]), d = g([x[0], x[1] - h1]);
  const H00 = (a[0] - b[0]) / (2 * h0), H10 = (a[1] - b[1]) / (2 * h0);
  const H01 = (c[0] - d[0]) / (2 * h1), H11 = (c[1] - d[1]) / (2 * h1);
  const off = (H01 + H10) / 2;
  return [H00, off, off, H11];
}

export function quadratic(kappa: number): Problem {
  const p: Problem = {
    id: "quad", label: "elongated bowl",
    formula: `f = ½ (x² + κ y²),  κ = ${fmtK(kappa)}`,
    xlab: "x", ylab: "y", dom: [-2, 2, -1.25, 1.25], start: [-1.6, 1],
    f: x => 0.5 * (x[0] * x[0] + kappa * x[1] * x[1]),
    g: x => [x[0], kappa * x[1]],
    H: () => [1, 0, 0, kappa],
    hessFD: false, xstar: [0, 0], gtol: 1e-8,
  };
  return p;
}
const fmtK = (k: number): string => k >= 10 ? k.toFixed(0) : k.toPrecision(2);

export const ROSEN: Problem = {
  id: "rosen", label: "Rosenbrock banana",
  formula: "f = (1 − x)² + 100 (y − x²)²",
  xlab: "x", ylab: "y", dom: [-2, 2, -1, 3], start: [-1.2, 1],
  f: x => (1 - x[0]) ** 2 + 100 * (x[1] - x[0] * x[0]) ** 2,
  g: x => [-2 * (1 - x[0]) - 400 * x[0] * (x[1] - x[0] * x[0]), 200 * (x[1] - x[0] * x[0])],
  H: x => {
    const a = 2 - 400 * (x[1] - 3 * x[0] * x[0]), b = -400 * x[0];
    return [a, b, b, 200];
  },
  hessFD: false, xstar: [1, 1], gtol: 1e-8,
};

function fpl(id: string, label: string, unit: "uM" | "nM" | "log"): Problem {
  const s = unit === "nM" ? 1000 : 1;   // x = s * EC50(uM)
  const toE = (x: number): number => unit === "log" ? 10 ** x : x / s;
  const p: Problem = {
    id, label,
    formula: unit === "log" ? "4PL: SSE(log₁₀ EC50, hill)" : `4PL: SSE(EC50, hill), EC50 in ${unit === "nM" ? "nM" : "µM"}`,
    xlab: unit === "log" ? "log₁₀ EC50 (EC50 in µM)" : unit === "nM" ? "EC50 (nM)" : "EC50 (µM)",
    ylab: "Hill slope",
    dom: unit === "log" ? [-1.5, 0.75, 0.2, 3.2] : [0.05 * s, 4 * s, 0.2, 3.2],
    start: [unit === "log" ? 0 : s, 1],
    f: x => sse4(toE(x[0]), x[1]),
    g: x => {
      const E = toE(x[0]), [gE, gh] = grad4(E, x[1]);
      const dEdx = unit === "log" ? E * Math.LN10 : 1 / s;
      return [gE * dEdx, gh];
    },
    H: x => fdHess(p.g, x),
    hessFD: true,
    xstar: [unit === "log" ? Math.log10(EC50) : EC50 * s, HILL],
    gtol: 1e-4,
  };
  return p;
}
export const FPL_UM = fpl("ec50", "4PL, EC50 in µM", "uM");
export const FPL_NM = fpl("nm", "4PL, EC50 in nM", "nM");
export const FPL_LOG = fpl("log", "4PL, log₁₀ EC50", "log");

// ---- small linear algebra ----
export const dot = (a: Vec, b: Vec): number => a[0] * b[0] + a[1] * b[1];
export const norm = (a: Vec): number => Math.hypot(a[0], a[1]);
export function eig2(H: Mat): [number, number] {   // symmetric 2x2, ascending
  const tr = H[0] + H[3], det = H[0] * H[3] - H[1] * H[2];
  const d = Math.sqrt(Math.max(0, tr * tr / 4 - det));
  return [tr / 2 - d, tr / 2 + d];
}
function solve2(H: Mat, b: Vec): Vec {
  const det = H[0] * H[3] - H[1] * H[2];
  return [(H[3] * b[0] - H[1] * b[1]) / det, (-H[2] * b[0] + H[0] * b[1]) / det];
}

/** Polish the minimiser with a few Newton steps; returns x* and f*. */
export function polish(p: Problem): { x: Vec; f: number } {
  let x: Vec = [p.xstar[0], p.xstar[1]];
  for (let i = 0; i < 8; i++) {
    const s = solve2(p.H(x), p.g(x));
    const xn: Vec = [x[0] - s[0], x[1] - s[1]];
    if (!(p.f(xn) <= p.f(x))) break;
    x = xn;
  }
  return { x, f: p.f(x) };
}

// ---- optimiser traces ----
export interface Trace {
  x: Vec[];          // iterate k (Nelder-Mead: best vertex)
  fx: number[];      // f at iterate k
  nf: number[];      // cumulative f calls after iteration k
  ng: number[];      // cumulative gradient calls
  nh: number[];      // cumulative Hessian calls
  simplex?: Vec[][]; // Nelder-Mead: the simplex after iteration k
  shifted?: number[]; // Newton: cumulative number of steps with H + mu I
  status: string;
}

const C1 = 1e-4, MAXHALF = 60;
function counter(p: Problem) {
  const c = { nf: 0, ng: 0, nh: 0 };
  return {
    c,
    f: (x: Vec) => { c.nf++; const v = p.f(x); return Number.isNaN(v) ? Infinity : v; },
    g: (x: Vec) => { c.ng++; return p.g(x); },
    H: (x: Vec) => { c.nh++; return p.H(x); },
  };
}
function newTrace(x: Vec, f: number, c: { nf: number; ng: number; nh: number }): Trace {
  return { x: [x], fx: [f], nf: [c.nf], ng: [c.ng], nh: [c.nh], status: "" };
}
function push(t: Trace, x: Vec, f: number, c: { nf: number; ng: number; nh: number }): void {
  t.x.push(x); t.fx.push(f); t.nf.push(c.nf); t.ng.push(c.ng); t.nh.push(c.nh);
}

/** Line search along d from x (Armijo, alpha = 1, 1/2, ...). Returns null on failure. */
function armijo(f: (x: Vec) => number, x: Vec, fx: number, g: Vec, d: Vec): { x: Vec; f: number } | null {
  const slope = dot(g, d);
  let a = 1;
  for (let i = 0; i < MAXHALF; i++) {
    const xn: Vec = [x[0] + a * d[0], x[1] + a * d[1]];
    const fn = f(xn);
    if (fn <= fx + C1 * a * slope) return { x: xn, f: fn };
    a /= 2;
  }
  return null;
}

export function steepest(p: Problem, x0: Vec, maxit = 20000): Trace {
  const k = counter(p);
  let x = x0, fx = k.f(x);
  const t = newTrace(x, fx, k.c);
  for (let it = 0; it < maxit; it++) {
    const g = k.g(x);
    if (!(norm(g) >= p.gtol)) { t.status = Number.isFinite(norm(g)) ? "converged" : "NaN gradient"; break; }
    const r = armijo(k.f, x, fx, g, [-g[0], -g[1]]);
    if (!r) { t.status = "stalled"; break; }
    x = r.x; fx = r.f;
    push(t, x, fx, k.c);
  }
  if (!t.status) t.status = "maxit";
  return t;
}

export function newton(p: Problem, x0: Vec, maxit = 200): Trace {
  const k = counter(p);
  let x = x0, fx = k.f(x);
  const t = newTrace(x, fx, k.c);
  t.shifted = [0];
  let nshift = 0;
  for (let it = 0; it < maxit; it++) {
    const g = k.g(x);
    if (!(norm(g) >= p.gtol)) { t.status = Number.isFinite(norm(g)) ? "converged" : "NaN gradient"; break; }
    let H = k.H(x);
    const [lmin] = eig2(H);
    if (!(lmin > 0)) {     // as in exercise C4(c): make H positive definite
      const mu = 1e-3 * Math.max(Math.abs(H[0]), Math.abs(H[3])) - lmin;
      H = [H[0] + mu, H[1], H[2], H[3] + mu];
      nshift++;
    }
    const s = solve2(H, [-g[0], -g[1]]);
    const r = armijo(k.f, x, fx, g, s);
    if (!r) { t.status = "stalled"; break; }
    x = r.x; fx = r.f;
    push(t, x, fx, k.c);
    t.shifted.push(nshift);
  }
  if (!t.status) t.status = "maxit";
  return t;
}

export function bfgs(p: Problem, x0: Vec, maxit = 2000): Trace {
  const k = counter(p);
  let x = x0, fx = k.f(x), g = k.g(x);
  let B: Mat = [1, 0, 0, 1];
  const t = newTrace(x, fx, k.c);
  for (let it = 0; it < maxit; it++) {
    if (!(norm(g) >= p.gtol)) { t.status = Number.isFinite(norm(g)) ? "converged" : "NaN gradient"; break; }
    const d: Vec = [-(B[0] * g[0] + B[1] * g[1]), -(B[2] * g[0] + B[3] * g[1])];
    const r = armijo(k.f, x, fx, g, d);
    if (!r) { t.status = "stalled"; break; }
    const s: Vec = [r.x[0] - x[0], r.x[1] - x[1]];
    x = r.x; fx = r.f;
    const gn = k.g(x);
    const y: Vec = [gn[0] - g[0], gn[1] - g[1]];
    g = gn;
    const sy = dot(s, y);
    if (sy > 0) {
      // B <- V^T B V + rho s s^T,  V = I - rho y s^T
      const rho = 1 / sy;
      const V: Mat = [1 - rho * y[0] * s[0], -rho * y[0] * s[1], -rho * y[1] * s[0], 1 - rho * y[1] * s[1]];
      const BV: Mat = [B[0] * V[0] + B[1] * V[2], B[0] * V[1] + B[1] * V[3], B[2] * V[0] + B[3] * V[2], B[2] * V[1] + B[3] * V[3]];
      B = [V[0] * BV[0] + V[2] * BV[2] + rho * s[0] * s[0], V[0] * BV[1] + V[2] * BV[3] + rho * s[0] * s[1],
        V[1] * BV[0] + V[3] * BV[2] + rho * s[1] * s[0], V[1] * BV[1] + V[3] * BV[3] + rho * s[1] * s[1]];
    }
    push(t, x, fx, k.c);
  }
  if (!t.status) t.status = "maxit";
  return t;
}

/** Nelder-Mead with R optim's defaults: reflection 1, contraction 0.5, expansion 2,
 *  initial step 0.1 * max|x0|, stop when f(worst) - f(best) <= reltol (|f(x0)| + reltol). */
export function nelderMead(p: Problem, x0: Vec, maxit = 5000, reltol = 1e-8): Trace {
  const k = counter(p);
  const f0 = k.f(x0);
  const step = 0.1 * Math.max(Math.abs(x0[0]), Math.abs(x0[1])) || 0.1;
  let S: Vec[] = [x0, [x0[0] + step, x0[1]], [x0[0], x0[1] + step]];
  let F = [f0, k.f(S[1]), k.f(S[2])];
  const tol = reltol * (Math.abs(f0) + reltol);
  const order = () => {
    const idx = [0, 1, 2].sort((a, b) => F[a] - F[b]);
    S = idx.map(i => S[i]); F = idx.map(i => F[i]);
  };
  order();
  const t = newTrace(S[0], F[0], k.c);
  t.simplex = [S.slice()];
  const lin = (a: Vec, b: Vec, c: number): Vec => [a[0] + c * (b[0] - a[0]), a[1] + c * (b[1] - a[1])];
  for (let it = 0; it < maxit; it++) {
    if (F[2] - F[0] <= tol) { t.status = "converged"; break; }
    const cen: Vec = [(S[0][0] + S[1][0]) / 2, (S[0][1] + S[1][1]) / 2];
    const xr = lin(cen, S[2], -1), fr = k.f(xr);
    if (fr < F[0]) {
      const xe = lin(cen, S[2], -2), fe = k.f(xe);
      if (fe < fr) { S[2] = xe; F[2] = fe; } else { S[2] = xr; F[2] = fr; }
    } else if (fr < F[1]) {
      S[2] = xr; F[2] = fr;
    } else {
      const outside = fr < F[2];
      const xc = outside ? lin(cen, xr, 0.5) : lin(cen, S[2], 0.5);
      const fc = k.f(xc);
      if (fc < (outside ? fr : F[2])) { S[2] = xc; F[2] = fc; }
      else {    // shrink towards the best vertex
        for (let i = 1; i < 3; i++) { S[i] = lin(S[0], S[i], 0.5); F[i] = k.f(S[i]); }
      }
    }
    order();
    push(t, S[0], F[0], k.c);
    t.simplex.push(S.slice());
  }
  if (!t.status) t.status = "maxit";
  return t;
}

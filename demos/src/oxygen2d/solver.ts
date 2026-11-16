// Steady oxygen in the chamber under the lid: D (c_xx + c_yy) = r(c), solved on a grid of
// nodes x_i = i h (i = 0..NX-1), y_j = j h, numbered with x running fastest (k = i + NX j,
// 0-based), by nonlinear SOR: each sweep visits the unknowns in that order and solves the
// five-point equation of one node exactly for its own value,
//     (D/h^2) (4 c - sum of the four neighbours) + r(c) = 0,
// then over-relaxes, c <- c + omega (c_point - c). For Michaelis-Menten uptake the point
// equation is a quadratic with one non-negative root (that is the clipping). Open edges are
// fixed at c0; a sealed edge is a no-flux wall: its nodes are unknowns, with the ghost point
// mirrored (c_{-1} = c_1), as in the lecture's closed trench end.

export type Uptake = "constant" | "linear" | "mm";

export const LX = 4000, LY = 2000;   // µm
export const D = 2000;               // µm^2/s
export const C0 = 200;               // µM
export const KM = 1;                 // µM
export const R_PER_1E5 = 0.5;        // µM/s per 1e5 cells/cm^2
export const HFINE = 12.5;           // µm, the finest grid; cell masks live on it

/** Edges open to air: left (x = 0), right, bottom (y = 0), top. */
export interface Edges { l: boolean; r: boolean; b: boolean; t: boolean }

export interface Grid {
  h: number; NX: number; NY: number;
  c: Float64Array;          // concentration at every node (µM)
  fixed: Uint8Array;        // 1 = boundary value (open edge)
  cells: Float64Array;      // 0/1: cells present at the node
  omega: number;
}

export interface Problem { uptake: Uptake; R0: number; edges: Edges }

export const nodeCount = (h: number) => [Math.round(LX / h) + 1, Math.round(LY / h) + 1];

/** A new grid, with a cell mask sampled from the fine mask and a starting field. */
export function makeGrid(h: number, edges: Edges, fineMask: Uint8Array, start?: Grid): Grid {
  const [NX, NY] = nodeCount(h);
  const [FX] = nodeCount(HFINE);
  const m = Math.round(h / HFINE);
  const c = new Float64Array(NX * NY), fixed = new Uint8Array(NX * NY), cells = new Float64Array(NX * NY);
  for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
    const k = i + NX * j;
    cells[k] = fineMask[i * m + FX * j * m];
    c[k] = start ? sample(start, i * h, j * h) : C0;
  }
  const g: Grid = { h, NX, NY, c, fixed, cells, omega: 1 };
  setEdges(g, edges);
  return g;
}

/** Bilinear interpolation of a grid's field at (x, y). */
export function sample(g: Grid, x: number, y: number): number {
  const fx = Math.min(g.NX - 1.000001, Math.max(0, x / g.h)), fy = Math.min(g.NY - 1.000001, Math.max(0, y / g.h));
  const i = Math.floor(fx), j = Math.floor(fy), u = fx - i, v = fy - j;
  const k = i + g.NX * j;
  return (1 - u) * (1 - v) * g.c[k] + u * (1 - v) * g.c[k + 1] + (1 - u) * v * g.c[k + g.NX] + u * v * g.c[k + g.NX + 1];
}

/** Mark the nodes of open edges as fixed (value c0) and choose omega. */
export function setEdges(g: Grid, e: Edges): void {
  const { NX, NY } = g;
  for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
    const k = i + NX * j;
    const f = (e.l && i === 0) || (e.r && i === NX - 1) || (e.b && j === 0) || (e.t && j === NY - 1);
    g.fixed[k] = f ? 1 : 0;
    if (f) g.c[k] = C0;
  }
  // Optimal SOR parameter of the Laplacian (Jacobi radius from the slowest mode in each
  // direction: half a wavelength between two open edges, a quarter with one sealed edge,
  // none with both sealed).
  const mode = (n: number, a: boolean, b: boolean) => (a && b) ? Math.cos(Math.PI / (n - 1)) : (a || b) ? Math.cos(Math.PI / (2 * (n - 1))) : 1;
  const rhoJ = Math.min(0.99999, (mode(NX, e.l, e.r) + mode(NY, e.b, e.t)) / 2);
  g.omega = 2 / (1 + Math.sqrt(1 - rhoJ * rhoJ));
}

/** One SOR sweep; returns the largest change. */
export function sweep(g: Grid, p: Problem): number {
  const { NX, NY, c, fixed, cells, omega } = g;
  const q = D / (g.h * g.h), a = 4 * q;
  const R0 = p.R0, kLin = p.R0 / C0;
  let dmax = 0;
  for (let j = 0; j < NY; j++) {
    const jd = j > 0 ? -NX : NX, ju = j < NY - 1 ? NX : -NX;
    for (let i = 0; i < NX; i++) {
      const k = i + NX * j;
      if (fixed[k]) continue;
      const il = i > 0 ? -1 : 1, ir = i < NX - 1 ? 1 : -1;
      const s = q * (c[k + il] + c[k + ir] + c[k + jd] + c[k + ju]);
      const R = R0 * cells[k];
      let cp: number;
      if (R === 0) cp = s / a;
      else if (p.uptake === "constant") cp = (s - R) / a;
      else if (p.uptake === "linear") cp = s / (a + kLin * cells[k]);
      else {
        // a c^2 + (a Km + R - s) c - s Km = 0, non-negative root, cancellation-free
        const B = s - a * KM - R, disc = Math.sqrt(B * B + 4 * a * s * KM);
        cp = s <= 0 ? 0 : B >= 0 ? (B + disc) / (2 * a) : 2 * s * KM / (disc - B);
      }
      let cn = c[k] + omega * (cp - c[k]);
      if (p.uptake === "mm" && cn < 0) cn = 0;
      const d = Math.abs(cn - c[k]);
      if (d > dmax) dmax = d;
      c[k] = cn;
    }
  }
  return dmax;
}

/** Largest residual |(D/h^2)(4c - neighbours) + r(c)| over the unknowns (µM/s). */
export function residual(g: Grid, p: Problem): number {
  const { NX, NY, c, fixed, cells } = g;
  const q = D / (g.h * g.h);
  let m = 0;
  for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
    const k = i + NX * j;
    if (fixed[k]) continue;
    const il = i > 0 ? -1 : 1, ir = i < NX - 1 ? 1 : -1, jd = j > 0 ? -NX : NX, ju = j < NY - 1 ? NX : -NX;
    const v = c[k];
    const r = p.R0 * cells[k] * (p.uptake === "constant" ? 1 : p.uptake === "linear" ? v / C0 : v / (KM + v));
    m = Math.max(m, Math.abs(q * (4 * v - c[k + il] - c[k + ir] - c[k + jd] - c[k + ju]) + r));
  }
  return m;
}

export const TOL = 1e-9;   // µM, largest change in a sweep

/** Sweep until converged or out of time. */
export function solve(g: Grid, p: Problem, maxMs = 1e9, maxSweeps = 200000): { sweeps: number; done: boolean } {
  const t0 = performance.now();
  let n = 0;
  while (n < maxSweeps) {
    const d = sweep(g, p);
    n++;
    if (d < TOL) return { sweeps: n, done: true };
    if ((n & 7) === 0 && performance.now() - t0 > maxMs) break;
  }
  return { sweeps: n, done: false };
}

export interface Stats { min: number; kmin: number; hyp: number; hypCells: number; neg: number; nUnk: number }

/** Smallest value and the fraction of unknowns (and of unknowns with cells) below thr. */
export function stats(g: Grid, thr: number): Stats {
  let min = Infinity, kmin = 0, nh = 0, nU = 0, nC = 0, nhc = 0, nn = 0;
  for (let k = 0; k < g.c.length; k++) {
    if (g.fixed[k]) continue;
    nU++;
    const v = g.c[k];
    if (v < min) { min = v; kmin = k; }
    if (v < thr) nh++;
    if (v < 0) nn++;
    if (g.cells[k] > 0) { nC++; if (v < thr) nhc++; }
  }
  return { min, kmin, hyp: nU ? nh / nU : 0, hypCells: nC ? nhc / nC : 0, neg: nU ? nn / nU : 0, nUnk: nU };
}

// Numerical core of the lpgeom demo (no DOM).
//  - 2-variable LPs by geometry: the feasible polygon is the box [0, BIG]^2 clipped by
//    each half-plane a.x <= b; vertices with a coordinate near BIG stand for "infinity".
//    The optimum is the best vertex; a simplex step walks to the better neighbour.
//  - n-variable LPs by a small dense two-phase tableau simplex (Bland's rule), used for
//    the toy FBA network: max c.x subject to rows (<=, >=, =) and x >= 0.

export type Vec = [number, number];
export interface HalfPlane { a: Vec; b: number }          // a.x <= b
export const BIG = 1e4;

const dot = (a: Vec, x: Vec): number => a[0] * x[0] + a[1] * x[1];
export const isFar = (v: Vec): boolean => v[0] > BIG / 2 || v[1] > BIG / 2;

function dedupe(p: Vec[]): Vec[] {
  const out: Vec[] = [];
  for (const v of p) {
    const last = out[out.length - 1];
    if (!last || Math.abs(last[0] - v[0]) + Math.abs(last[1] - v[1]) > 1e-7) out.push(v);
  }
  while (out.length > 1 && Math.abs(out[0][0] - out[out.length - 1][0]) + Math.abs(out[0][1] - out[out.length - 1][1]) <= 1e-7) out.pop();
  return out;
}

/** Sutherland-Hodgman: clip a convex polygon (counter-clockwise) by a.x <= b. */
export function clip(poly: Vec[], h: HalfPlane): Vec[] {
  const out: Vec[] = [];
  const n = poly.length;
  const tol = 1e-9 * (1 + Math.abs(h.b));
  const f = (p: Vec) => dot(h.a, p) - h.b;
  for (let i = 0; i < n; i++) {
    const p = poly[i], q = poly[(i + 1) % n];
    const fp = f(p), fq = f(q);
    if (fp <= tol) out.push(p);
    if (n > 1 && ((fp < -tol && fq > tol) || (fp > tol && fq < -tol))) {
      const t = fp / (fp - fq);
      out.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
    }
  }
  return dedupe(out);
}

/** Feasible polygon of x >= 0 and the half-planes (counter-clockwise, may be empty). */
export function polygon(hs: HalfPlane[]): Vec[] {
  let p: Vec[] = [[0, 0], [BIG, 0], [BIG, BIG], [0, BIG]];
  for (const h of hs) { p = clip(p, h); if (!p.length) break; }
  return p;
}

export type Status = "optimal" | "tie" | "all" | "unbounded" | "infeasible";
export interface Sol2 { status: Status; poly: Vec[]; z: number; opt: number[] }   // opt = vertex indices

export function solve2(hs: HalfPlane[], c: Vec): Sol2 {
  const poly = polygon(hs);
  if (!poly.length) return { status: "infeasible", poly, z: NaN, opt: [] };
  if (c[0] === 0 && c[1] === 0) return { status: "all", poly, z: 0, opt: poly.map((_, i) => i) };
  const zs = poly.map(v => dot(c, v));
  const zmax = Math.max(...zs);
  const tol = 1e-9 * Math.max(1, Math.abs(zmax));
  const opt = zs.map((z, i) => z >= zmax - tol ? i : -1).filter(i => i >= 0);
  const real = opt.filter(i => !isFar(poly[i]));
  if (!real.length) return { status: "unbounded", poly, z: Infinity, opt };
  return { status: opt.length > 1 ? "tie" : "optimal", poly, z: zmax, opt };
}

export function active(h: HalfPlane, v: Vec): boolean {
  return Math.abs(dot(h.a, v) - h.b) <= 1e-7 * (1 + Math.abs(h.b) + Math.abs(v[0]) + Math.abs(v[1]));
}

/** One simplex step on the polygon from vertex index i. `cons` are ALL constraints,
 *  including x1 >= 0 and x2 >= 0 written as -x <= 0. Returns the best improving neighbour
 *  (Dantzig: largest objective gain per unit of the constraint slack that is freed), the
 *  rates of both neighbours, and which constraints are freed along each edge. */
export interface Edge { to: number; rate: number; freed: number[]; along: number[]; far: boolean }
export function edges(poly: Vec[], i: number, cons: HalfPlane[], c: Vec): Edge[] {
  const n = poly.length;
  if (n < 2) return [];
  const nb = n === 2 ? [1 - i] : [(i + n - 1) % n, (i + 1) % n];
  const v = poly[i];
  return nb.map(j => {
    const w = poly[j];
    const d: Vec = [w[0] - v[0], w[1] - v[1]];
    const actV = cons.map((h, k) => active(h, v) ? k : -1).filter(k => k >= 0);
    const along = actV.filter(k => active(cons[k], w));
    const freed = actV.filter(k => !along.includes(k));
    let s = 0;
    for (const k of freed) s = Math.max(s, -dot(cons[k].a, d));
    if (s <= 0) s = Math.hypot(d[0], d[1]);
    return { to: j, rate: dot(c, d) / s, freed, along, far: isFar(w) };
  });
}

// ---------------------------------------------------------------- dense simplex
export type Sense = "<=" | ">=" | "=";
export interface LPResult { status: "optimal" | "infeasible" | "unbounded"; x: number[]; z: number }

const TOL = 1e-9;

/** max c.x subject to A x (sense) b, x >= 0. Two-phase tableau, Bland's rule. */
export function solveLP(c: number[], A: number[][], sense: Sense[], b: number[]): LPResult {
  const n = c.length, m = A.length;
  const rows = A.map(r => r.slice()), rhs = b.slice(), sn = sense.slice();
  for (let i = 0; i < m; i++) if (rhs[i] < 0) {
    rows[i] = rows[i].map(v => -v); rhs[i] = -rhs[i];
    sn[i] = sn[i] === "<=" ? ">=" : sn[i] === ">=" ? "<=" : "=";
  }
  // columns: x (n), one slack/surplus per inequality, one artificial per >= or = row
  const nS = sn.filter(s => s !== "=").length, nA = sn.filter(s => s !== "<=").length;
  const N = n + nS + nA;
  const T: number[][] = [], basis: number[] = [];
  const isArt = new Array(N).fill(false);
  let si = n, ai = n + nS;
  for (let i = 0; i < m; i++) {
    const row = new Array(N + 1).fill(0);
    for (let j = 0; j < n; j++) row[j] = rows[i][j];
    row[N] = rhs[i];
    if (sn[i] === "<=") { row[si] = 1; basis.push(si++); }
    else {
      if (sn[i] === ">=") row[si++] = -1;
      row[ai] = 1; isArt[ai] = true; basis.push(ai++);
    }
    T.push(row);
  }
  const pivot = (r: number, q: number) => {
    const p = T[r][q];
    for (let j = 0; j <= N; j++) T[r][j] /= p;
    for (let i = 0; i < m; i++) if (i !== r && T[i][q] !== 0) {
      const f = T[i][q];
      for (let j = 0; j <= N; j++) T[i][j] -= f * T[r][j];
    }
    basis[r] = q;
  };
  const run = (cost: number[], allowed: boolean[]): "optimal" | "unbounded" => {
    for (let it = 0; it < 5000; it++) {
      let q = -1;
      for (let j = 0; j < N && q < 0; j++) {
        if (!allowed[j] || basis.includes(j)) continue;
        let rc = cost[j];
        for (let i = 0; i < m; i++) rc -= cost[basis[i]] * T[i][j];
        if (rc > TOL) q = j;
      }
      if (q < 0) return "optimal";
      let r = -1, best = Infinity;
      for (let i = 0; i < m; i++) if (T[i][q] > TOL) {
        const ratio = T[i][N] / T[i][q];
        if (ratio < best - TOL || (Math.abs(ratio - best) <= TOL && basis[i] < basis[r])) { best = ratio; r = i; }
      }
      if (r < 0) return "unbounded";
      pivot(r, q);
    }
    return "optimal";
  };
  const all = new Array(N).fill(true);
  if (nA > 0) {
    run(isArt.map(a => a ? -1 : 0), all);
    let infeas = 0;
    for (let i = 0; i < m; i++) if (isArt[basis[i]]) infeas += T[i][N];
    if (infeas > 1e-7) return { status: "infeasible", x: [], z: NaN };
    for (let i = 0; i < m; i++) if (isArt[basis[i]]) {
      for (let j = 0; j < N; j++) if (!isArt[j] && Math.abs(T[i][j]) > TOL) { pivot(i, j); break; }
    }
  }
  const cost = new Array(N).fill(0);
  for (let j = 0; j < n; j++) cost[j] = c[j];
  const st = run(cost, isArt.map(a => !a));
  if (st === "unbounded") return { status: "unbounded", x: [], z: Infinity };
  const x = new Array(n).fill(0);
  for (let i = 0; i < m; i++) if (basis[i] < n) x[basis[i]] = T[i][N];
  let z = 0;
  for (let j = 0; j < n; j++) z += c[j] * x[j];
  return { status: "optimal", x, z };
}

// ---------------------------------------------------------------- toy FBA network
export const RXN = ["EX_glc", "EX_o2", "GLY", "BYP1", "BYP2", "RES", "FER", "EX_ac", "BIO"];
export const MET = ["G", "O2", "X", "P", "ATP", "Ac"];
export const S = [
  [1, 0, -1, -1, 0, 0, 0, 0, 0],
  [0, 1, 0, 0, 0, -2, 0, 0, 0],
  [0, 0, 0, 1, -1, 0, 0, 0, 0],
  [0, 0, 2, 0, 2, -1, -1, 0, -1],
  [0, 0, 2, 0, 2, 11, 1, 0, -5],
  [0, 0, 0, 0, 0, 0, 1, -1, 0],
];
const BIO = 8;

export interface FBAOpts { glc: number; o2: number; ko: boolean[] }

function fbaRows(o: FBAOpts, extra: { row: number[]; s: Sense; b: number }[] = []) {
  const A: number[][] = S.map(r => r.slice()), sense: Sense[] = S.map(() => "=" as Sense), b: number[] = S.map(() => 0);
  const unit = (j: number) => RXN.map((_, k) => k === j ? 1 : 0);
  A.push(unit(0)); sense.push("<="); b.push(o.glc);
  A.push(unit(1)); sense.push("<="); b.push(o.o2);
  o.ko.forEach((k, j) => { if (k) { A.push(unit(j)); sense.push("<="); b.push(0); } });
  for (const e of extra) { A.push(e.row); sense.push(e.s); b.push(e.b); }
  return { A, sense, b };
}

/** Maximise BIO. */
export function fba(o: FBAOpts): LPResult {
  const { A, sense, b } = fbaRows(o);
  return solveLP(RXN.map((_, j) => j === BIO ? 1 : 0), A, sense, b);
}

/** Flux variability: with BIO fixed at its maximum, min and max of every flux. */
export function fva(o: FBAOpts, zbest: number): [number, number][] {
  const fix = { row: RXN.map((_, j) => j === BIO ? 1 : 0), s: ">=" as Sense, b: zbest - 1e-9 * (1 + zbest) };
  const { A, sense, b } = fbaRows(o, [fix]);
  return RXN.map((_, j) => {
    const e = RXN.map((_, k) => k === j ? 1 : 0);
    const hi = solveLP(e, A, sense, b), lo = solveLP(e.map(v => -v), A, sense, b);
    return [lo.status === "optimal" ? -lo.z : NaN, hi.status === "optimal" ? hi.z : Infinity];
  });
}

export function maxResidual(v: number[]): number {
  let r = 0;
  for (const row of S) { let s = 0; row.forEach((a, j) => s += a * v[j]); r = Math.max(r, Math.abs(s)); }
  return r;
}

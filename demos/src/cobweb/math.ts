// Population maps and the two root/steady-state methods of the bisection-fixedpoint unit.
// K = 100 throughout (thousands of fish for the harvested stock).
//   bh:      Beverton–Holt  g(N) = R N / (1 + (R-1) N / K)
//   ricker:  Ricker         g(N) = N exp(r (1 - N/K))
//   harvest: Ricker minus a constant catch H: g(N) = N exp(r (1 - N/K)) - H
export type Model = "bh" | "ricker" | "harvest";
export interface Params { model: Model; r: number; H: number }
export const K = 100;

export function g(p: Params, N: number): number {
  if (p.model === "bh") return p.r * N / (1 + (p.r - 1) * N / K);
  const v = N * Math.exp(p.r * (1 - N / K));
  return p.model === "harvest" ? v - p.H : v;
}

export function dg(p: Params, N: number): number {
  if (p.model === "bh") { const q = 1 + (p.r - 1) * N / K; return p.r / (q * q); }
  return Math.exp(p.r * (1 - N / K)) * (1 - p.r * N / K);
}

export const f = (p: Params, N: number): number => g(p, N) - N;

/** Upper end of the plotting range: 150, or more if the Ricker hump is higher. */
export function range(p: Params): number {
  if (p.model === "bh") return 150;
  if (p.r <= 1) return 150;
  const top = K / p.r * Math.exp(p.r - 1);
  return Math.max(150, Math.ceil(top * 1.04 / 50) * 50);
}

/** Bisection to full precision on a bracket with a sign change. */
export function refine(p: Params, a: number, b: number): number {
  let fa = f(p, a);
  for (let i = 0; i < 200; i++) {
    const m = 0.5 * (a + b);
    if (m <= a || m >= b) break;
    const fm = f(p, m);
    if (fm === 0) return m;
    if ((fm < 0) === (fa < 0)) { a = m; fa = fm; } else b = m;
  }
  return 0.5 * (a + b);
}

export interface Steady { x: number; d: number }
/** Steady states in [0, hi]: N = 0 for the unharvested maps, and sign changes of f on a fine grid. */
export function steadyStates(p: Params, hi: number): Steady[] {
  const out: number[] = [];
  if (p.model !== "harvest") out.push(0);
  const n = 3000;
  let xa = p.model === "harvest" ? 0 : hi / n * 0.5, fa = f(p, xa);
  for (let i = 1; i <= n; i++) {
    const xb = hi * i / n, fb = f(p, xb);
    if (fa === 0 && xa > 0) out.push(xa);
    else if ((fa < 0) !== (fb < 0) && fb !== 0) out.push(refine(p, xa, xb));
    xa = xb; fa = fb;
  }
  return out.map(x => ({ x, d: dg(p, x) }));
}

/** Largest catch with a steady state: max over N of N (exp(r(1-N/K)) - 1). Returns [H, N]. */
export function maxCatch(r: number): [number, number] {
  const s = (N: number) => N * (Math.exp(r * (1 - N / K)) - 1);
  let a = 0, b = K;                     // golden-section search; s is unimodal on [0, K]
  const phi = (Math.sqrt(5) - 1) / 2;
  for (let i = 0; i < 120; i++) {
    const c = b - phi * (b - a), d = a + phi * (b - a);
    if (s(c) > s(d)) b = d; else a = c;
  }
  const N = 0.5 * (a + b);
  return [s(N), N];
}

export interface Orbit { x: number[]; collapse: number }   // collapse: year the stock hit 0, or -1
/** N_0 ... N_n. A harvested stock that would go negative is set to 0 and stays there. */
export function orbit(p: Params, N0: number, n: number): Orbit {
  const x = [N0];
  let collapse = -1;
  for (let t = 1; t <= n; t++) {
    const prev = x[t - 1];
    let v = collapse >= 0 ? 0 : g(p, prev);
    if (!(v > 0) && collapse < 0 && prev > 0) { collapse = t; v = 0; }
    if (!(v > 0)) v = 0;
    x.push(v);
  }
  return { x, collapse };
}

export type LongRun =
  | { kind: "collapse"; year: number }
  | { kind: "cycle"; values: number[] }      // period = values.length (1 = converges)
  | { kind: "none" };
/** Iterate 3000 years, then look at the next 64: how many distinct values (6 decimals)? */
export function longRun(p: Params, N0: number): LongRun {
  const o = orbit(p, N0, 3000 + 64);
  if (o.collapse >= 0) return { kind: "collapse", year: o.collapse };
  const tail = o.x.slice(3001);
  const seen = new Map<number, number>();
  for (const v of tail) { const k = Math.round(v * 1e6); if (!seen.has(k)) seen.set(k, v); }
  if (seen.size > 32) return { kind: "none" };
  // values in the order of the orbit, starting from the largest
  const per = seen.size, start = tail.indexOf(Math.max(...tail.slice(0, per)));
  return { kind: "cycle", values: tail.slice(start, start + per) };
}

export interface BisStep { a: number; b: number; m: number; fm: number }
/** Bisection steps from [a, b]; null if f(a), f(b) have the same sign. */
export function bisect(p: Params, a: number, b: number, n: number): BisStep[] | null {
  let fa = f(p, a);
  const fb = f(p, b);
  if (fa === 0 || fb === 0) return [];
  if ((fa < 0) === (fb < 0)) return null;
  const out: BisStep[] = [];
  for (let k = 0; k < n; k++) {
    const m = 0.5 * (a + b), fm = f(p, m);
    out.push({ a, b, m, fm });
    if (fm === 0) break;
    if ((fm < 0) === (fa < 0)) { a = m; fa = fm; } else b = m;
  }
  return out;
}

// Problems and one-dimensional minimisers for the golden demo (optimisation-1d unit).
// Every minimiser is a step-by-step state machine: one call of step() = one new
// evaluation of f (the first golden step evaluates two points, the first parabolic step
// three), so the demo can show every evaluation.

export const EPS = 2.220446049250313e-16;
export const G = (Math.sqrt(5) - 1) / 2;          // 0.618...
export const C = (3 - Math.sqrt(5)) / 2;          // 0.382... = 1 - G

// ---------------------------------------------------------------- the three problems
// mRNA decay (lecture): SSE(k) = sum (y - exp(-k t))^2
const TD = [0, 0.5, 1, 2, 4, 6, 8];
const YD = [1.000, 0.804, 0.684, 0.484, 0.217, 0.094, 0.083];
export const sse = (k: number): number => {
  let s = 0;
  for (let i = 0; i < TD.length; i++) { const r = YD[i] - Math.exp(-k * TD[i]); s += r * r; }
  return s;
};
const sse1 = (k: number): number => {      // SSE'
  let s = 0;
  for (let i = 0; i < TD.length; i++) { const e = Math.exp(-k * TD[i]); s += 2 * (YD[i] - e) * TD[i] * e; }
  return s;
};
const sse2 = (k: number): number => {      // SSE''
  let s = 0;
  for (let i = 0; i < TD.length; i++) { const e = Math.exp(-k * TD[i]); s += 2 * TD[i] * TD[i] * e * (2 * e - YD[i]); }
  return s;
};

// zero-truncated Poisson (Newton unit): minus log-likelihood, as on the lecture slide
export const nll = (l: number): number => -(173 * Math.log(l) - 100 * l - 100 * Math.log(-Math.expm1(-l)));
const score = (l: number): number => 173 / l - 100 / (-Math.expm1(-l));
const dscore = (l: number): number => -173 / (l * l) + 100 * Math.exp(-l) / (Math.expm1(-l) ** 2);

// circadian bioluminescence (lecture, set.seed(15)): RSS of the best sine wave with period T
const TC = Array.from({ length: 36 }, (_, i) => 4 * i);
const YC = [2.10, 3.07, 2.42, 2.29, 1.64, 1.08, 1.95, 2.81, 2.52, 1.77, 1.80, 1.32,
  1.92, 2.06, 3.02, 2.17, 1.59, 1.42, 1.44, 2.22, 3.01, 2.67, 1.54, 1.27,
  2.03, 2.44, 2.82, 2.11, 1.82, 1.67, 1.52, 2.24, 3.03, 2.60, 2.08, 1.55];
/** Residual sum of squares of y on (1, cos, sin), by modified Gram-Schmidt (a QR fit). */
export const rss = (T: number): number => {
  const n = TC.length;
  const cols = [TC.map(() => 1), TC.map(t => Math.cos(2 * Math.PI * t / T)), TC.map(t => Math.sin(2 * Math.PI * t / T))];
  const r = YC.slice();
  const q: number[][] = [];
  for (const c of cols) {
    const v = c.slice();
    for (const qq of q) { let d = 0; for (let i = 0; i < n; i++) d += qq[i] * v[i]; for (let i = 0; i < n; i++) v[i] -= d * qq[i]; }
    let nn = 0; for (let i = 0; i < n; i++) nn += v[i] * v[i];
    nn = Math.sqrt(nn);
    if (!(nn > 1e-12)) continue;                 // a degenerate column (T dividing 4): skip it
    for (let i = 0; i < n; i++) v[i] /= nn;
    q.push(v);
    let d = 0; for (let i = 0; i < n; i++) d += v[i] * r[i];
    for (let i = 0; i < n; i++) r[i] -= d * v[i];
  }
  let s = 0; for (let i = 0; i < n; i++) s += r[i] * r[i];
  return s;
};

export type ProbKey = "decay" | "circadian" | "ztp";
export interface Problem {
  key: ProbKey;
  label: string;          // button
  f: (x: number) => number;
  xname: string;          // axis label
  fname: string;
  view: [number, number]; // x range of the plot = default bracket
  ylim: [number, number];
  xstar: number;          // global minimiser (for the error)
  fstar: number;
  f2: number;             // f'' at xstar
  other?: number;         // a second local minimiser
  digits: number;         // decimals for x in tables
}

function newton(g: (x: number) => number, dg: (x: number) => number, x: number): number {
  for (let i = 0; i < 60; i++) {
    const s = g(x) / dg(x);
    x -= s;
    if (Math.abs(s) < 1e-16 * Math.abs(x)) break;
  }
  return x;
}
/** Minimiser of a smooth f by the secant method on a central-difference f' (only for RSS). */
function minByDeriv(f: (x: number) => number, x0: number, x1: number): number {
  const h = 1e-5;
  const d = (x: number) => (f(x + h) - f(x - h)) / (2 * h);
  let d0 = d(x0), d1 = d(x1);
  for (let i = 0; i < 60 && d1 !== d0; i++) {
    const x2 = x1 - d1 * (x1 - x0) / (d1 - d0);
    x0 = x1; d0 = d1; x1 = x2; d1 = d(x1);
    if (Math.abs(x1 - x0) < 1e-14 * Math.abs(x1)) break;
  }
  return x1;
}
const f2num = (f: (x: number) => number, x: number, h: number) => (f(x + h) - 2 * f(x) + f(x - h)) / (h * h);

const kStar = newton(sse1, sse2, 0.375);
const lStar = newton(score, dscore, 1.2);
const TStar = minByDeriv(rss, 24.6, 24.65);
const TLocal = minByDeriv(rss, 32.4, 32.41);

export const PROBLEMS: Record<ProbKey, Problem> = {
  decay: {
    key: "decay", label: "mRNA decay SSE(k)", f: sse, xname: "decay rate k (1/h)", fname: "SSE(k)",
    view: [0, 1], ylim: [-0.03, 0.62], xstar: kStar, fstar: sse(kStar), f2: sse2(kStar), digits: 4,
  },
  circadian: {
    key: "circadian", label: "circadian RSS(T)", f: rss, xname: "period T (h)", fname: "RSS(T)",
    view: [10, 40], ylim: [0, 12.5], xstar: TStar, fstar: rss(TStar), f2: f2num(rss, TStar, 1e-3),
    other: TLocal, digits: 3,
  },
  ztp: {
    key: "ztp", label: "ZTP −log L(λ)", f: nll, xname: "Poisson rate λ", fname: "−log L(λ)",
    view: [0.2, 4], ylim: [45, 175], xstar: lStar, fstar: nll(lStar), f2: -dscore(lStar), digits: 4,
  },
};

/** Half-width of the flat bottom: sqrt(2 eps |f*| / f''). */
export const flatWidth = (p: Problem): number => Math.sqrt(2 * EPS * Math.abs(p.fstar) / p.f2);

// ---------------------------------------------------------------- minimisers
export type Method = "golden" | "parabolic" | "brent";
export interface Pt { x: number; f: number }
export interface Row {
  n: number;              // evaluations so far
  kind: string;           // "start", "golden", "parabolic"
  x: number; f: number;   // the new point (the last one, if several)
  a: number; b: number;   // bracket after the step (NaN for parabolic)
  best: number;           // current estimate
  parab?: Pt[];           // the three points the parabola went through
  note?: string;
}
export interface Run {
  method: Method;
  rows: Row[];
  evals: Pt[];
  done: boolean;
  why: string;
  step(): void;
  // for drawing the current state
  bracket(): [number, number] | null;
  inner(): Pt[];          // golden: x1, x2; Brent: x, w, v; parabolic: the three points
  next(): string;         // what the next step will do
}

const MAXEV = 100;

/** Golden-section search, exactly as on the lecture slide. */
export function golden(f: (x: number) => number, a0: number, b0: number, tol: number): Run {
  let a = a0, b = b0, x1 = NaN, x2 = NaN, f1 = NaN, f2 = NaN;
  const run: Run = {
    method: "golden", rows: [], evals: [], done: false, why: "",
    step() {
      if (run.done) return;
      const wPrev = b - a;
      if (run.evals.length === 0) {
        x1 = b - G * (b - a); f1 = f(x1); x2 = a + G * (b - a); f2 = f(x2);
        run.evals.push({ x: x1, f: f1 }, { x: x2, f: f2 });
        run.rows.push({ n: 2, kind: "start", x: x2, f: f2, a, b, best: f1 < f2 ? x1 : x2, note: "x₁, x₂" });
      } else {
        let xn: number, fn: number, keep: string;
        if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - G * (b - a); f1 = f(x1); xn = x1; fn = f1; keep = "[a, x₂]"; }
        else { a = x1; x1 = x2; f1 = f2; x2 = a + G * (b - a); f2 = f(x2); xn = x2; fn = f2; keep = "[x₁, b]"; }
        run.evals.push({ x: xn, f: fn });
        run.rows.push({ n: run.evals.length, kind: "golden", x: xn, f: fn, a, b, best: f1 < f2 ? x1 : x2, note: keep });
      }
      if (!(b - a > tol)) { run.done = true; run.why = `width ${fmtSci(b - a)} ≤ tol`; }
      else if (run.rows.length > 1 && !(b - a < wPrev)) { run.done = true; run.why = "the bracket stopped shrinking: a and b are (nearly) neighbouring doubles"; }
      else if (run.evals.length >= MAXEV) { run.done = true; run.why = `${MAXEV} evaluations`; }
    },
    bracket: () => [a, b],
    inner: () => (run.evals.length ? [{ x: x1, f: f1 }, { x: x2, f: f2 }] : []),
    next: () => {
      if (run.done) return "";
      if (!run.evals.length) return `evaluate x₁ = a + 0.382(b − a) and x₂ = a + 0.618(b − a)`;
      return f1 < f2 ? `f(x₁) < f(x₂): keep [a, x₂] = [${fx(a)}, ${fx(x2)}]` : `f(x₁) ≥ f(x₂): keep [x₁, b] = [${fx(x1)}, ${fx(b)}]`;
    },
  };
  return run;
}

/** Vertex of the parabola through three points. */
export function vertex(p: Pt[]): number {
  const [A, B, Cc] = p;
  const num = (B.x - A.x) ** 2 * (B.f - Cc.f) - (B.x - Cc.x) ** 2 * (B.f - A.f);
  const den = (B.x - A.x) * (B.f - Cc.f) - (B.x - Cc.x) * (B.f - A.f);
  return B.x - 0.5 * num / den;
}
/** The parabola through three points, as a function. */
export function parabola(p: Pt[]): (x: number) => number {
  const [A, B, Cc] = p;
  return x => A.f * (x - B.x) * (x - Cc.x) / ((A.x - B.x) * (A.x - Cc.x)) +
    B.f * (x - A.x) * (x - Cc.x) / ((B.x - A.x) * (B.x - Cc.x)) +
    Cc.f * (x - A.x) * (x - B.x) / ((Cc.x - A.x) * (Cc.x - B.x));
}

/** Successive parabolic interpolation from a, (a+b)/2, b: the new vertex replaces the oldest point. */
export function parabolic(f: (x: number) => number, a0: number, b0: number, tol: number): Run {
  let pts: Pt[] = [];
  const run: Run = {
    method: "parabolic", rows: [], evals: [], done: false, why: "",
    step() {
      if (run.done) return;
      if (!pts.length) {
        pts = [a0, (a0 + b0) / 2, b0].map(x => ({ x, f: f(x) }));
        run.evals.push(...pts);
        run.rows.push({ n: 3, kind: "start", x: b0, f: pts[2].f, a: NaN, b: NaN, best: pts[1].x, note: "a, mid, b" });
        return;
      }
      const used = pts.slice();
      const xn = vertex(used);
      if (!Number.isFinite(xn)) { run.done = true; run.why = "the three points lie on a line: no vertex"; return; }
      const fn = f(xn);
      run.evals.push({ x: xn, f: fn });
      const step = Math.abs(xn - pts[2].x);
      pts = [pts[1], pts[2], { x: xn, f: fn }];
      const conc = parabola(used)(xn + 1) - parabola(used)(xn);   // > 0 if it opens upwards
      run.rows.push({ n: run.evals.length, kind: "parabolic", x: xn, f: fn, a: NaN, b: NaN, best: xn, parab: used,
        note: conc < 0 ? "opens down: a maximum!" : "" });
      if (!Number.isFinite(fn)) { run.done = true; run.why = `f(${fx(xn)}) is not a number: left the domain`; }
      else if (step <= tol + Math.sqrt(EPS) * Math.abs(xn)) { run.done = true; run.why = `step ${fmtSci(step)} ≤ tol + √ε·|x|`; }
      else if (Math.abs(xn) > 1e6) { run.done = true; run.why = "ran away"; }
      else if (run.evals.length >= MAXEV) { run.done = true; run.why = `${MAXEV} evaluations`; }
    },
    bracket: () => null,
    inner: () => pts,
    next: () => run.done ? "" : pts.length ? "jump to the vertex of the parabola through the last three points" : "evaluate a, (a + b)/2 and b",
  };
  return run;
}

/** Brent's minimiser, as in R's optimize (Brent_fmin in R's optimize.c). */
export function brent(f: (x: number) => number, ax: number, bx: number, tol: number): Run {
  const eps = Math.sqrt(EPS), tol3 = tol / 3;
  let a = ax, b = bx;
  let v = a + C * (b - a), w = v, x = v, d = 0, e = 0, fx_ = NaN, fv = NaN, fw = NaN;
  const run: Run = {
    method: "brent", rows: [], evals: [], done: false, why: "",
    step() {
      if (run.done) return;
      if (!run.evals.length) {
        fx_ = f(x); fv = fx_; fw = fx_;
        run.evals.push({ x, f: fx_ });
        run.rows.push({ n: 1, kind: "start", x, f: fx_, a, b, best: x, note: "a + 0.382(b − a)" });
        checkStop();
        return;
      }
      const xm = (a + b) / 2, tol1 = eps * Math.abs(x) + tol3, t2 = 2 * tol1;
      let p = 0, q = 0, r = 0, u: number;
      if (Math.abs(e) > tol1) {
        r = (x - w) * (fx_ - fv);
        q = (x - v) * (fx_ - fw);
        p = (x - v) * q - (x - w) * r;
        q = (q - r) * 2;
        if (q > 0) p = -p; else q = -q;
        r = e; e = d;
      }
      let kind: string;
      let parab: Pt[] | undefined;
      if (Math.abs(p) >= Math.abs(q * 0.5 * r) || p <= q * (a - x) || p >= q * (b - x)) {
        e = x < xm ? b - x : a - x;
        d = C * e;
        kind = "golden";
      } else {
        d = p / q; u = x + d;
        if (u - a < t2 || b - u < t2) { d = tol1; if (x >= xm) d = -d; }
        kind = "parabolic";
        parab = [{ x: v, f: fv }, { x: w, f: fw }, { x, f: fx_ }];
      }
      if (Math.abs(d) >= tol1) u = x + d;
      else if (d > 0) u = x + tol1;
      else u = x - tol1;
      const fu = f(u);
      run.evals.push({ x: u, f: fu });
      if (fu <= fx_) {
        if (u < x) b = x; else a = x;
        v = w; w = x; x = u;
        fv = fw; fw = fx_; fx_ = fu;
      } else {
        if (u < x) a = u; else b = u;
        if (fu <= fw || w === x) { v = w; fv = fw; w = u; fw = fu; }
        else if (fu <= fv || v === x || v === w) { v = u; fv = fu; }
      }
      run.rows.push({ n: run.evals.length, kind, x: u, f: fu, a, b, best: x, parab });
      checkStop();
    },
    bracket: () => [a, b],
    inner: () => run.evals.length ? [{ x, f: fx_ }, { x: w, f: fw }, { x: v, f: fv }] : [],
    next: () => run.done ? "" : run.evals.length ? "parabolic step if it behaves, otherwise golden" : "evaluate x = a + 0.382(b − a)",
  };
  function checkStop(): void {
    const xm = (a + b) / 2, tol1 = eps * Math.abs(x) + tol3, t2 = 2 * tol1;
    if (Math.abs(x - xm) <= t2 - (b - a) / 2) { run.done = true; run.why = "converged (R's stopping rule)"; }
    else if (run.evals.length >= MAXEV) { run.done = true; run.why = `${MAXEV} evaluations`; }
  }
  return run;
}

export function makeRun(m: Method, f: (x: number) => number, a: number, b: number, tol: number): Run {
  return m === "golden" ? golden(f, a, b, tol) : m === "parabolic" ? parabolic(f, a, b, tol) : brent(f, a, b, tol);
}

/** Run to the end. */
export function runAll(r: Run): Run { while (!r.done) r.step(); return r; }

// ---------------------------------------------------------------- formatting
export function fmtSci(v: number, d = 2): string {
  if (!Number.isFinite(v)) return String(v);
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e-3 && a < 1e4) return v.toPrecision(d);
  const [m, e] = v.toExponential(d - 1).split("e");
  return `${m}e${Number(e)}`;
}
export function fx(v: number, dec = 4): string {
  if (!Number.isFinite(v)) return String(v);
  if (Math.abs(v) >= 1e4) return v.toExponential(2).replace("e+", "e");
  return v.toFixed(dec);
}

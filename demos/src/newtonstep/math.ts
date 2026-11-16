// The three problems of the newton unit and one step of Newton, secant or bisection.
// Pure functions, no DOM, so they can be checked with node.

export type ProbKey = "dose" | "score" | "mm";
export type Method = "newton" | "secant" | "bisection";

export interface Problem {
  key: ProbKey;
  label: string;           // button text
  sym: string;             // variable name: c, λ, u
  fname: string;           // y-axis label
  xlab: string;            // x-axis label
  rootName: string;        // name of the wanted root
  f: (x: number) => number;
  df: (x: number) => number;
  /** Newton step written in closed form (Michaelis–Menten: u - sinh u, as in the lecture). */
  newtonStep?: (x: number) => number;
  root: number;
  positive: string;        // non-empty: x < 0 is meaningless, and this is the word for x
  xr: [number, number]; yr: [number, number];
  xticks: number[]; yticks: number[];
  res: number;             // dragging snaps to multiples of this
  defaults: Record<Method, number[]>;
}

// Two-site dose-response: R(c) = 0.3 c/(0.1 + c) + 0.7 c^2/(25 + c^2), f = R - 0.5.
const dose: Problem = {
  key: "dose", label: "dose–response f(c)", sym: "c", fname: "f(c) = R(c) − 0.5", xlab: "dose c (µM)",
  rootName: "c₅₀",
  f: c => 0.3 * c / (0.1 + c) + 0.7 * c * c / (25 + c * c) - 0.5,
  df: c => 0.03 / ((0.1 + c) * (0.1 + c)) + 35 * c / ((25 + c * c) * (25 + c * c)),
  root: 3.261305265055212,
  positive: "dose",
  xr: [-8, 12], yr: [-0.6, 0.6],
  xticks: [-8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12], yticks: [-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6],
  res: 0.05,
  defaults: { newton: [1], secant: [1, 2], bisection: [1, 10] },
};

// Zero-truncated Poisson score: 100 cells, sum of counts 173.
const N = 100, SX = 173;
const score: Problem = {
  key: "score", label: "score S(λ)", sym: "λ", fname: "score S(λ)", xlab: "Poisson rate λ",
  rootName: "λ̂",
  f: l => SX / l - N / (-Math.expm1(-l)),
  df: l => -SX / (l * l) + N * Math.exp(-l) / (Math.expm1(-l) ** 2),
  root: 1.2184695258301268,
  positive: "rate",
  xr: [-1.5, 4.5], yr: [-100, 150],
  xticks: [-1, 0, 1, 2, 3, 4], yticks: [-100, -50, 0, 50, 100, 150],
  res: 0.01,
  defaults: { newton: [1.73], secant: [2, 2.5], bisection: [0.5, 4] },
};

// Michaelis-Menten in log substrate u = ln(S/Km): sigma(u) - 1/2, root u = 0 (S = Km).
const sig = (u: number) => 1 / (1 + Math.exp(-u));
const mm: Problem = {
  key: "mm", label: "Michaelis–Menten σ(u) − ½", sym: "u", fname: "v/Vmax − 0.5", xlab: "u = ln(S / Km)",
  rootName: "u* = 0",
  f: u => sig(u) - 0.5,
  df: u => { const s = sig(u); return s * (1 - s); },
  newtonStep: u => u - Math.sinh(u),
  root: 0,
  positive: "",
  xr: [-4.5, 4.5], yr: [-0.6, 0.6],
  xticks: [-4, -3, -2, -1, 0, 1, 2, 3, 4], yticks: [-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6],
  res: 0.01,
  defaults: { newton: [1], secant: [1, 2], bisection: [-3, 4] },
};

export const PROBLEMS: Record<ProbKey, Problem> = { dose, score, mm };

export const TOL = 1e-12;   // as the lecture's newton() and secant(): relative step
export const MAXIT = 60;

export interface Row {
  k: number;               // 0 = start value(s); steps are 1, 2, ...
  x: number; fx: number;
  /** points the step was built from: Newton [x_{k-1}], secant [x_{k-2}, x_{k-1}] */
  from: number[];
  a?: number; b?: number;  // bisection: the bracket after this step
}

export interface Run {
  prob: Problem; method: Method;
  rows: Row[];             // start rows (k = 0) then one row per step
  evals: number;
  done: boolean;
  status: string;          // last message
  level: "" | "good" | "warn";
  cycling: boolean;
  fa?: number; fb?: number; a?: number; b?: number;
}

export function start(prob: Problem, method: Method, pts: number[]): Run {
  const run: Run = { prob, method, rows: [], evals: 0, done: false, status: "", level: "", cycling: false };
  if (method === "newton") {
    run.rows.push({ k: 0, x: pts[0], fx: prob.f(pts[0]), from: [] });
  } else {
    const [p, q] = method === "bisection" ? [Math.min(pts[0], pts[1]), Math.max(pts[0], pts[1])] : [pts[0], pts[1]];
    run.rows.push({ k: 0, x: p, fx: prob.f(p), from: [] });
    run.rows.push({ k: 0, x: q, fx: prob.f(q), from: [] });
    run.evals = 2;
    if (method === "bisection") {
      run.a = p; run.b = q; run.fa = run.rows[0].fx; run.fb = run.rows[1].fx;
      if (!(run.fa * run.fb < 0) && run.fa !== 0 && run.fb !== 0) {
        run.done = true; run.level = "warn";
        run.status = `f(a) and f(b) have the same sign: no bracket, bisection cannot start.`;
      }
    } else if (p === q) {
      run.done = true; run.level = "warn"; run.status = "The two start values are equal.";
    }
  }
  return run;
}

export const steps = (run: Run) => run.rows.filter(r => r.k > 0).length;

/** Take one step. Returns false if the run had already stopped. */
export function step(run: Run): boolean {
  if (run.done) return false;
  const P = run.prob, last = run.rows[run.rows.length - 1], k = steps(run) + 1;
  let x: number, row: Row;
  run.status = ""; run.level = "";
  if (run.method === "newton") {
    const d = P.df(last.x);
    run.evals += 2;
    if (d === 0 && !P.newtonStep) { run.done = true; run.level = "warn"; run.status = `f′(${P.sym}) = 0: the tangent is flat, the step is a division by zero.`; return true; }
    x = P.newtonStep ? P.newtonStep(last.x) : last.x - last.fx / d;
    row = { k, x, fx: P.f(x), from: [last.x] };
  } else if (run.method === "secant") {
    const prev = run.rows[run.rows.length - 2];
    if (last.fx === prev.fx) { run.done = true; run.level = "warn"; run.status = `f(${P.sym}ₖ) = f(${P.sym}ₖ₋₁): the chord is flat, division by zero.`; return true; }
    x = last.x - last.fx * (last.x - prev.x) / (last.fx - prev.fx);
    run.evals += 1;
    row = { k, x, fx: P.f(x), from: [prev.x, last.x] };
  } else {
    const a = run.a!, b = run.b!;
    x = a + (b - a) / 2;
    const fm = P.f(x);
    run.evals += 1;
    if (Math.sign(fm) === Math.sign(run.fa!)) { run.a = x; run.fa = fm; } else { run.b = x; run.fb = fm; }
    row = { k, x, fx: fm, from: [a, b], a: run.a, b: run.b };
  }
  run.rows.push(row);
  const xp = run.method === "bisection" ? NaN : last.x;
  if (!Number.isFinite(x) || !Number.isFinite(row.fx)) {
    run.done = true; run.level = "warn";
    run.status = `${P.sym}${sub(k)} = ${x}${Number.isFinite(x) ? `, f = ${row.fx}` : ""}: the iteration has diverged.`;
    return true;
  }
  const conv = run.method === "bisection"
    ? (run.b! - run.a!) <= TOL * Math.max(Math.abs(x), 1) || row.fx === 0
    : Math.abs(x - xp) <= TOL * Math.max(Math.abs(x), 1);
  if (conv) {
    run.done = true;
    const wrong = Math.abs(x - P.root) > 1e-6 * Math.max(Math.abs(P.root), 1);
    if (wrong) {
      run.level = "warn";
      run.status = `Converged in ${k} steps to ${P.sym} = ${x.toPrecision(6)}: a root of f, but not ${P.rootName}` +
        (P.positive && x < 0 ? ` — a negative ${P.positive}!` : ".");
    } else {
      run.level = "good";
      run.status = `Converged: ${run.method === "bisection" ? "bracket" : "step"} below 10⁻¹²·|${P.sym}| after ${k} steps, ${run.evals} evaluations.`;
    }
    return true;
  }
  if (k >= MAXIT) { run.done = true; run.level = "warn"; run.status = `No convergence after ${MAXIT} steps (maxit).`; return true; }
  // warnings that do not stop the iteration
  const n = run.rows.length;
  if (run.method !== "bisection" && k >= 3) {
    const x2 = run.rows[n - 3].x;
    const cyc = Math.abs(x - x2) <= 1e-6 * Math.max(Math.abs(x), 1) && Math.abs(x - xp) > 1e-3;
    if (cyc) { run.level = "warn"; run.status = `Cycling: ${P.sym}${sub(k)} ≈ ${P.sym}${sub(k - 2)}, the iterates bounce back and forth.`; }
    run.cycling = cyc;
  }
  if (!run.status && P.positive && x < 0) { run.level = "warn"; run.status = `${P.sym}${sub(k)} = ${x.toPrecision(4)} < 0: a negative ${P.positive} has no meaning.`; }
  return true;
}

/** Error against the wanted root, and correct digits -log10(|e| / max(|root|, 1)). */
export function errDigits(prob: Problem, x: number): { e: number; d: number } {
  const e = Math.abs(x - prob.root);
  return { e, d: e === 0 ? Infinity : -Math.log10(e / Math.max(Math.abs(prob.root), 1)) };
}

const SUB = "₀₁₂₃₄₅₆₇₈₉";
export const sub = (k: number): string => String(k).split("").map(c => SUB[Number(c)]).join("");

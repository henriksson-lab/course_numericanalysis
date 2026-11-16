// Numerics for the rkstep demo (runge-kutta unit): explicit Runge-Kutta tableaux,
// one step with all stages recorded, fixed-step global errors, and an adaptive
// embedded-pair solver with the same controller as the unit's adapt_heun / dp5.

export type Vec = number[];
export type F = (t: number, y: Vec) => Vec;

// ---------------------------------------------------------------- fixed-step methods

export interface Tableau { name: string; c: number[]; A: number[][]; b: number[]; order: number }

export const METHODS: Record<string, Tableau> = {
  euler: { name: "Euler", c: [0], A: [[]], b: [1], order: 1 },
  midpoint: { name: "Midpoint", c: [0, 0.5], A: [[], [0.5]], b: [0, 1], order: 2 },
  heun: { name: "Heun", c: [0, 1], A: [[], [1]], b: [0.5, 0.5], order: 2 },
  rk4: { name: "RK4", c: [0, 0.5, 0.5, 1], A: [[], [0.5], [0, 0.5], [0, 0, 1]], b: [1 / 6, 1 / 3, 1 / 3, 1 / 6], order: 4 },
};

export interface Stage { t: number; y: number; k: number; via: number }
/** One step of a scalar ODE, recording every stage (point, slope, and the slope used to reach it). */
export function stepStages(tb: Tableau, f: (t: number, y: number) => number, t0: number, y0: number, h: number) {
  const st: Stage[] = [];
  const ks: number[] = [];
  for (let i = 0; i < tb.c.length; i++) {
    let incr = 0;
    for (let j = 0; j < i; j++) incr += (tb.A[i][j] ?? 0) * ks[j];
    const y = y0 + h * incr;
    const k = f(t0 + tb.c[i] * h, y);
    ks.push(k);
    st.push({ t: t0 + tb.c[i] * h, y, k, via: tb.c[i] > 0 ? incr / tb.c[i] : k });
  }
  let slope = 0;
  for (let i = 0; i < ks.length; i++) slope += tb.b[i] * ks[i];
  return { stages: st, slope, y1: y0 + h * slope };
}

/** n fixed steps of a scalar ODE from (t0, y0) with step h; returns y(t0 + n h). */
export function fixedSolve(tb: Tableau, f: (t: number, y: number) => number, t0: number, y0: number, h: number, n: number): number {
  let y = y0;
  for (let i = 0; i < n; i++) y = stepStages(tb, f, t0 + i * h, y, h).y1;
  return y;
}

// ---------------------------------------------------------------- scalar running examples

export interface Scalar {
  key: string; label: string; f: (t: number, y: number) => number; exact: (t: number) => number;
  y0: number; T: number; yname: string; unit: string; tunit: string;
}
const K = 1.2, N0 = 0.05, R = 0.8;
export const LOGISTIC: Scalar = {
  key: "logistic", label: "logistic OD", f: (_t, N) => R * N * (1 - N / K),
  exact: t => K / (1 + (K / N0 - 1) * Math.exp(-R * t)), y0: N0, T: 10, yname: "N", unit: "", tunit: "h",
};
export const CLEARANCE: Scalar = {
  key: "clearance", label: "drug clearance", f: (_t, C) => -0.8 * C,
  exact: t => 10 * Math.exp(-0.8 * t), y0: 10, T: 10, yname: "C", unit: " mg/L", tunit: "h",
};

// ---------------------------------------------------------------- systems

export const LV = { a: 1, b: 0.1, c: 0.075, d: 1.5 };
export const lvF: F = (_t, y) => [LV.a * y[0] - LV.b * y[0] * y[1], LV.c * y[0] * y[1] - LV.d * y[1]];
export const lvV = (x: number, y: number): number => LV.c * x - LV.d * Math.log(x) + LV.b * y - LV.a * Math.log(y);

const RA = 216, RA0 = 0.216, RB = 5, RN = 2;
export const repF: F = (_t, y) => {
  const m = y.slice(0, 3), p = y.slice(3, 6);
  const pr = [p[2], p[0], p[1]];
  return [
    ...m.map((mi, i) => -mi + RA / (1 + pr[i] ** RN) + RA0),
    ...p.map((pi, i) => -RB * (pi - m[i])),
  ];
};
export const REP_P1_50 = 47.335989;

// ---------------------------------------------------------------- embedded pairs

export interface Pair {
  key: string; name: string; c: number[]; A: number[][]; bHi: number[]; bLo: number[];
  q: number;       // order of the lower member: exponent -1/(q+1)
  fsal: boolean;   // last stage = f at the new point
}
export const PAIRS: Record<string, Pair> = {
  he: {
    key: "he", name: "Heun–Euler 2(1)", c: [0, 1], A: [[], [1]], bHi: [0.5, 0.5], bLo: [1, 0], q: 1, fsal: false,
  },
  bs: {
    key: "bs", name: "Bogacki–Shampine 3(2)", c: [0, 0.5, 0.75, 1],
    A: [[], [0.5], [0, 0.75], [2 / 9, 1 / 3, 4 / 9]],
    bHi: [2 / 9, 1 / 3, 4 / 9, 0], bLo: [7 / 24, 1 / 4, 1 / 3, 1 / 8], q: 2, fsal: true,
  },
  dp: {
    key: "dp", name: "Dormand–Prince 5(4)", c: [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1],
    A: [[], [1 / 5], [3 / 40, 9 / 40], [44 / 45, -56 / 15, 32 / 9],
      [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
      [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
      [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84]],
    bHi: [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0],
    bLo: [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40], q: 4, fsal: true,
  },
};

export interface AdaptResult {
  t: number[]; y: Vec[]; h: number[];       // accepted steps (h[i] = step that ended at t[i+1])
  rejT: number[]; rejH: number[];           // rejected attempts: start time and attempted h
  nfev: number; stopped: boolean;
}

/**
 * Adaptive embedded-pair solver. Keeps the higher-order answer; scaled RMS error with
 * sc = atol + rtol max(|y|, |ynew|); accept if err <= 1; h <- h min(5, max(0.2, 0.9 err^(-1/(q+1)))).
 */
export function adapt(p: Pair, f: F, t0: number, y0: Vec, tend: number, rtol: number, atol: number,
  h0: number, maxSteps = 400000): AdaptResult {
  const n = y0.length, s = p.c.length;
  let t = t0, y = y0.slice(), h = h0;
  const out: AdaptResult = { t: [t0], y: [y0.slice()], h: [], rejT: [], rejH: [], nfev: 0, stopped: false };
  let k1: Vec | null = null;
  while (t < tend) {
    if (out.h.length + out.rejT.length >= maxSteps) { out.stopped = true; break; }
    h = Math.min(h, tend - t);
    const ks: Vec[] = [];
    for (let i = 0; i < s; i++) {
      if (i === 0 && k1) { ks.push(k1); continue; }
      const ys = y.slice();
      for (let j = 0; j < i; j++) {
        const a = p.A[i][j] ?? 0;
        if (a !== 0) for (let m = 0; m < n; m++) ys[m] += h * a * ks[j][m];
      }
      ks.push(f(t + p.c[i] * h, ys));
      out.nfev++;
    }
    const yhi = y.slice(), ylo = y.slice();
    for (let i = 0; i < s; i++) for (let m = 0; m < n; m++) {
      yhi[m] += h * p.bHi[i] * ks[i][m];
      ylo[m] += h * p.bLo[i] * ks[i][m];
    }
    let ss = 0;
    for (let m = 0; m < n; m++) {
      const sc = atol + rtol * Math.max(Math.abs(y[m]), Math.abs(yhi[m]));
      ss += ((yhi[m] - ylo[m]) / sc) ** 2;
    }
    const err = Math.sqrt(ss / n);
    // on rejection an FSAL pair reuses k1; Heun-Euler recomputes it, as adapt_heun in the unit does
    k1 = p.fsal ? ks[0] : null;
    if (err <= 1 && Number.isFinite(err)) {
      t = t + h; y = yhi;
      out.t.push(t); out.y.push(y.slice()); out.h.push(h);
      k1 = p.fsal ? ks[s - 1] : null;
    } else { out.rejT.push(t); out.rejH.push(h); }
    const fac = Number.isFinite(err) ? (err === 0 ? 5 : Math.min(5, Math.max(0.2, 0.9 * err ** (-1 / (p.q + 1))))) : 0.2;
    h *= fac;
    if (h < 1e-14 * Math.max(1, Math.abs(t))) { out.stopped = true; break; }
  }
  return out;
}

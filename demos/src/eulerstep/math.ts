// Models and Euler runs for the eulerstep demo. Same numbers as euler/lec.qmd.

/** Euler for a scalar autonomous ODE y' = f(y): values y_0 ... y_n. */
export function euler1(f: (y: number) => number, y0: number, h: number, n: number): number[] {
  const y = new Array<number>(n + 1);
  y[0] = y0;
  for (let i = 0; i < n; i++) y[i + 1] = y[i] + h * f(y[i]);
  return y;
}

/** Euler for an autonomous system y' = f(y): one array per time point. */
export function eulerSys(f: (y: number[]) => number[], y0: number[], h: number, n: number): number[][] {
  const Y: number[][] = [y0.slice()];
  for (let i = 0; i < n; i++) {
    const y = Y[i], d = f(y);
    Y.push(y.map((v, j) => v + h * d[j]));
  }
  return Y;
}

/** Classical RK4 for a system, used only for the smooth reference curves. */
export function rk4Sys(f: (y: number[]) => number[], y0: number[], h: number, n: number): number[][] {
  const Y: number[][] = [y0.slice()];
  const add = (a: number[], b: number[], s: number) => a.map((v, j) => v + s * b[j]);
  for (let i = 0; i < n; i++) {
    const y = Y[i];
    const k1 = f(y), k2 = f(add(y, k1, h / 2)), k3 = f(add(y, k2, h / 2)), k4 = f(add(y, k3, h));
    Y.push(y.map((v, j) => v + h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])));
  }
  return Y;
}

/** Number of steps of size h that land exactly on T, or 0 if T is not a multiple of h. */
export function stepsTo(T: number, h: number): number {
  const n = Math.round(T / h);
  return n >= 1 && Math.abs(n * h - T) < 1e-9 * T ? n : 0;
}

// ---------------------------------------------------------------- one equation
export interface Scalar {
  key: "logistic" | "clearance";
  label: string;
  f: (y: number) => number;
  /** exact solution through (t0, y0), evaluated at t >= t0; NaN where undefined */
  exactFrom: (t: number, t0: number, y0: number) => number;
  y0: number;
  Te: number;            // time at which the global error is measured
  spans: number[];       // time windows to choose from
  ylim: [number, number];
  rate: number;          // k (clearance) or r (logistic, near K): stable for h < 2/rate
  rateName: string;
  hList: number[];
  yName: string;         // axis label
  unit: string;          // time unit
}

const K = 1.2, N0 = 0.05, R = 0.8, KC = 0.8, C0 = 10;

export const LOGISTIC: Scalar = {
  key: "logistic",
  label: "logistic growth",
  f: N => R * N * (1 - N / K),
  exactFrom: (t, t0, y0) => (y0 > 0 ? K / (1 + (K / y0 - 1) * Math.exp(-R * (t - t0))) : NaN),
  y0: N0, Te: 5, spans: [10, 60], ylim: [0, 1.7], rate: R, rateName: "r",
  hList: [0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 1.25, 2, 2.4, 2.5, 2.6, 3],
  yName: "OD  N(t)", unit: "h",
};

export const CLEARANCE: Scalar = {
  key: "clearance",
  label: "drug clearance",
  f: C => -KC * C,
  exactFrom: (t, t0, y0) => y0 * Math.exp(-KC * (t - t0)),
  y0: C0, Te: 6, spans: [12, 24], ylim: [-15, 15], rate: KC, rateName: "k",
  hList: [0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 1.25, 1.5, 2, 2.4, 2.5, 2.6, 3],
  yName: "C(t)  (mg/L)", unit: "h",
};

/** Global error at Te for a run with step h (NaN if Te is not a step point). */
export function globalErr(m: Scalar, h: number): number {
  const n = stepsTo(m.Te, h);
  if (!n) return NaN;
  const y = euler1(m.f, m.y0, h, n);
  return y[n] - m.exactFrom(m.Te, 0, m.y0);
}

/** Error-vs-h curve: h = Te/n for log-spaced integers n, down to hMin. */
export function errCurve(m: Scalar, hMin: number): { h: number; e: number }[] {
  const ns = new Set<number>();
  const nMax = Math.floor(m.Te / hMin);
  for (let l = 0; l <= Math.log10(nMax) + 1e-9; l += 0.05) ns.add(Math.max(1, Math.round(10 ** l)));
  return [...ns].sort((a, b) => a - b).map(n => {
    const h = m.Te / n;
    const y = euler1(m.f, m.y0, h, n);
    return { h, e: y[n] - m.exactFrom(m.Te, 0, m.y0) };
  });
}

// ---------------------------------------------------------------- systems
export const SIR = { beta: 0.3, gamma: 0.1, N: 10000, y0: [9990, 10, 0], T: 160 };
export const sirF = (y: number[]): number[] => {
  const inf = SIR.beta * y[0] * y[1] / SIR.N, rec = SIR.gamma * y[1];
  return [-inf, inf - rec, rec];
};
const R0 = SIR.beta / SIR.gamma;
/** Exact peak of I, from dI/dS: 3007.960706 */
export const IMAX = SIR.y0[0] + SIR.y0[1] - (SIR.N / R0) * (1 + Math.log(R0 * SIR.y0[0] / SIR.N));

/** Peak of I over the Euler grid points, and its day. */
export function sirPeak(Y: number[][], h: number): { peak: number; day: number } {
  let best = -Infinity, at = 0;
  Y.forEach((y, i) => { if (y[1] > best) { best = y[1]; at = i; } });
  return { peak: best, day: at * h };
}

export const OSC = { w: 2 * Math.PI / 24, y0: [1, 0], T: 240 };
export const oscF = (y: number[]): number[] => [y[1], -OSC.w * OSC.w * y[0]];
export const amp = (y: number[]): number => Math.hypot(y[0], y[1] / OSC.w);

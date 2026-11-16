// Models and numerics for the ssa demo (no DOM): the lecture's constitutive gene
// (birth-death, k = 10, delta = 1 per h) and telegraph gene (k_on, k_off = 9 k_on,
// k_m = 100, delta = 1 per h, so the promoter is on 10 % of the time and the mean mRNA
// is 10 for every k_on). Gillespie steps, tau-leaps, populations of cells, and the exact
// mean, variance and distribution at time t (moment equations, master equation).

export const K = 10, DELTA = 1, KM = 100, OFF_PER_ON = 9;

export type Kind = "const" | "tele";
export interface Model { kind: Kind; kon: number; koff: number; }
export const makeModel = (kind: Kind, kon: number): Model => ({ kind, kon, koff: OFF_PER_ON * kon });

export interface State { g: number; m: number; }
export const NAMES: Record<Kind, string[]> = {
  const: ["make", "degrade"],
  tele: ["on", "off", "make", "degrade"],
};

export function propensities(md: Model, s: State): number[] {
  return md.kind === "const"
    ? [K, DELTA * s.m]
    : [md.kon * (1 - s.g), md.koff * s.g, KM * s.g, DELTA * s.m];
}

export function apply(md: Model, s: State, j: number): State {
  if (md.kind === "const") return { g: 1, m: s.m + (j === 0 ? 1 : -1) };
  if (j === 0) return { g: 1, m: s.m };
  if (j === 1) return { g: 0, m: s.m };
  return { g: s.g, m: s.m + (j === 2 ? 1 : -1) };
}

export interface SsaStep {
  a: number[]; a0: number; u1: number; u2: number; tau: number; j: number;
  t0: number; before: State; after: State;
}

/** One step of Gillespie's direct method. u1 = 1 - unif() lies in (0, 1], so log(u1) is finite. */
export function ssaStep(md: Model, s: State, t0: number, unif: () => number): SsaStep {
  const a = propensities(md, s);
  const a0 = a.reduce((x, y) => x + y, 0);
  const u1 = 1 - unif(), u2 = unif();
  const tau = -Math.log(u1) / a0;
  const r = u2 * a0;
  let j = 0, c = a[0];
  while (j < a.length - 1 && (r >= c || a[j] === 0)) { j++; c += a[j]; }
  return { a, a0, u1, u2, tau, j, t0, before: s, after: apply(md, s, j) };
}

/** Poisson random number: multiplication method, splitting large means (exact). */
export function poisson(lam: number, unif: () => number): number {
  if (lam <= 0) return 0;
  if (lam > 40) { const h = lam / 2; return poisson(h, unif) + poisson(lam - h, unif); }
  const L = Math.exp(-lam);
  let k = 0, p = unif();
  while (p > L) { k++; p *= unif(); }
  return k;
}

export interface Leap { tau: number; t0: number; m0: number; lamMake: number; lamDeg: number; nMake: number; nDeg: number; raw: number; m: number; }

/** One tau-leap of the birth-death model; a negative count is recorded and clipped to 0. */
export function tauLeap(m0: number, t0: number, tau: number, unif: () => number): Leap {
  const lamMake = K * tau, lamDeg = DELTA * m0 * tau;
  const nMake = poisson(lamMake, unif), nDeg = poisson(lamDeg, unif);
  const raw = m0 + nMake - nDeg;
  return { tau, t0, m0, lamMake, lamDeg, nMake, nDeg, raw, m: Math.max(0, raw) };
}

export interface Population { counts: Int32Array; negFrac: number; events: number; }

/** N cells from g = 0, m = 0 to time T: SSA, or tau-leaping with step tau (constitutive only). */
export function population(md: Model, N: number, T: number, unif: () => number, tau = 0): Population {
  const counts = new Int32Array(N);
  let neg = 0, events = 0;
  for (let c = 0; c < N; c++) {
    let m = 0;
    if (tau > 0 && md.kind === "const") {
      let wentNeg = false;
      const n = Math.round(T / tau);
      for (let i = 0; i < n; i++) {
        m = m + poisson(K * tau, unif) - poisson(DELTA * m * tau, unif);
        if (m < 0) { wentNeg = true; m = 0; }
      }
      if (wentNeg) neg++;
      events += n;
    } else if (md.kind === "const") {
      let t = 0;
      for (;;) {
        const a1 = DELTA * m, a0 = K + a1;
        t -= Math.log(1 - unif()) / a0;
        if (t > T) break;
        if (unif() * a0 < K) m++; else m--;
        events++;
      }
    } else {
      let t = 0, g = 0;
      const kon = md.kon, koff = md.koff;
      for (;;) {
        const s1 = g ? koff : kon, s2 = g ? KM : 0, a3 = DELTA * m, a0 = s1 + s2 + a3;
        t -= Math.log(1 - unif()) / a0;
        if (t > T) break;
        const r = unif() * a0;
        if (r < s1) g = 1 - g; else if (r < s1 + s2) m++; else m--;
        events++;
      }
    }
    counts[c] = m;
  }
  return { counts, negFrac: neg / N, events };
}

export interface Stats { mean: number; variance: number; fano: number; se: number; zero: number; max: number; }
export function stats(x: Int32Array): Stats {
  const n = x.length;
  let s = 0, z = 0, mx = 0;
  for (let i = 0; i < n; i++) { s += x[i]; if (x[i] === 0) z++; if (x[i] > mx) mx = x[i]; }
  const mean = s / n;
  let ss = 0;
  for (let i = 0; i < n; i++) ss += (x[i] - mean) ** 2;
  const variance = n > 1 ? ss / (n - 1) : NaN;
  return { mean, variance, fano: variance / mean, se: Math.sqrt(variance / n), zero: z / n, max: mx };
}

/** Exact mean and variance of the mRNA at time T from g = m = 0 (closed moment equations, RK4). */
export function moments(md: Model, T: number): { mean: number; variance: number } {
  if (md.kind === "const") {
    const mean = (K / DELTA) * (1 - Math.exp(-DELTA * T));
    return { mean, variance: mean };   // Poisson at all times when m(0) = 0
  }
  const { kon, koff } = md;
  // y = [E g, E m, E gm, E m^2]
  const f = (y: number[]): number[] => [
    kon * (1 - y[0]) - koff * y[0],
    KM * y[0] - DELTA * y[1],
    kon * (y[1] - y[2]) - koff * y[2] + KM * y[0] - DELTA * y[2],
    KM * (2 * y[2] + y[0]) + DELTA * (y[1] - 2 * y[3]),
  ];
  const n = Math.max(50, Math.ceil(T * (kon + koff + DELTA) / 0.5));
  const h = T / n;
  let y = [0, 0, 0, 0];
  for (let i = 0; i < n; i++) y = rk4(f, y, h);
  return { mean: y[1], variance: y[3] - y[1] * y[1] };
}

function rk4(f: (y: number[]) => number[], y: number[], h: number): number[] {
  const k1 = f(y);
  const k2 = f(y.map((v, i) => v + h / 2 * k1[i]));
  const k3 = f(y.map((v, i) => v + h / 2 * k2[i]));
  const k4 = f(y.map((v, i) => v + h * k3[i]));
  return y.map((v, i) => v + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
}

/** Steady-state Fano factor of the telegraph mRNA. */
export const fanoSS = (kon: number): number => {
  const koff = OFF_PER_ON * kon;
  return 1 + KM * koff / ((kon + koff) * (kon + koff + DELTA));
};

/** Tau-leap mean and variance at T (unclipped recursion), leap tau. */
export function tauMoments(T: number, tau: number): { mean: number; variance: number } {
  let mu = 0, v = 0;
  const n = Math.round(T / tau), d = DELTA * tau;
  for (let i = 0; i < n; i++) { v = (1 - d) ** 2 * v + K * tau + d * mu; mu = mu + K * tau - d * mu; }
  return { mean: mu, variance: v };
}

/** Exact distribution of the mRNA count at T (Poisson, or the telegraph master equation by RK4). */
export function pmf(md: Model, T: number): Float64Array {
  if (md.kind === "const") {
    const lam = moments(md, T).mean, nmax = Math.max(40, Math.ceil(lam + 12 * Math.sqrt(lam) + 5));
    const p = new Float64Array(nmax + 1);
    p[0] = Math.exp(-lam);
    for (let n = 1; n <= nmax; n++) p[n] = p[n - 1] * lam / n;
    return p;
  }
  const { kon, koff } = md;
  const sd = Math.sqrt(10 * fanoSS(kon));
  const nmax = Math.min(600, Math.max(60, Math.ceil(10 + 20 * sd)));
  const S = nmax + 1;
  // P[0..S-1]: promoter off; P[S..2S-1]: on
  const f = (P: Float64Array, out: Float64Array): void => {
    for (let n = 0; n <= nmax; n++) {
      const p0 = P[n], p1 = P[S + n];
      let d0 = -(kon + DELTA * n) * p0 + koff * p1;
      let d1 = -(koff + (n < nmax ? KM : 0) + DELTA * n) * p1 + kon * p0;
      if (n < nmax) { d0 += DELTA * (n + 1) * P[n + 1]; d1 += DELTA * (n + 1) * P[S + n + 1]; }
      if (n > 0) d1 += KM * P[S + n - 1];
      out[n] = d0; out[S + n] = d1;
    }
  };
  const maxExit = koff + KM + DELTA * nmax + kon;
  const steps = Math.max(100, Math.ceil(T * maxExit / 1.2));
  const h = T / steps;
  let P = new Float64Array(2 * S), Q = new Float64Array(2 * S);
  P[0] = 1;
  const k1 = new Float64Array(2 * S), k2 = new Float64Array(2 * S), k3 = new Float64Array(2 * S), k4 = new Float64Array(2 * S);
  const tmp = new Float64Array(2 * S);
  for (let s = 0; s < steps; s++) {
    f(P, k1);
    for (let i = 0; i < 2 * S; i++) tmp[i] = P[i] + h / 2 * k1[i];
    f(tmp, k2);
    for (let i = 0; i < 2 * S; i++) tmp[i] = P[i] + h / 2 * k2[i];
    f(tmp, k3);
    for (let i = 0; i < 2 * S; i++) tmp[i] = P[i] + h * k3[i];
    f(tmp, k4);
    for (let i = 0; i < 2 * S; i++) Q[i] = P[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    [P, Q] = [Q, P];
  }
  const p = new Float64Array(S);
  for (let n = 0; n <= nmax; n++) p[n] = Math.max(0, P[n] + P[S + n]);
  return p;
}

/** ODE for the means (the lecture's ODE) on a grid 0..tEnd. */
export function odeMean(md: Model, tEnd: number, n = 400): { t: number[]; m: number[] } {
  const t: number[] = [], m: number[] = [];
  if (md.kind === "const") {
    for (let i = 0; i <= n; i++) { const ti = tEnd * i / n; t.push(ti); m.push(K / DELTA * (1 - Math.exp(-DELTA * ti))); }
    return { t, m };
  }
  const { kon, koff } = md;
  const f = (y: number[]): number[] => [kon * (1 - y[0]) - koff * y[0], KM * y[0] - DELTA * y[1]];
  const sub = Math.max(1, Math.ceil((tEnd / n) * (kon + koff + DELTA) / 0.5));
  const h = tEnd / n / sub;
  let y = [0, 0];
  t.push(0); m.push(0);
  for (let i = 1; i <= n; i++) {
    for (let s = 0; s < sub; s++) y = rk4(f, y, h);
    t.push(tEnd * i / n); m.push(y[1]);
  }
  return { t, m };
}

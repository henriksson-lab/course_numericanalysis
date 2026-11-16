// Numbers behind the Monte Carlo lab. Kept free of DOM code, so the values quoted on the
// exercise page can be checked in node.
import { rng } from "../lib/ui";

export type Rand = ReturnType<typeof rng>;
export type ProbKey = "auc" | "motif" | "rare";

// ---- the OD curve of the running example (K = 1.2, N0 = 0.05, r = 0.8), area on [0, 10] h ----
const K = 1.2, N0 = 0.05, R = 0.8;
export const od = (t: number): number => K / (1 + (K / N0 - 1) * Math.exp(-R * t));
export const AUC_EXACT = (K / R) * Math.log((Math.exp(10 * R) + K / N0 - 1) / (K / N0));   // 7.244448

/** Simpson's rule, for the true sd of 10 od(U). */
function simpson(f: (x: number) => number, a: number, b: number, n = 2000): number {
  const h = (b - a) / n;
  let s = f(a) + f(b);
  for (let i = 1; i < n; i++) s += (i % 2 ? 4 : 2) * f(a + i * h);
  return s * h / 3;
}
// Var(10 od(U)) = 10 * int_0^10 od^2 dt - AUC^2
const AUC_SD = Math.sqrt(10 * simpson(t => od(t) ** 2, 0, 10) - AUC_EXACT ** 2);

// ---- motif hits in promoters: G = 1e6, L = 10, 50 promoters of 1 kb every 20 kb from 10001 ----
export const G = 1e6, L = 10, NHITS = 200, KOBS = 18;
export const P_HIT = 50 * (1000 + L - 1) / (G - L + 1);          // 0.05045
export const MOTIF_EXACT = 0.01314003;                           // pbinom(17, 200, P_HIT, lower = FALSE)
/** Number of the 200 uniformly placed hits that overlap a promoter. */
export function overlaps(r: Rand): number {
  let k = 0;
  for (let h = 0; h < NHITS; h++) {
    const s = 1 + Math.floor(r.unif() * (G - L + 1));           // start position 1 .. G-L+1
    const t = s - (10001 - (L - 1));                            // 0 at the first start that overlaps promoter 1
    if (t < 0) continue;
    const j = Math.floor(t / 20000);
    if (j < 50 && t - 20000 * j < 1000 + L - 1) k++;
  }
  return k;
}

// ---- rare event P(X > 5), X ~ N(0, 1); importance sampling with g = N(5, 1) ----
export const RARE_EXACT = 2.866515718791933e-7;
const Q10 = 7.61985302416047e-24;                                 // P(X > 10)
// E_g[(h w)^2] = e^25 P(X > 10), so sd(h w) = sqrt(e^25 Q10 - p^2)
const IS_SD = Math.sqrt(Math.exp(25) * Q10 - RARE_EXACT ** 2);   // relative sd 2.38

export interface Problem {
  key: ProbKey;
  label: string;
  exact: number;
  sd: number;                         // true sd of one draw: the SE is sd / sqrt(N)
  draw: (r: Rand) => number;          // one Monte Carlo term
  hit?: (y: number) => boolean;       // counts "events" (indicator > 0) for the readout
}

export function problem(key: ProbKey, is: boolean): Problem {
  if (key === "auc") return {
    key, label: "area under the OD curve, 0–10 h", exact: AUC_EXACT, sd: AUC_SD,
    draw: r => 10 * od(10 * r.unif()),
  };
  if (key === "motif") return {
    key, label: "P(K ≥ 18) of 200 random motif hits", exact: MOTIF_EXACT,
    sd: Math.sqrt(MOTIF_EXACT * (1 - MOTIF_EXACT)),
    draw: r => (overlaps(r) >= KOBS ? 1 : 0), hit: y => y > 0,
  };
  if (!is) return {
    key, label: "P(X > 5), plain draws from N(0, 1)", exact: RARE_EXACT,
    sd: Math.sqrt(RARE_EXACT * (1 - RARE_EXACT)),
    draw: r => (r.norm() > 5 ? 1 : 0), hit: y => y > 0,
  };
  return {
    key, label: "P(X > 5), importance sampling from N(5, 1)", exact: RARE_EXACT, sd: IS_SD,
    draw: r => { const x = 5 + r.norm(); return x > 5 ? Math.exp(12.5 - 5 * x) : 0; },
    hit: y => y > 0,
  };
}

/** Running mean and variance (Welford), with checkpoints on a geometric grid of N. */
export class Running {
  n = 0; mean = 0; m2 = 0; hits = 0; sumw = 0; sumw2 = 0;
  hist: { n: number; est: number; se: number }[] = [];
  private next = 1;
  add(y: number, hit: boolean): void {
    this.n++;
    const d = y - this.mean;
    this.mean += d / this.n;
    this.m2 += d * (y - this.mean);
    if (hit) this.hits++;
    if (this.n >= this.next) { this.mark(); this.next = Math.max(this.n + 1, Math.floor(this.n * 1.06)); }
  }
  get se(): number { return this.n > 1 ? Math.sqrt(this.m2 / (this.n - 1) / this.n) : 0; }
  mark(): void {
    const last = this.hist[this.hist.length - 1];
    if (this.n > 0 && (!last || last.n !== this.n)) this.hist.push({ n: this.n, est: this.mean, se: this.se });
  }
}

// ---- curse of dimensionality: the toy marginal likelihood of the lecture ----
// d parameters, prior N(0, 1), one measurement y_i = 1 with noise sd 1: Z = z1^d.
const SQ2PI = Math.sqrt(2 * Math.PI);
const dnorm = (x: number): number => Math.exp(-0.5 * x * x) / SQ2PI;
export const Z1 = Math.exp(-0.25) / (2 * Math.sqrt(Math.PI));    // dnorm(1, 0, sqrt 2) = 0.2196956
/** E[L^2]/Z^2 per dimension = 2 e^(1/6) / sqrt(3) = 1.36413. */
export const VAR_RATIO = 2 * Math.exp(1 / 6) / Math.sqrt(3);

/** Inverse normal CDF (Acklam's rational approximation, relative error < 1.2e-9). */
export function qnorm(p: number): number {
  const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
  const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
  const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
  const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
  const pl = 0.02425;
  if (p < pl) {
    const q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  }
  if (p > 1 - pl) return -qnorm(1 - p);
  const q = p - 0.5, r = q * q;
  return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
}

/** Per-axis relative error of the grid (n prior quantiles), e1(n) = mean L / z1 - 1. Memoised. */
export const NMAX = 3000;
const e1cache = new Map<number, number>();
export function e1(n: number): number {
  let v = e1cache.get(n);
  if (v === undefined) {
    let s = 0;
    for (let i = 1; i <= n; i++) s += dnorm(1 - qnorm((i - 0.5) / n));
    v = s / n / Z1 - 1;
    e1cache.set(n, v);
  }
  return v;
}

/** Points per axis that a budget N allows in d dimensions. */
export const perAxis = (N: number, d: number): number => Math.max(1, Math.floor(Math.pow(N, 1 / d) * (1 + 1e-12)));

/** Relative error of the grid estimate of Z with n points per axis (the grid is a tensor product, so Z_grid = (z1 (1 + e1))^d). */
export const gridErr = (n: number, d: number): number => Math.abs(Math.pow(1 + e1(n), d) - 1);

/** Expected relative SE of plain Monte Carlo with N prior draws. */
export const mcRelSE = (N: number, d: number): number => Math.sqrt(Math.pow(VAR_RATIO, d) - 1) / Math.sqrt(N);

/** One Monte Carlo run with N draws from the prior: relative error and relative SE. */
export function mcRun(N: number, d: number, seed: number): { err: number; relse: number } {
  const r = rng(1000 * seed + d);
  const c = -d * Math.log(SQ2PI);
  let mean = 0, m2 = 0;
  for (let k = 1; k <= N; k++) {
    let ll = c;
    for (let i = 0; i < d; i++) { const e = 1 - r.norm(); ll -= 0.5 * e * e; }
    const y = Math.exp(ll);
    const dd = y - mean; mean += dd / k; m2 += dd * (y - mean);
  }
  const Z = Math.pow(Z1, d);
  return { err: Math.abs(mean / Z - 1), relse: Math.sqrt(m2 / (N - 1) / N) / Z };
}

/** log10 of the grid evaluations needed for a relative error <= tol (Infinity if n > NMAX). */
export function gridNeed(d: number, tol: number): { n: number; log10N: number } {
  for (let n = 1; n <= NMAX; n++) if (gridErr(n, d) <= tol) return { n, log10N: d * Math.log10(n) };
  return { n: Infinity, log10N: Infinity };
}

/** log10 of the Monte Carlo draws needed for a relative SE <= tol. */
export const mcNeed = (d: number, tol: number): number => Math.log10((Math.pow(VAR_RATIO, d) - 1) / (tol * tol));

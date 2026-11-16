// Taylor polynomials around 0 for three functions, with the remainder-term bound.
// Kept free of DOM code, so the numbers on the exercise page can be checked in node.

export type FnKey = "exp" | "log" | "sqrt";

export interface TaylorFn {
  key: FnKey;
  label: string;              // as shown on buttons and axes
  f: (x: number) => number;   // computed the plain way, as a student would in R
  coef: number[];             // c_0 .. c_7, P_n(x) = sum c_k x^k
  terms: string[];            // the k-th term as text, with its sign ("+ x²/2")
  /** g(ξ) in R_n = c_{n+1} g(ξ) x^{n+1}; g(0) = 1 and g is monotone in ξ. */
  g: (xi: number, n: number) => number;
  xmin: number;               // left end of the domain we allow for x
  ylim: [number, number];     // y range of panel (a)
}

const SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹";
export const sup = (k: number): string => String(k).split("").map(d => SUP[Number(d)]).join("");
const xk = (k: number): string => (k === 0 ? "" : k === 1 ? "x" : "x" + sup(k));

function fact(k: number): number { let p = 1; for (let i = 2; i <= k; i++) p *= i; return p; }

const expCoef = Array.from({ length: 8 }, (_, k) => 1 / fact(k));
const logCoef = Array.from({ length: 8 }, (_, k) => (k === 0 ? 0 : (k % 2 ? 1 : -1) / k));
const sqrtCoef: number[] = [1];
for (let k = 0; k < 7; k++) sqrtCoef.push(sqrtCoef[k] * (0.5 - k) / (k + 1));   // binom(1/2, k)
const sqrtFrac: [number, number][] = [[1, 1], [1, 2], [-1, 8], [1, 16], [-5, 128], [7, 256], [-21, 1024]];

function signed(first: boolean, neg: boolean, body: string): string {
  return first ? (neg ? "−" + body : body) : (neg ? " − " : " + ") + body;
}

export const FNS: Record<FnKey, TaylorFn> = {
  exp: {
    key: "exp", label: "eˣ", f: Math.exp, coef: expCoef,
    terms: Array.from({ length: 7 }, (_, k) => signed(k === 0, false, k === 0 ? "1" : k === 1 ? "x" : `${xk(k)}/${fact(k)}`)),
    g: xi => Math.exp(xi),
    xmin: -1, ylim: [-1, 7.5],
  },
  log: {
    key: "log", label: "log(1+x)", f: x => Math.log(1 + x), coef: logCoef,
    terms: Array.from({ length: 7 }, (_, k) => (k === 0 ? "0" : signed(k === 1, k % 2 === 0, k === 1 ? "x" : `${xk(k)}/${k}`))),
    g: (xi, n) => Math.pow(1 + xi, -(n + 1)),
    xmin: -0.9, ylim: [-2.5, 1.75],
  },
  sqrt: {
    key: "sqrt", label: "√(1+x)", f: x => Math.sqrt(1 + x), coef: sqrtCoef,
    terms: sqrtFrac.map(([p, q], k) => {
      if (k === 0) return "1";
      const a = Math.abs(p);
      const body = (a === 1 ? xk(k) : `${a}${xk(k)}`) + (q === 1 ? "" : `/${q}`);
      return signed(false, p < 0, body);
    }),
    g: (xi, n) => Math.pow(1 + xi, 0.5 - (n + 1)),
    xmin: -0.9, ylim: [-0.25, 2],
  },
};

/** P_n(x) by Horner's rule. */
export function taylor(fn: TaylorFn, n: number, x: number): number {
  let p = fn.coef[n];
  for (let k = n - 1; k >= 0; k--) p = p * x + fn.coef[k];
  return p;
}

/** The formula of P_n as text. For log(1+x), the 0 term is dropped once n ≥ 1. */
export function formula(fn: TaylorFn, n: number): string {
  if (fn.key === "log" && n >= 1) return fn.terms.slice(1, n + 1).join("");
  return fn.terms.slice(0, n + 1).join("");
}

/** Computed error f(x) − P_n(x), in double precision. */
export function err(fn: TaylorFn, n: number, x: number): number {
  return fn.f(x) - taylor(fn, n, x);
}

/** Bounds on |R_n(x)| from the remainder term, ξ between 0 and x: [lo, hi]. */
export function bound(fn: TaylorFn, n: number, x: number): [number, number] {
  const t = Math.abs(fn.coef[n + 1]) * Math.pow(Math.abs(x), n + 1);   // first omitted term
  const a = fn.g(0, n), b = fn.g(x, n);
  return [t * Math.min(a, b), t * Math.max(a, b)];
}

/** First omitted term |c_{n+1}| |x|^{n+1}: a line with slope n+1 on log–log axes. */
export function nextTerm(fn: TaylorFn, n: number, x: number): number {
  return Math.abs(fn.coef[n + 1]) * Math.pow(Math.abs(x), n + 1);
}

/** Measured log–log slope of |error| between x/10 and x. */
export function slope(fn: TaylorFn, n: number, x: number): number {
  const e1 = Math.abs(err(fn, n, x)), e0 = Math.abs(err(fn, n, x / 10));
  return Math.log10(e1 / e0);
}

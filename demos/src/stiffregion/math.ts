// Numbers for the stiffregion demo: amplification factors R(z) of four one-step methods,
// stability limits along a ray, and the closed-form numerical solution of a 2x2 linear
// system y' = A y whose eigenvalues are either real (slow, fast) or a complex pair.
// Because the methods are linear, n steps multiply each eigen-mode by R(h lambda)^n:
// no time stepping loop is needed, so blow-up and step counts are exact and instant.

export type C = [number, number];
export const cadd = (a: C, b: C): C => [a[0] + b[0], a[1] + b[1]];
export const cmul = (a: C, b: C): C => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
export const cdiv = (a: C, b: C): C => {
  const d = b[0] * b[0] + b[1] * b[1];
  return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d];
};
export const cabs = (a: C): number => Math.hypot(a[0], a[1]);

export type Method = "euler" | "rk4" | "beuler" | "trap";
export const METHODS: Method[] = ["euler", "rk4", "beuler", "trap"];
export const LABEL: Record<Method, string> = {
  euler: "Euler", rk4: "RK4", beuler: "backward Euler", trap: "trapezoid",
};
export const RFORM: Record<Method, string> = {
  euler: "R = 1 + z", rk4: "R = 1 + z + z²/2 + z³/6 + z⁴/24",
  beuler: "R = 1/(1 − z)", trap: "R = (1 + z/2)/(1 − z/2)",
};

/** Amplification factor R(z) per step on y' = lambda y, z = h lambda. */
export function R(m: Method, z: C): C {
  switch (m) {
    case "euler": return [1 + z[0], z[1]];
    case "rk4": {
      // Horner: 1 + z(1 + z/2(1 + z/3(1 + z/4)))
      let p: C = [1 + z[0] / 4, z[1] / 4];
      p = cadd([1, 0], cmul([z[0] / 3, z[1] / 3], p));
      p = cadd([1, 0], cmul([z[0] / 2, z[1] / 2], p));
      return cadd([1, 0], cmul(z, p));
    }
    case "beuler": return cdiv([1, 0], [1 - z[0], -z[1]]);
    case "trap": return cdiv([1 + z[0] / 2, z[1] / 2], [1 - z[0] / 2, -z[1] / 2]);
  }
}
/** Real R for real z (all four methods map real z to real R). */
export const Rr = (m: Method, z: number): number => R(m, [z, 0])[0];

/**
 * Largest s > 0 with |R(s d)| <= 1 for all of (0, s], d a unit direction.
 * Infinity if the whole ray is stable; 0 if it is unstable from the start.
 */
export function rayLimit(m: Method, d: C): number {
  const bad = (s: number) => cabs(R(m, [s * d[0], s * d[1]])) > 1 + 1e-9;
  if (m === "beuler" || m === "trap") {
    // A-stable: the whole left half-plane (and the imaginary axis) is stable
    if (d[0] <= 1e-12) return Infinity;
  }
  let lo = 0, s = 1e-3;
  if (bad(s)) return 0;
  while (s < 1e3 && !bad(s)) { lo = s; s *= 1.05; }
  if (s >= 1e3) return Infinity;
  let hi = s;
  for (let i = 0; i < 60; i++) { const mid = (lo + hi) / 2; if (bad(mid)) hi = mid; else lo = mid; }
  return lo;
}

/** Boundary of a bounded stability region, polar around a centre (for Euler and RK4). */
export function regionPolygon(m: Method, cx: number, nAng = 360): C[] {
  const pts: C[] = [];
  for (let k = 0; k < nAng; k++) {
    const a = 2 * Math.PI * k / nAng, d: C = [Math.cos(a), Math.sin(a)];
    const ok = (r: number) => cabs(R(m, [cx + r * d[0], r * d[1]])) <= 1;
    let lo = 0, r = 0.01;
    while (r < 10 && ok(r)) { lo = r; r += 0.01; }
    let hi = r;
    for (let i = 0; i < 40; i++) { const mid = (lo + hi) / 2; if (ok(mid)) lo = mid; else hi = mid; }
    pts.push([cx + lo * d[0], lo * d[1]]);
  }
  return pts;
}

// ---------------------------------------------------------------- the test problem
/**
 * Real case: eigenvalues ls (slow) and lf (fast), y1 = e^{ls t} + e^{lf t},
 * y2 = e^{ls t} − e^{lf t} (A = V diag(ls, lf) V^{-1}, V = [[1, 1], [1, −1]], y0 = (2, 0)):
 * y2 rises fast from 0, then decays slowly, like the enzyme complex C.
 * Complex case: lambda = a ± i w, y1 + i y2 = e^{lambda t} (A = [[a, −w], [w, a]], y0 = (1, 0)).
 */
export interface Problem {
  key: string;
  complex: boolean;
  ls: number;       // real case: slow eigenvalue; complex case: real part a
  lf: number;       // real case: fast eigenvalue; complex case: imaginary part w
  T: number;        // time span
  unit: string;     // time unit
}

/** Eigenvalues shown in the plane (two points). */
export function eigs(p: Problem): [C, C] {
  return p.complex ? [[p.ls, p.lf], [p.ls, -p.lf]] : [[p.ls, 0], [p.lf, 0]];
}

export function exact(p: Problem, t: number): [number, number] {
  if (p.complex) {
    const r = Math.exp(p.ls * t);
    return [r * Math.cos(p.lf * t), r * Math.sin(p.lf * t)];
  }
  const a = Math.exp(p.ls * t), b = Math.exp(p.lf * t);
  return [a + b, a - b];
}

/** Numerical solution after n steps of size h (closed form). */
export function numeric(p: Problem, m: Method, h: number, n: number): [number, number] {
  if (p.complex) {
    const r = R(m, [h * p.ls, h * p.lf]);
    const mod = Math.pow(cabs(r), n), ang = n * Math.atan2(r[1], r[0]);
    return [mod * Math.cos(ang), mod * Math.sin(ang)];
  }
  const a = Math.pow(Rr(m, h * p.ls), n), b = Math.pow(Rr(m, h * p.lf), n);
  return [a + b, a - b];
}

/** Error at t = T with N equal steps (max over the two components). */
export function errAtT(p: Problem, m: Method, N: number): number {
  const h = p.T / N;
  const y = numeric(p, m, h, N), e = exact(p, p.T);
  const v = Math.max(Math.abs(y[0] - e[0]), Math.abs(y[1] - e[1]));
  return Number.isFinite(v) ? v : Infinity;
}

export const NMAX = 1e8;
/** Smallest number of equal steps N with error at T <= tol (scan, then integer bisection). */
export function stepsNeeded(p: Problem, m: Method, tol: number): number {
  let prev = 0;
  for (let N = 1; N <= NMAX; N = Math.max(N + 1, Math.floor(N * 1.02))) {
    if (errAtT(p, m, N) <= tol) {
      let lo = prev, hi = N;           // lo fails (or 0), hi passes
      while (hi - lo > 1) {
        const mid = Math.floor((lo + hi) / 2);
        if (errAtT(p, m, mid) <= tol) hi = mid; else lo = mid;
      }
      return hi;
    }
    prev = N;
  }
  return Infinity;
}

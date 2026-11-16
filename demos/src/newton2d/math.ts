// The toggle switch of the newton-systems unit and Newton's method for it.
//   F1(u, v) = alpha / (1 + v^n) - u,   F2(u, v) = alpha / (1 + u^n) - v
//   J = [ -1  b ; c  -1 ],  b = -alpha n v^(n-1) / (1 + v^n)^2,  c = -alpha n u^(n-1) / (1 + u^n)^2

export interface Par { alpha: number; n: number }
export type Vec = [number, number];

export function F(p: Par, x: Vec): Vec {
  return [p.alpha / (1 + x[1] ** p.n) - x[0], p.alpha / (1 + x[0] ** p.n) - x[1]];
}

/** Jacobian as [J11, J12, J21, J22]. */
export function J(p: Par, x: Vec): [number, number, number, number] {
  const [u, v] = x, a = p.alpha, n = p.n;
  const b = -a * n * v ** (n - 1) / (1 + v ** n) ** 2;
  const c = -a * n * u ** (n - 1) / (1 + u ** n) ** 2;
  return [-1, b, c, -1];
}

export const det = (j: number[]): number => j[0] * j[3] - j[1] * j[2];

/** Solve J d = -F for a 2x2 system (Cramer's rule). */
export function step(p: Par, x: Vec): { f: Vec; j: number[]; d: Vec; det: number } {
  const f = F(p, x), j = J(p, x), D = det(j);
  const d: Vec = [(-f[0] * j[3] + f[1] * j[1]) / D, (-f[1] * j[0] + f[0] * j[2]) / D];
  return { f, j, d, det: D };
}

export const norm = (a: Vec): number => Math.hypot(a[0], a[1]);

export const TOL = 1e-12, MAXIT = 50;

export interface Path { xs: Vec[]; status: "converged" | "maxit" | "nan" | "escaped" }

/** Newton iterates from x0; stops when ||F(x_k)|| < TOL, after MAXIT steps, on NaN, or if |x| > 1e6. */
export function newton(p: Par, x0: Vec): Path {
  const xs: Vec[] = [x0];
  let x = x0;
  for (let k = 0; ; k++) {
    const f = F(p, x);
    if (!Number.isFinite(f[0]) || !Number.isFinite(f[1])) return { xs, status: "nan" };
    if (norm(f) < TOL) return { xs, status: "converged" };
    if (Math.abs(x[0]) > 1e6 || Math.abs(x[1]) > 1e6) return { xs, status: "escaped" };
    if (k >= MAXIT) return { xs, status: "maxit" };
    const s = step(p, x);
    x = [x[0] + s.d[0], x[1] + s.d[1]];
    xs.push(x);
  }
}

export interface Root { x: Vec; kind: "hi" | "mid" | "lo" }

/** All steady states with u, v >= 0: roots of u = alpha/(1 + g(u)^n), g(u) = alpha/(1 + u^n). */
export function roots(p: Par): Root[] {
  const g = (u: number) => p.alpha / (1 + u ** p.n);
  const h = (u: number) => g(g(u)) - u;
  const out: Root[] = [];
  const N = 4000, hi = p.alpha;
  let u0 = 0, h0 = h(0);
  for (let i = 1; i <= N; i++) {
    const u1 = hi * i / N, h1 = h(u1);
    if (h1 === 0) { out.push(mk(u1)); }
    else if (h0 * h1 < 0) {
      let a = u0, b = u1, ha = h0;
      for (let it = 0; it < 200 && b - a > 1e-16 * Math.max(1, b); it++) {
        const m = 0.5 * (a + b), hm = h(m);
        if (hm === 0) { a = b = m; break; }
        if (ha * hm < 0) b = m; else { a = m; ha = hm; }
      }
      out.push(mk(0.5 * (a + b)));
    }
    u0 = u1; h0 = h1;
  }
  function mk(u: number): Root {
    const v = g(u);
    // polish with a few Newton steps on the full system (exact to the last digits)
    let x: Vec = [u, v];
    for (let k = 0; k < 4; k++) {
      const s = step(p, x);
      if (!Number.isFinite(s.d[0]) || !Number.isFinite(s.d[1]) || norm(s.d) > 1e-6) break;
      x = [x[0] + s.d[0], x[1] + s.d[1]];
    }
    const kind = Math.abs(x[0] - x[1]) < 1e-7 ? "mid" : x[0] > x[1] ? "hi" : "lo";
    return { x, kind };
  }
  // drop duplicates (sign changes at the same root), sort U high, middle, V high
  const uniq: Root[] = [];
  for (const r of out) if (!uniq.some(q => norm([q.x[0] - r.x[0], q.x[1] - r.x[1]]) < 1e-7)) uniq.push(r);
  return uniq.sort((a, b) => (b.x[0] - b.x[1]) - (a.x[0] - a.x[1]));
}

/** Index of the root that Newton from x0 reaches, or -1 (no convergence, NaN, escape, other). */
export function basinOf(p: Par, rs: Root[], x0: Vec): number {
  let u = x0[0], v = x0[1];
  const a = p.alpha, n = p.n, n1 = n - 1;
  for (let k = 0; k <= MAXIT; k++) {
    const vn = v ** n, un = u ** n;
    const f1 = a / (1 + vn) - u, f2 = a / (1 + un) - v;
    if (!(Math.abs(f1) < 1e300 && Math.abs(f2) < 1e300)) return -1;   // NaN or Inf
    if (f1 * f1 + f2 * f2 < TOL * TOL) {
      let best = -1, bd = 1e-3;
      for (let i = 0; i < rs.length; i++) {
        const d = Math.hypot(u - rs[i].x[0], v - rs[i].x[1]);
        if (d < bd) { bd = d; best = i; }
      }
      return best;
    }
    if (Math.abs(u) > 1e6 || Math.abs(v) > 1e6) return -1;
    const b = -a * n * v ** n1 / ((1 + vn) * (1 + vn));
    const c = -a * n * u ** n1 / ((1 + un) * (1 + un));
    const D = 1 - b * c;                     // det of [-1 b; c -1]
    const du = (f1 + b * f2) / D, dv = (c * f1 + f2) / D;
    u += du; v += dv;
  }
  return -1;
}

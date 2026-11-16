// Models, fixed points, Jacobians, eigenvalues and RK4 for the phase-plane demo.
//   toggle:  du/dt = au / (1 + v^2) - u,   dv/dt = av / (1 + u^2) - v    (n = 2)
//   selkov:  du/dt = -u + a v + u^2 v,     dv/dt = b - a v - u^2 v

export type Vec = [number, number];
export type Model = "toggle" | "selkov";
export interface Par { model: Model; au: number; av: number; a: number; b: number }

export function F(p: Par, x: Vec): Vec {
  const [u, v] = x;
  if (p.model === "toggle") return [p.au / (1 + v * v) - u, p.av / (1 + u * u) - v];
  const r = u * u * v;
  return [-u + p.a * v + r, p.b - p.a * v - r];
}

/** Jacobian [j11, j12, j21, j22]. */
export function J(p: Par, x: Vec): [number, number, number, number] {
  const [u, v] = x;
  if (p.model === "toggle") {
    const dv = -2 * p.au * v / (1 + v * v) ** 2, du = -2 * p.av * u / (1 + u * u) ** 2;
    return [-1, dv, du, -1];
  }
  return [-1 + 2 * u * v, p.a + u * u, -2 * u * v, -p.a - u * u];
}

export interface Eig { tr: number; det: number; re: [number, number]; im: number; type: string; stable: boolean; saddle: boolean }

/** Eigenvalues of a 2x2 and the type of the fixed point. */
export function eig(j: [number, number, number, number]): Eig {
  const tr = j[0] + j[3], det = j[0] * j[3] - j[1] * j[2];
  const disc = tr * tr / 4 - det;
  let re: [number, number], im = 0;
  if (disc >= 0) {
    const s = Math.sqrt(disc);
    re = [tr / 2 + s, tr / 2 - s];
  } else {
    re = [tr / 2, tr / 2]; im = Math.sqrt(-disc);
  }
  const ZERO = 1e-3;
  let type: string;
  if (det < 0) type = "saddle";
  else if (Math.abs(tr / 2) < ZERO && im > 0) type = "centre? (Re λ ≈ 0)";
  else if (Math.abs(re[0]) < ZERO || Math.abs(re[1]) < ZERO) type = "undecided (λ ≈ 0)";
  else if (im > 0) type = tr < 0 ? "stable spiral" : "unstable spiral";
  else type = tr < 0 ? "stable node" : "unstable node";
  return { tr, det, re, im, type, stable: det > 0 && tr < 0, saddle: det < 0 };
}

export interface Fixed { x: Vec; e: Eig; name: string }

/** All fixed points (toggle: scan + bisection along u; Sel'kov: closed form). */
export function fixedPoints(p: Par): Fixed[] {
  const out: Vec[] = [];
  if (p.model === "selkov") {
    out.push([p.b, p.b / (p.a + p.b * p.b)]);
  } else {
    const vOf = (u: number) => p.av / (1 + u * u);
    const h = (u: number) => { const v = vOf(u); return u - p.au / (1 + v * v); };
    const N = 6000, hi = p.au + 0.5;
    let u0 = 0, h0 = h(0);
    for (let i = 1; i <= N; i++) {
      const u1 = hi * i / N, h1 = h(u1);
      if (h0 === 0) out.push([u0, vOf(u0)]);
      else if (h0 * h1 < 0) {
        let lo = u0, up = u1, hl = h0;
        for (let k = 0; k < 80; k++) {
          const m = (lo + up) / 2, hm = h(m);
          if (hm * hl <= 0) up = m; else { lo = m; hl = hm; }
        }
        const u = (lo + up) / 2;
        out.push([u, vOf(u)]);
      }
      u0 = u1; h0 = h1;
    }
  }
  const fs = out.map(x => ({ x, e: eig(J(p, x)), name: "" }));
  if (p.model === "toggle") {
    if (fs.length === 1) fs[0].name = "the only state";
    else fs.forEach((f, i) => { f.name = i === 0 ? "V on" : i === fs.length - 1 ? "U on" : "middle"; });
  } else fs[0].name = "fixed point";
  return fs;
}

/** Eigenvector of a 2x2 for a real eigenvalue lam, normalised. */
export function evec(j: [number, number, number, number], lam: number): Vec {
  let x: Vec = Math.abs(j[1]) > Math.abs(j[2]) ? [j[1], lam - j[0]] : [lam - j[3], j[2]];
  if (Math.hypot(x[0], x[1]) === 0) x = [1, 0];
  const n = Math.hypot(x[0], x[1]);
  return [x[0] / n, x[1] / n];
}

/** One classical RK4 step of dx/dt = s F(x). */
export function rk4(p: Par, x: Vec, h: number, s = 1): Vec {
  const k1 = F(p, x);
  const k2 = F(p, [x[0] + h / 2 * s * k1[0], x[1] + h / 2 * s * k1[1]]);
  const k3 = F(p, [x[0] + h / 2 * s * k2[0], x[1] + h / 2 * s * k2[1]]);
  const k4 = F(p, [x[0] + h * s * k3[0], x[1] + h * s * k3[1]]);
  return [x[0] + s * h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
          x[1] + s * h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])];
}

export interface Traj { t: number[]; x: Vec[]; settled: number; end: number }

/**
 * RK4 with step h to tmax. `settled` is the first index from which the trajectory stays
 * within tol of a fixed point `end` (or -1). A run stops early at a stable fixed point; an
 * unstable one counts only if the trajectory is still there at tmax (passing by is not settling).
 */
export function trajectory(p: Par, x0: Vec, h: number, tmax: number, fps: Fixed[], tol = 1e-3): Traj {
  const t = [0], x: Vec[] = [x0];
  let cur = x0, settled = -1, end = -1;
  const near = (y: Vec) => fps.findIndex(f => Math.hypot(y[0] - f.x[0], y[1] - f.x[1]) < tol);
  const n = Math.round(tmax / h);
  for (let i = 1; i <= n; i++) {
    cur = rk4(p, cur, h);
    if (!Number.isFinite(cur[0]) || !Number.isFinite(cur[1])) break;
    t.push(i * h); x.push(cur);
    const k = near(cur);
    if (k >= 0 && fps[k].e.stable) return { t, x, settled: i, end: k };
  }
  const k = near(x[x.length - 1]);
  if (k >= 0) {
    end = k; settled = x.length - 1;
    const f = fps[k].x;
    while (settled > 0 && Math.hypot(x[settled - 1][0] - f[0], x[settled - 1][1] - f[1]) < tol) settled--;
  }
  return { t, x, settled, end };
}

/** Separatrix of a saddle: two backward runs from saddle +- 1e-6 v_stable. */
export function separatrix(p: Par, s: Fixed, box: [number, number, number, number]): Vec[] {
  const j = J(p, s.x), lam = Math.min(s.e.re[0], s.e.re[1]);
  const v = evec(j, lam);
  const arms: Vec[][] = [];
  for (const sg of [1, -1]) {
    let x: Vec = [s.x[0] + sg * 1e-6 * v[0], s.x[1] + sg * 1e-6 * v[1]];
    const arm: Vec[] = [x];
    for (let i = 0; i < 3000; i++) {
      x = rk4(p, x, 0.005, -1);
      if (!Number.isFinite(x[0]) || !Number.isFinite(x[1])) break;
      arm.push(x);
      if (x[0] < box[0] || x[0] > box[1] || x[1] < box[2] || x[1] > box[3]) break;
    }
    arms.push(arm);
  }
  return [...arms[0].reverse(), s.x, ...arms[1]];
}

/** Fate of a start for the basin map: index of the stable state reached (-1 none; saddle index if stuck). */
export function fate(p: Par, x0: Vec, fps: Fixed[]): number {
  let x = x0;
  const h = 0.05;
  for (let i = 0; i < 2000; i++) {
    x = rk4(p, x, h);
    if (i % 4 === 3) {
      for (let k = 0; k < fps.length; k++) {
        const tol = fps[k].e.stable ? 0.02 : 1e-6;
        if (Math.abs(x[0] - fps[k].x[0]) < tol && Math.abs(x[1] - fps[k].x[1]) < tol) return k;
      }
    }
  }
  let best = -1, bd = Infinity;
  fps.forEach((f, k) => { const d = Math.hypot(x[0] - f.x[0], x[1] - f.x[1]); if (d < bd) { bd = d; best = k; } });
  return best;
}

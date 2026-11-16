// Numerics for the turing demo, kept free of DOM code so it can be tested in node.
//  - a mixed-radix complex FFT (any N; the demo uses 40, 80, 160 = 2^k * 5),
//  - the unit's 2D Schnakenberg IMEX step on a periodic grid: reaction explicit,
//    five-point diffusion implicit by FFT (divide by 1 + dt D kappa),
//  - the unit's Park-Miller generator lcg() and spot counter,
//  - the 1D synthesis-diffusion-degradation model: IMEX step with Thomas, exact steady state.

// ------------------------------------------------------------------ FFT
export interface Plan { n: number; fac: number[]; cos: Float64Array; sin: Float64Array; tr: Float64Array; ti: Float64Array }

export function makePlan(n: number): Plan {
  const fac: number[] = [];
  let m = n;
  for (const p of [5, 3, 4, 2]) while (m % p === 0) { fac.push(p); m /= p; }
  for (let p = 7; m > 1; p += 2) while (m % p === 0) { fac.push(p); m /= p; }
  const cos = new Float64Array(n), sin = new Float64Array(n);
  for (let k = 0; k < n; k++) { cos[k] = Math.cos(2 * Math.PI * k / n); sin[k] = Math.sin(2 * Math.PI * k / n); }
  const pm = Math.max(...fac, 1);
  return { n, fac, cos, sin, tr: new Float64Array(pm), ti: new Float64Array(pm) };
}

// x (stride sx) -> y (contiguous), length n, sign -1 forward, +1 inverse (unnormalised)
function rec(P: Plan, xr: Float64Array, xi: Float64Array, xo: number, sx: number, n: number,
  yr: Float64Array, yi: Float64Array, yo: number, fi: number, sg: number): void {
  if (n === 1) { yr[yo] = xr[xo]; yi[yo] = xi[xo]; return; }
  const p = P.fac[fi], m = n / p, N = P.n, st = N / n;
  for (let j = 0; j < p; j++) rec(P, xr, xi, xo + j * sx, sx * p, m, yr, yi, yo + j * m, fi + 1, sg);
  const tr = P.tr, ti = P.ti, C = P.cos, S = P.sin;
  if (p === 2) {
    for (let k = 0; k < m; k++) {
      const a = yo + k, b = a + m, w = k * st;
      const c = C[w], s = sg * S[w];
      const br = yr[b] * c - yi[b] * s, bi = yr[b] * s + yi[b] * c;
      const ar = yr[a], ai = yi[a];
      yr[a] = ar + br; yi[a] = ai + bi; yr[b] = ar - br; yi[b] = ai - bi;
    }
    return;
  }
  if (p === 4) {
    for (let k = 0; k < m; k++) {
      const i0 = yo + k, i1 = i0 + m, i2 = i1 + m, i3 = i2 + m;
      const w1 = k * st, w2 = 2 * w1, w3 = 3 * w1;
      const c1 = C[w1], s1 = sg * S[w1], c2 = C[w2], s2 = sg * S[w2], c3 = C[w3], s3 = sg * S[w3];
      const t0r = yr[i0], t0i = yi[i0];
      const t1r = yr[i1] * c1 - yi[i1] * s1, t1i = yr[i1] * s1 + yi[i1] * c1;
      const t2r = yr[i2] * c2 - yi[i2] * s2, t2i = yr[i2] * s2 + yi[i2] * c2;
      const t3r = yr[i3] * c3 - yi[i3] * s3, t3i = yr[i3] * s3 + yi[i3] * c3;
      const ar = t0r + t2r, ai = t0i + t2i, br = t0r - t2r, bi = t0i - t2i;
      const cr = t1r + t3r, ci = t1i + t3i, er = t1r - t3r, ei = t1i - t3i;
      // sg * i * e = (-sg ei, sg er)
      yr[i0] = ar + cr; yi[i0] = ai + ci; yr[i2] = ar - cr; yi[i2] = ai - ci;
      yr[i1] = br - sg * ei; yi[i1] = bi + sg * er; yr[i3] = br + sg * ei; yi[i3] = bi - sg * er;
    }
    return;
  }
  const NP = N / p;
  for (let k = 0; k < m; k++) {
    for (let j = 0; j < p; j++) {          // twiddle Y_j[k] by W_n^{jk}
      const id = yo + j * m + k, w = j * k * st;   // < N since j < p, k < m
      const c = C[w], s = sg * S[w];
      tr[j] = yr[id] * c - yi[id] * s; ti[j] = yr[id] * s + yi[id] * c;
    }
    for (let q = 0; q < p; q++) {          // length-p DFT
      let sr = 0, si = 0;
      let w = 0;
      const dw = q * NP;
      for (let j = 0; j < p; j++) {
        const c = C[w], s = sg * S[w];
        sr += tr[j] * c - ti[j] * s; si += tr[j] * s + ti[j] * c;
        w += dw; if (w >= N) w -= N;
      }
      yr[yo + q * m + k] = sr; yi[yo + q * m + k] = si;
    }
  }
}

/** 2D FFT in place of an N x N complex array (column-major), via scratch lines. */
export function fft2(P: Plan, re: Float64Array, im: Float64Array, sg: number,
  br: Float64Array, bi: Float64Array, cr: Float64Array, ci: Float64Array): void {
  const n = P.n;
  for (let j = 0; j < n; j++) {               // columns: contiguous
    for (let i = 0; i < n; i++) { br[i] = re[i + n * j]; bi[i] = im[i + n * j]; }
    rec(P, br, bi, 0, 1, n, cr, ci, 0, 0, sg);
    for (let i = 0; i < n; i++) { re[i + n * j] = cr[i]; im[i + n * j] = ci[i]; }
  }
  for (let i = 0; i < n; i++) {               // rows: stride n
    for (let j = 0; j < n; j++) { br[j] = re[i + n * j]; bi[j] = im[i + n * j]; }
    rec(P, br, bi, 0, 1, n, cr, ci, 0, 0, sg);
    for (let j = 0; j < n; j++) { re[i + n * j] = cr[j]; im[i + n * j] = ci[j]; }
  }
}

// ------------------------------------------------------------------ the unit's generator
export function lcg(n: number, seed = 1): Float64Array {
  const z = new Float64Array(n);
  let s = seed;
  for (let i = 0; i < n; i++) { s = (16807 * s) % 2147483647; z[i] = s / 2147483647; }
  return z;
}

// ------------------------------------------------------------------ Schnakenberg 2D
export interface Turing {
  N: number; dx: number; dt: number; d: number; a: number; b: number; t: number; steps: number;
  u: Float64Array; v: Float64Array;
  P: Plan; inv1: Float64Array; invd: Float64Array;
  wr: Float64Array; wi: Float64Array; br: Float64Array; bi: Float64Array; cr: Float64Array; ci: Float64Array;
}

/** New simulation on an N x N periodic grid of side Lbox, noise 0.02 (lcg - 0.5) as in the unit. */
export function makeTuring(N: number, Lbox: number, dt: number, d: number, a: number, b: number, seed: number): Turing {
  const dx = Lbox / N, P = makePlan(N);
  const s = new Float64Array(N);
  for (let k = 0; k < N; k++) s[k] = 4 / (dx * dx) * Math.sin(Math.PI * k / N) ** 2;
  const inv1 = new Float64Array(N * N), invd = new Float64Array(N * N);
  for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
    const kap = s[i] + s[j];
    inv1[i + N * j] = 1 / (1 + dt * kap); invd[i + N * j] = 1 / (1 + dt * d * kap);
  }
  const us = a + b, vs = b / (us * us);
  const z = lcg(2 * N * N, seed);
  const u = new Float64Array(N * N), v = new Float64Array(N * N);
  for (let k = 0; k < N * N; k++) { u[k] = us + 0.02 * (z[k] - 0.5); v[k] = vs + 0.02 * (z[N * N + k] - 0.5); }
  const f = () => new Float64Array(N * N), g = () => new Float64Array(N);
  return { N, dx, dt, d, a, b, t: 0, steps: 0, u, v, P, inv1, invd, wr: f(), wi: f(), br: g(), bi: g(), cr: g(), ci: g() };
}

/** One IMEX Euler step. u and v share one complex FFT: w = u + i v. */
export function stepTuring(S: Turing): void {
  const { N, dt, a, b, u, v, wr, wi, inv1, invd } = S;
  const NN = N * N;
  for (let k = 0; k < NN; k++) {
    const uu = u[k], vv = v[k], q = uu * uu * vv;
    wr[k] = uu + dt * (a - uu + q);
    wi[k] = vv + dt * (b - q);
  }
  fft2(S.P, wr, wi, -1, S.br, S.bi, S.cr, S.ci);
  // split W into U (from real u) and V (from real v) using W(-k) symmetry, divide, recombine
  for (let j = 0; j < N; j++) {
    const jm = (N - j) % N;
    for (let i = 0; i < N; i++) {
      const k = i + N * j, km = (N - i) % N + N * jm;
      if (km < k) continue;
      const ar = wr[k], ai = wi[k], cr = wr[km], ci = wi[km];
      // U_k = (W_k + conj W_-k)/2, V_k = (W_k - conj W_-k)/(2i)
      const Ur = 0.5 * (ar + cr), Ui = 0.5 * (ai - ci);
      const Vr = 0.5 * (ai + ci), Vi = -0.5 * (ar - cr);
      const p1 = inv1[k], pd = invd[k];     // symmetric: same at -k
      const ur = Ur * p1, ui = Ui * p1, vr = Vr * pd, vi = Vi * pd;
      // Z_k = U'_k + i V'_k ; Z_-k = conj(U'_k) + i conj(V'_k)
      wr[k] = ur - vi; wi[k] = ui + vr;
      wr[km] = ur + vi; wi[km] = -ui + vr;
    }
  }
  fft2(S.P, wr, wi, 1, S.br, S.bi, S.cr, S.ci);
  const sc = 1 / NN;
  for (let k = 0; k < NN; k++) { u[k] = wr[k] * sc; v[k] = wi[k] * sc; }
  S.steps++;
  S.t = S.steps * dt;
}

/** Local maxima above the mean among the 8 periodic neighbours (the unit's spots()). */
export function countSpots(u: Float64Array, N: number): number {
  let m = 0;
  for (let k = 0; k < N * N; k++) m += u[k];
  m /= N * N;
  let c = 0;
  for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
    const x = u[i + N * j];
    if (!(x > m)) continue;
    let ok = true;
    for (let dj = -1; dj <= 1 && ok; dj++) for (let di = -1; di <= 1; di++) {
      if (!di && !dj) continue;
      if (!(x > u[(i + di + N) % N + N * ((j + dj + N) % N)])) { ok = false; break; }
    }
    if (ok) c++;
  }
  return c;
}

/** Reaction Jacobian at the uniform state. */
export function jac(a: number, b: number): [number, number, number, number] {
  const us = a + b, vs = b / (us * us);
  return [-1 + 2 * us * vs, us * us, -2 * us * vs, -us * us];   // J11 J12 J21 J22
}

/** Largest real part of the eigenvalues of [[p, q], [r, s]]. */
export function maxRe(p: number, q: number, r: number, s: number): number {
  const tr = p + s, det = p * s - q * r, disc = tr * tr / 4 - det;
  return disc >= 0 ? tr / 2 + Math.sqrt(disc) : tr / 2;
}

/** Continuous growth rate of wavenumber q. */
export function growth(q: number, d: number, a: number, b: number): number {
  const [p, qq, r, s] = jac(a, b), q2 = q * q;
  return maxRe(p - q2, qq, r, s - d * q2);
}

/** Spectral radius of the 2x2 matrix M. */
function specRad(p: number, q: number, r: number, s: number): number {
  const tr = p + s, det = p * s - q * r, disc = tr * tr / 4 - det;
  if (disc >= 0) { const sq = Math.sqrt(disc); return Math.max(Math.abs(tr / 2 + sq), Math.abs(tr / 2 - sq)); }
  return Math.sqrt(Math.abs(det));
}

/** Growth rate per unit time of mode q under the grid's IMEX step: ln(rho(M))/dt,
 *  M = diag(1/(1+dt k2), 1/(1+dt d k2)) (I + dt J), k2 = (4/dx^2) sin^2(q dx/2). */
export function growthIMEX(q: number, d: number, a: number, b: number, dx: number, dt: number): number {
  const [p, qq, r, s] = jac(a, b);
  const k2 = 4 / (dx * dx) * Math.sin(q * dx / 2) ** 2;
  const g1 = 1 / (1 + dt * k2), g2 = 1 / (1 + dt * d * k2);
  return Math.log(specRad(g1 * (1 + dt * p), g1 * dt * qq, g2 * dt * r, g2 * (1 + dt * s))) / dt;
}

// ------------------------------------------------------------------ 1D morphogen gradient
export interface Grad { n: number; dx: number; L: number; c: Float64Array; t: number; lo: Float64Array; di: Float64Array; up: Float64Array; cp: Float64Array; dp: Float64Array }

export function makeGrad(L: number, n: number): Grad {
  const f = () => new Float64Array(n + 1);
  return { n, dx: L / n, L, c: f(), t: 0, lo: f(), di: f(), up: f(), cp: f(), dp: f() };
}

/** Solve tridiagonal (lo, di, up) x = rhs in place of rhs (Thomas). */
function thomas(lo: Float64Array, di: Float64Array, up: Float64Array, rhs: Float64Array, cp: Float64Array): void {
  const n = di.length;
  cp[0] = up[0] / di[0]; rhs[0] /= di[0];
  for (let i = 1; i < n; i++) {
    const m = di[i] - lo[i] * cp[i - 1];
    cp[i] = up[i] / m;
    rhs[i] = (rhs[i] - lo[i] * rhs[i - 1]) / m;
  }
  for (let i = n - 2; i >= 0; i--) rhs[i] -= cp[i] * rhs[i + 1];
}

/** Fill lo, di, up with A = I + dt K_D (ghost points: influx at 0, wall at L). */
function setMatrix(G: Grad, D: number, dt: number): void {
  const n = G.n, r = D / (G.dx * G.dx);
  for (let i = 0; i <= n; i++) { G.lo[i] = -dt * r; G.di[i] = 1 + 2 * dt * r; G.up[i] = -dt * r; }
  G.lo[0] = 0; G.up[0] = -2 * dt * r; G.lo[n] = -2 * dt * r; G.up[n] = 0;
}

/** nsteps IMEX Euler steps: (I + dt K_D) c+ = c + dt (-k(x) c + src). */
export function stepGrad(G: Grad, D: number, kx: Float64Array, j: number, dt: number, nsteps: number): void {
  setMatrix(G, D, dt);
  const c = G.c, n = G.n;
  for (let s = 0; s < nsteps; s++) {
    for (let i = 0; i <= n; i++) c[i] = c[i] - dt * kx[i] * c[i];
    c[0] += dt * 2 * j / G.dx;
    thomas(G.lo, G.di, G.up, c, G.cp);
    G.t += dt;
  }
}

/** Exact steady state with k = k1 on [0, m) and k2 on [m, L] (m = L/2 or k1 = k2). */
export function steadyExact(D: number, k1: number, k2: number, j: number, L: number): (x: number) => number {
  const l1 = Math.sqrt(D / k1), l2 = Math.sqrt(D / k2), m = L / 2;
  if (k1 === k2) return x => j * l1 / D * Math.cosh((L - x) / l1) / Math.sinh(L / l1);
  // c1 = P cosh(x/l1) + Q sinh(x/l1), Q = -j l1 / D ; c2 = B cosh((L - x)/l2)
  const Q = -j * l1 / D;
  const ch = Math.cosh(m / l1), sh = Math.sinh(m / l1), C2 = Math.cosh((L - m) / l2), S2 = Math.sinh((L - m) / l2);
  // P ch - B C2 = -Q sh ;  P sh / l1 + B S2 / l2 = -Q ch / l1
  const a11 = ch, a12 = -C2, a21 = sh / l1, a22 = S2 / l2, r1 = -Q * sh, r2 = -Q * ch / l1;
  const det = a11 * a22 - a12 * a21;
  const P = (r1 * a22 - a12 * r2) / det, B = (a11 * r2 - a21 * r1) / det;
  return x => x < m ? P * Math.cosh(x / l1) + Q * Math.sinh(x / l1) : B * Math.cosh((L - x) / l2);
}

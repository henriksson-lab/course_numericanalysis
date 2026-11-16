// Flow plus diffusion explorer for the advection-diffusion unit. Three tabs:
//  "layer"   the steady back-diffusion test problem of the lecture (L = 500 µm,
//            u = 100 µm/s, D = 1000 µm²/s, c(0) = 0, c(L) = 1). A slider sets the grid,
//            i.e. the cell Péclet number P = u Δx / D. (a) central differences against the
//            exact solution: wiggles when P > 2; (b) upwind differences, the exact solution
//            and the exact solution with D + u Δx / 2 (numerical diffusion).
//  "pulse"   the stimulus pulse of the lecture (u = 200 µm/s, D = 100 µm²/s, 1 s pulse,
//            channel 3000 µm, cells at 2000 µm), explicit Euler in time with upwind and
//            with central differences. Sliders for Δx and the Courant number C = u Δt / Δx.
//            (a) profile at a chosen time; (b) what the cells see. Exact solution
//            (Ogata-Banks, semi-infinite channel) and the exact solution with D + D_num.
//  "tsensor" the T-sensor / H-filter seen from above: w = 100 µm, Poiseuille flow,
//            sample in y < 50 µm, implicit marching in x (backward Euler, cell-centred grid
//            of M = 50 cells, no-flux walls), as in the lecture. Heat map, 10-90 % contour
//            lines of the interdiffusion zone, the mixing length (all |c - 1/2| < 0.01),
//            the outlet profile and the H-filter split for three species.
// Auto-mounts on <div data-demo="advect">. Optional attributes:
//   data-tabs="layer pulse tsensor"  (which tabs, in this order)  data-tab="layer"
//   data-p="5" (starting cell Péclet number of the layer tab, nearest grid)
//   data-dx="10" data-c="0.5" data-t="6" (pulse tab)
//   data-d="1000" data-ubar="1000" data-l="10" (T-sensor: µm²/s, µm/s, mm)
import { svgEl, htmlEl, slider, presets, button, checkbox, scale, autoMount } from "../lib/ui";

const NAME = "advect";
// Okabe-Ito, as in make-advection-diffusion-figs.R (central = red, upwind = blue)
const COL = {
  cen: "#D55E00", up: "#0072B2", ink: "#000000", grey: "#999999", orange: "#E69F00",
  green: "#009E73", sky: "#56B4E9", purple: "#CC79A7", warn: "#b00020",
};
// hcl.colors(9, "Blues 3", rev = TRUE), the palette of fig01
const BLUES = ["#F9F9F9", "#E1EEFF", "#C3DBFD", "#A1C4F1", "#79ABE2", "#468FD0", "#0072B4", "#00538E", "#00366C"];

const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .35em; }
[data-demo="${NAME}"] .ad-row { display: flex; flex-wrap: wrap; gap: .3em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .ad-group { display: flex; gap: .35em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .ad-row .qb-slider { flex: 1 1 330px; min-width: 260px; max-width: 520px; grid-template-columns: auto 1fr 7.4em; }
[data-demo="${NAME}"] .ad-plots { display: flex; flex-wrap: wrap; gap: 8px; width: 100%; }
[data-demo="${NAME}"] .ad-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] svg .wt { font-size: 18px; fill: #fff; font-weight: 600; }
[data-demo="${NAME}"] .ad-read { width: 100%; box-sizing: border-box; display: flex; flex-wrap: wrap; gap: .1em 1.6em; }
[data-demo="${NAME}"] .ad-read b { font-weight: 600; }
[data-demo="${NAME}"] .ad-warn { color: ${COL.warn}; font-weight: 600; }
[data-demo="${NAME}"] .ad-ok { color: ${COL.green}; font-weight: 600; }
[data-demo="${NAME}"] .ad-tab { display: flex; flex-direction: column; gap: .35em; width: 100%; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; padding: .2em .6em; }
.reveal [data-demo="${NAME}"].qb-demo .qb-btn { font-size: .9em; padding: .12em .5em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---------------------------------------------------------------- numerics

/** Thomas algorithm, same argument order as the course: sub-, main, super-diagonal, rhs. */
function thomas(a: number[], b: number[], c: number[], d: number[]): number[] {
  const n = d.length, bb = b.slice(), dd = d.slice();
  for (let i = 1; i < n; i++) { const w = a[i] / bb[i - 1]; bb[i] -= w * c[i - 1]; dd[i] -= w * dd[i - 1]; }
  const x = new Array<number>(n);
  x[n - 1] = dd[n - 1] / bb[n - 1];
  for (let i = n - 2; i >= 0; i--) x[i] = (dd[i] - c[i] * x[i + 1]) / bb[i];
  return x;
}

/** exp(z^2) erfc(z) for z >= 0 (Numerical Recipes erfcc, relative error < 1.2e-7). */
function erfcx(z: number): number {
  const t = 1 / (1 + 0.5 * z);
  return t * Math.exp(-1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
    t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
    t * (-0.82215223 + t * 0.17087277)))))))));
}
function erfc(z: number): number {
  return z >= 0 ? Math.exp(-z * z) * erfcx(z) : 2 - Math.exp(-z * z) * erfcx(-z);
}

// ---- layer tab: steady back-diffusion ----
const LY = { L: 500, u: 100, D: 1000 };
const LAYER_N = [400, 200, 100, 50, 40, 25, 20, 16, 10, 8, 5];       // P = 50 / N
const exact1 = (x: number, D = LY.D): number => Math.expm1(LY.u * x / D) / Math.expm1(LY.u * LY.L / D);

function steady(N: number, scheme: "central" | "upwind"): number[] {
  const dx = LY.L / N, P = LY.u * dx / LY.D, m = N - 1;
  const [lo, di, up] = scheme === "central" ? [1 + P / 2, -2, 1 - P / 2] : [1 + P, -2 - P, 1];
  const a = new Array(m).fill(lo), b = new Array(m).fill(di), c = new Array(m).fill(up), d = new Array(m).fill(0);
  d[m - 1] = -up;                                                     // c(L) = 1
  return [0, ...thomas(a, b, c, d), 1];
}

// ---- pulse tab: explicit Euler, upwind and central ----
const PU = { u: 200, D: 100, Lc: 3000, xobs: 2000, tp: 1, tend: 16 };
const PULSE_DX = [1, 2, 2.5, 5, 10, 20];
const BLOW = 1e12;

function stepExact(x: number, t: number, D: number): number {      // c(0, t) = 1, semi-infinite
  if (t <= 0) return 0;
  const s = 2 * Math.sqrt(D * t), u = PU.u;
  const a = (x - u * t) / s, b = (x + u * t) / s;
  return 0.5 * erfc(a) + 0.5 * Math.exp(-a * a) * erfcx(b);       // exp(ux/D) erfc(b) = exp(-a^2) erfcx(b)
}
const pulseExact = (x: number, t: number, D: number): number => stepExact(x, t, D) - stepExact(x, t - PU.tp, D);

interface PulseRun {
  t: number[]; obs: number[];            // time series at the cells
  snaps: Float64Array[];                 // profile at snapshot times SNAP_T
  blownAt: number;                       // time where |c| exceeded BLOW (Infinity if never)
  maxAbs: number[];                      // max |c| at each snapshot
}
const SNAP_DT = 0.1;
const SNAP_T = Array.from({ length: Math.round(PU.tend / SNAP_DT) + 1 }, (_, i) => i * SNAP_DT);

/** Same update as up_pulse() in the figure script (and its central-difference twin). */
function runPulse(scheme: "upwind" | "central", dx: number, C: number): PulseRun {
  const n = Math.round(PU.Lc / dx) + 1, dt = C * dx / PU.u, r = PU.D * dt / (dx * dx);
  const ns = Math.round(PU.tend / dt), iobs = Math.round(PU.xobs / dx);
  let c = new Float64Array(n), old = new Float64Array(n);
  const t: number[] = [], obs: number[] = [], snaps: Float64Array[] = [], maxAbs: number[] = [];
  let si = 0, blownAt = Infinity;
  const snapStep = SNAP_T.map(ts => Math.round(ts / dt));
  const takeSnaps = (k: number) => {
    while (si < SNAP_T.length && snapStep[si] <= k) {
      snaps.push(c.slice());
      let m = 0; for (let i = 0; i < n; i++) m = Math.max(m, Math.abs(c[i]));
      maxAbs.push(m); si++;
    }
  };
  takeSnaps(0);
  const half = C / 2;
  for (let k = 1; k <= ns; k++) {
    [old, c] = [c, old];
    if (scheme === "upwind") {
      for (let i = 1; i < n - 1; i++) c[i] = old[i] - C * (old[i] - old[i - 1]) + r * (old[i - 1] - 2 * old[i] + old[i + 1]);
    } else {
      for (let i = 1; i < n - 1; i++) c[i] = old[i] - half * (old[i + 1] - old[i - 1]) + r * (old[i - 1] - 2 * old[i] + old[i + 1]);
    }
    c[n - 1] = c[n - 2];                                              // outflow: zero slope
    c[0] = k * dt <= PU.tp ? 1 : 0;                                    // inlet: a 1 s pulse
    t.push(k * dt); obs.push(c[iobs]);
    takeSnaps(k);
    if (k % 20 === 0 || k === ns) {
      let m = 0; for (let i = 0; i < n; i++) { const v = Math.abs(c[i]); if (!(v <= m)) m = v; }
      if (!(m < BLOW)) { blownAt = k * dt; break; }
    }
  }
  while (snaps.length < SNAP_T.length) { snaps.push(c.slice().fill(NaN)); maxAbs.push(Infinity); }
  return { t, obs, snaps, blownAt, maxAbs };
}

// ---- T-sensor tab ----
const TS = { w: 100, M: 50 };
const TS_DY = TS.w / TS.M;
const TS_Y = Array.from({ length: TS.M }, (_, j) => (j + 0.5) * TS_DY);
const pois = (y: number, ubar: number): number => 6 * ubar * (y / TS.w) * (1 - y / TS.w);
const SPECIES = [
  { name: "small molecule", D: 1000, col: COL.green },
  { name: "protein", D: 100, col: COL.up },
  { name: "cell", D: 0.05, col: COL.orange },
];

/** One backward-Euler marching step of length dx (µm). */
function marchStep(c: number[], u: number[], D: number, dx: number): number[] {
  const M = c.length, a = u.map(v => D * dx / (v * TS_DY * TS_DY));
  const di = a.map(v => 1 + 2 * v); di[0] = 1 + a[0]; di[M - 1] = 1 + a[M - 1];
  const neg = a.map(v => -v);
  return thomas(neg, di, neg, c);
}
const c0 = (): number[] => TS_Y.map(y => (y < TS.w / 2 ? 1 : 0));
/** Fraction of the solute flux in the buffer half (y > w/2). */
function crossFrac(c: number[], u: number[]): number {
  let top = 0, all = 0;
  const s0 = c0();
  for (let j = 0; j < c.length; j++) { all += u[j] * s0[j]; if (TS_Y[j] > TS.w / 2) top += u[j] * c[j]; }
  return top / all;
}
const mixed = (c: number[]): boolean => c.every(v => Math.abs(v - 0.5) <= 0.01);
/** Where the 10 % or 90 % level crosses, by linear interpolation in y (c decreases with y). */
function levelY(c: number[], lev: number): number {
  if (c[0] < lev) return 0;
  for (let j = 0; j < c.length - 1; j++) {
    if (c[j] >= lev && c[j + 1] < lev) return TS_Y[j] + (c[j] - lev) / (c[j] - c[j + 1]) * TS_DY;
  }
  return TS.w;
}

// ---------------------------------------------------------------- formatting

function num(v: number, d = 3): string {
  if (!Number.isFinite(v)) return v > 0 ? "∞" : v < 0 ? "−∞" : "—";
  if (v === 0) return "0";
  const a = Math.abs(v);
  let s: string;
  if (a >= 1e5 || a < 1e-3) { const [m, e] = v.toExponential(d - 1).split("e"); s = `${m}e${Number(e)}`; }
  else if (a >= 10 ** d) s = v.toFixed(0);
  else s = String(Number(v.toPrecision(d)));
  return s.replace("-", "−");
}
const fix = (v: number, d: number): string => v.toFixed(d).replace("-", "−");
function lenStr(um: number): string {
  if (um >= 1e6) return `${num(um / 1e6, 3)} m`;
  if (um >= 1000) return `${num(um / 1000, 3)} mm`;
  return `${num(um, 3)} µm`;
}

function frame(svg: SVGSVGElement, x: (v: number) => number, y: (v: number) => number,
  xt: number[], yt: number[], x0: number, x1: number, y0: number, y1: number,
  xlab: string, ylab: string, H: number, xfmt = (v: number) => String(v), yfmt = (v: number) => String(v)): void {
  for (const v of yt) {
    svgEl("line", { x1: x(x0), x2: x(x1), y1: y(v), y2: y(v), class: v === 0 ? "axis" : "grid" }, svg);
    svgEl("text", { x: x(x0) - 6, y: y(v) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = yfmt(v);
  }
  for (const v of xt) {
    svgEl("line", { x1: x(v), x2: x(v), y1: y(y0), y2: y(y0) + 5, class: "axis" }, svg);
    svgEl("text", { x: x(v), y: y(y0) + 23, "text-anchor": "middle", class: "sm" }, svg).textContent = xfmt(v);
  }
  svgEl("line", { x1: x(x0), x2: x(x1), y1: y(y0), y2: y(y0), class: "axis" }, svg);
  svgEl("line", { x1: x(x0), x2: x(x0), y1: y(y0), y2: y(y1), class: "axis" }, svg);
  svgEl("text", { x: (x(x0) + x(x1)) / 2, y: H - 6, "text-anchor": "middle" }, svg).textContent = xlab;
  const ym = (y(y0) + y(y1)) / 2;
  svgEl("text", { x: 16, y: ym, "text-anchor": "middle", transform: `rotate(-90 16 ${ym})` }, svg).textContent = ylab;
}

/** Polyline with values clamped to [lo, hi] (blown-up values show as a band at the edge). */
function poly(svg: SVGSVGElement, xs: ArrayLike<number>, ys: ArrayLike<number>, x: (v: number) => number,
  y: (v: number) => number, lo: number, hi: number, attrs: Record<string, string | number>): void {
  const pts: string[] = [];
  for (let i = 0; i < xs.length; i++) {
    const v = ys[i];
    if (!Number.isFinite(v)) { continue; }
    pts.push(`${x(xs[i]).toFixed(1)},${y(Math.max(lo, Math.min(hi, v))).toFixed(1)}`);
  }
  if (pts.length > 1) svgEl("polyline", { points: pts.join(" "), fill: "none", ...attrs }, svg);
}

function readItem(ro: HTMLElement, html: string, cls = ""): void {
  const s = htmlEl("span", cls, "", ro);
  s.innerHTML = html;
}

function nearestIdx(list: number[], v: number): number {
  let best = 0;
  list.forEach((x, i) => { if (Math.abs(Math.log(x / v)) < Math.abs(Math.log(list[best] / v))) best = i; });
  return best;
}

/** Coalesce redraws into one per animation frame (the pulse runs can take a moment). */
function scheduler(fn: () => void): () => void {
  let pending = false;
  return () => {
    if (pending) return;
    pending = true;
    requestAnimationFrame(() => { pending = false; fn(); });
  };
}

// ---------------------------------------------------------------- layer tab

function buildLayer(box: HTMLElement, ds: DOMStringMap): void {
  const P0 = Number(ds.p) > 0 ? Number(ds.p) : 5;
  const start = nearestIdx(LAYER_N.map(N => LY.u * LY.L / N / LY.D), P0);
  const st = { k: start, whole: false };
  const Pof = (k: number) => LY.u * (LY.L / LAYER_N[k]) / LY.D;

  const ctr = htmlEl("div", "ad-row", "", box);
  const pS = slider(ctr, "cell Péclet P = uΔx/D", 0, LAYER_N.length - 1, 1, st.k,
    v => { st.k = v; draw(); }, v => `${num(Pof(v))} (Δx ${num(LY.L / LAYER_N[v])})`);
  const g = htmlEl("div", "ad-group", "", ctr);
  htmlEl("span", "", "lecture grids:", g);
  const lect = [50, 25, 10];                      // N = 50, 25, 10: P = 1, 2, 5
  const pBtns = presets(g, lect.map(N => ({
    label: `P = ${num(LY.u * LY.L / N / LY.D)}`, apply: () => { st.k = LAYER_N.indexOf(N); pS.set(st.k); draw(); },
  })));
  const ctr2 = htmlEl("div", "ad-row", "", box);
  const wc = checkbox(ctr2, "whole channel (0–500 µm)", st.whole, v => { st.whole = v; draw(); });
  button(ctr2, "reset", () => { st.k = start; st.whole = false; wc.checked = false; pS.set(st.k); draw(); });

  const plots = htmlEl("div", "ad-plots", "", box);
  const W = 560, H = 330;
  const svgs = [0, 1].map(i => {
    const s = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": i ? "upwind" : "central" }, plots);
    s.style.width = "calc(50% - 4px)";
    return s;
  });
  const ro = htmlEl("div", "qb-readout ad-read", "", box);

  const Y0 = -0.6, Y1 = 1.1, top = 34, bot = H - 50, left = 62, right = W - 24;
  const xx = Array.from({ length: 501 }, (_, i) => i);

  function draw(): void {
    const N = LAYER_N[st.k], dx = LY.L / N, P = Pof(st.k);
    pBtns.forEach((b, i) => b.classList.toggle("active", lect[i] === N));
    const xg = Array.from({ length: N + 1 }, (_, i) => i * dx);
    const cc = steady(N, "central"), cu = steady(N, "upwind");
    const ex = xg.map(v => exact1(v));
    const Dn = LY.u * dx / 2;
    const x0 = st.whole ? 0 : 250;
    const x = scale(x0, 500, left, right), y = scale(Y0, Y1, bot, top);
    const xt = st.whole ? [0, 100, 200, 300, 400, 500] : [250, 300, 350, 400, 450, 500];
    const inView = (v: number) => v >= x0 - 1e-9;

    svgs.forEach((svg, i) => {
      svg.replaceChildren();
      const sch = i === 0 ? "central" : "upwind";
      svgEl("text", { x: 8, y: 22, class: "ttl" }, svg).textContent =
        `(${i ? "b" : "a"}) ${sch}, Δx = ${num(dx)} µm, P = ${num(P)}`;
      frame(svg, x, y, xt, [-0.5, 0, 0.5, 1], x0, 500, Y0, Y1, "x (µm)", "c", H, String, v => fix(v, 1));
      const xs = xx.filter(inView);
      if (i === 1) poly(svg, xs, xs.map(v => exact1(v, LY.D + Dn)), x, y, Y0, Y1,
        { stroke: COL.grey, "stroke-width": 3, "stroke-dasharray": "8 5" });
      poly(svg, xs, xs.map(v => exact1(v)), x, y, Y0, Y1, { stroke: COL.ink, "stroke-width": 2.5 });
      const cs = i === 0 ? cc : cu, col = i === 0 ? COL.cen : COL.up;
      const idx = xg.map((v, j) => j).filter(j => inView(xg[j]));
      poly(svg, idx.map(j => xg[j]), idx.map(j => cs[j]), x, y, Y0, Y1, { stroke: col, "stroke-width": 2.5 });
      const rad = idx.length > 60 ? 3 : 5.5;
      for (const j of idx) svgEl("circle", { cx: x(xg[j]), cy: y(Math.max(Y0, Math.min(Y1, cs[j]))), r: rad, fill: col }, svg);
      // legend
      const items: [string, string, string][] = i === 0
        ? [["exact", COL.ink, ""], ["central", COL.cen, ""]]
        : [["exact, D", COL.ink, ""], [`exact, D + uΔx/2 = ${num(LY.D + Dn)}`, COL.grey, "8 5"], ["upwind", COL.up, ""]];
      items.forEach(([lab, c, dash], k) => {
        const ly = top + 14 + k * 24, lx = left + 12;
        svgEl("line", { x1: lx, x2: lx + 30, y1: ly, y2: ly, stroke: c, "stroke-width": 3, "stroke-dasharray": dash }, svg);
        svgEl("text", { x: lx + 38, y: ly + 6 }, svg).textContent = lab;
      });
    });

    ro.replaceChildren();
    const errC = Math.max(...cc.map((v, j) => Math.abs(v - ex[j]))), errU = Math.max(...cu.map((v, j) => Math.abs(v - ex[j])));
    const minC = Math.min(...cc);
    const wig = cc.some((v, j) => j > 0 && v < cc[j - 1] - 1e-12);
    const rho = P === 2 ? Infinity : (1 + P / 2) / (1 - P / 2);
    readItem(ro, `<b>Δx</b> = ${num(dx)} µm, ${N} cells`);
    readItem(ro, `<b>central</b>: ρ = ${P === 2 ? "∞" : num(rho)}, min c = ${num(minC, 3)}, max error ${num(errC, 2)}`,
      wig ? "ad-warn" : "");
    readItem(ro, wig ? "wiggles (P > 2)" : "no wiggles", wig ? "ad-warn" : "ad-ok");
    readItem(ro, `<b>upwind</b>: ρ = 1 + P = ${num(1 + P)}, max error ${num(errU, 2)}`);
    readItem(ro, `<b>D<sub>num</sub></b> = uΔx/2 = ${num(Dn)} µm²/s = ${num(Dn / LY.D, 2)} D`);
    readItem(ro, `exact growth per cell e<sup>P</sup> = ${num(Math.exp(P))}`);
  }
  draw();
}

// ---------------------------------------------------------------- pulse tab

function buildPulse(box: HTMLElement, ds: DOMStringMap): void {
  const dx0 = PULSE_DX[nearestIdx(PULSE_DX, Number(ds.dx) > 0 ? Number(ds.dx) : 10)];
  const C0 = Number(ds.c) > 0 ? Math.min(1.2, Math.max(0.1, Number(ds.c))) : 0.5;
  const t0 = Number.isFinite(Number(ds.t)) && ds.t !== undefined ? Math.min(16, Math.max(0, Number(ds.t))) : 6;
  const cen0 = ds.central === "1";
  const st = { dx: dx0, C: C0, t: t0, central: cen0 };

  const ctr = htmlEl("div", "ad-row", "", box);
  const dxS = slider(ctr, "Δx (µm)", 0, PULSE_DX.length - 1, 1, PULSE_DX.indexOf(st.dx),
    v => { st.dx = PULSE_DX[v]; rerun(); }, v => `${PULSE_DX[v]} (P = ${num(PU.u * PULSE_DX[v] / PU.D)})`);
  const cS = slider(ctr, "Courant C = uΔt/Δx", 0.1, 1.2, 0.05, st.C, v => { st.C = v; rerun(); }, v => fix(v, 2));
  const ctr2 = htmlEl("div", "ad-row", "", box);
  const tS = slider(ctr2, "snapshot time t (s)", 0, 16, 0.1, st.t, v => { st.t = v; schedule(); }, v => fix(v, 1));
  const g = htmlEl("div", "ad-group", "", ctr2);
  const PRE = [
    { label: "lecture: Δx 10, C 0.5", dx: 10, C: 0.5, t: 10.5 },
    { label: "C 0.9", dx: 10, C: 0.9, t: 6 },
    { label: "C 0.95", dx: 10, C: 0.95, t: 6 },
    { label: "Δx 1", dx: 1, C: 0.5, t: 10.5 },
  ];
  const pBtns = presets(g, PRE.map(p => ({
    label: p.label, apply: () => {
      st.dx = p.dx; st.C = p.C; st.t = p.t;
      dxS.set(PULSE_DX.indexOf(p.dx)); cS.set(p.C); tS.set(p.t); rerun();
    },
  })));
  const cb = checkbox(g, "central", st.central, v => { st.central = v; schedule(); });
  button(g, "reset", () => {
    st.dx = dx0; st.C = C0; st.t = t0; st.central = cen0; cb.checked = cen0;
    dxS.set(PULSE_DX.indexOf(st.dx)); cS.set(st.C); tS.set(st.t); rerun();
  });

  const plots = htmlEl("div", "ad-plots", "", box);
  const H = 320, WA = 620, WB = 500, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "profile along the channel" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 4px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "concentration at the cells against time" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 4px)`;
  const legRow = htmlEl("div", "ad-row", "", box);
  const ro = htmlEl("div", "qb-readout ad-read", "", box);

  let runs: { up: PulseRun; cen: PulseRun; key: string } | null = null;
  // exact time series at the cells (fixed)
  const TT = Array.from({ length: 801 }, (_, i) => i * 0.02);
  const exObs = TT.map(t => pulseExact(PU.xobs, t, PU.D));
  const exPeak = Math.max(...exObs), exTpk = TT[exObs.indexOf(exPeak)];

  function rerun(): void { schedule(); }
  const schedule = scheduler(draw);

  const Y0 = -0.5, Y1 = 1.4, top = 34, bot = H - 50;

  function peak(r: PulseRun): [number, number] {
    let m = -Infinity, tm = 0;
    r.obs.forEach((v, i) => { if (v > m) { m = v; tm = r.t[i]; } });
    return [m, tm];
  }

  function draw(): void {
    const { dx, C } = st;
    const key = `${dx}|${C}`;
    if (!runs || runs.key !== key) runs = { up: runPulse("upwind", dx, C), cen: runPulse("central", dx, C), key };
    pBtns.forEach((b, i) => b.classList.toggle("active", PRE[i].dx === dx && Math.abs(PRE[i].C - C) < 1e-9 && Math.abs(PRE[i].t - st.t) < 1e-9));
    const dt = C * dx / PU.u, r = PU.D * dt / (dx * dx), P = PU.u * dx / PU.D;
    const Dn = PU.u * dx * (1 - C) / 2, Deff = PU.D + Dn;
    const si = Math.round(st.t / SNAP_DT);
    const n = Math.round(PU.Lc / dx) + 1;
    const xg = Array.from({ length: n }, (_, i) => i * dx);

    // ---------- (a) profile at time t ----------
    svgA.replaceChildren();
    const ax = scale(0, PU.Lc, 62, WA - 28), ay = scale(Y0, Y1, bot, top);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgA).textContent = `(a) profile at t = ${fix(st.t, 1)} s`;
    frame(svgA, ax, ay, [0, 500, 1000, 1500, 2000, 2500, 3000], [0, 0.5, 1], 0, PU.Lc, Y0, Y1, "x (µm)", "c", H, String, v => fix(v, 1));
    svgEl("line", { x1: ax(PU.xobs), x2: ax(PU.xobs), y1: ay(Y0), y2: ay(Y1), stroke: COL.purple, "stroke-width": 2, "stroke-dasharray": "3 4" }, svgA);
    svgEl("text", { x: ax(PU.xobs) + 6, y: ay(Y1) + 16, fill: COL.purple }, svgA).textContent = "cells";
    const xe = Array.from({ length: 601 }, (_, i) => i * 5);
    if (Deff > 0) poly(svgA, xe, xe.map(v => pulseExact(v, st.t, Deff)), ax, ay, Y0, Y1, { stroke: COL.grey, "stroke-width": 3, "stroke-dasharray": "8 5" });
    poly(svgA, xe, xe.map(v => pulseExact(v, st.t, PU.D)), ax, ay, Y0, Y1, { stroke: COL.ink, "stroke-width": 2.5 });
    if (st.central) poly(svgA, xg, runs.cen.snaps[si], ax, ay, Y0, Y1, { stroke: COL.cen, "stroke-width": 2 });
    poly(svgA, xg, runs.up.snaps[si], ax, ay, Y0, Y1, { stroke: COL.up, "stroke-width": 2.5 });
    legRow.replaceChildren();
    const sw = (col: string, dash: boolean, lab: string) =>
      readItem(legRow, `<span style="display:inline-block;width:30px;border-top:3px ${dash ? "dashed" : "solid"} ${col};vertical-align:middle;margin-right:.35em"></span>${lab}`);
    sw(COL.ink, false, "exact, D = 100");
    if (Deff > 0) sw(COL.grey, true, `exact, D + D<sub>num</sub> = ${num(Deff)}`);
    sw(COL.up, false, "upwind");
    if (st.central) sw(COL.cen, false, "central");
    const blownNote = (run: PulseRun, name: string) => {
      const m = run.maxAbs[si];
      if (m > 1.5) readItem(legRow, `${name} at t = ${fix(st.t, 1)} s: max |c| = ${Number.isFinite(m) ? num(m, 2) : "> 1e12 (stopped)"}`, "ad-warn");
    };
    blownNote(runs.up, "upwind");
    if (st.central) blownNote(runs.cen, "central");

    // ---------- (b) at the cells ----------
    svgB.replaceChildren();
    const bx = scale(0, 16, 62, WB - 20), by = scale(Y0, Y1, bot, top);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgB).textContent = "(b) what the cells see (x = 2000 µm)";
    frame(svgB, bx, by, [0, 4, 8, 12, 16], [0, 0.5, 1], 0, 16, Y0, Y1, "time (s)", "c at the cells", H, String, v => fix(v, 1));
    svgEl("line", { x1: bx(st.t), x2: bx(st.t), y1: by(Y0), y2: by(Y1), stroke: "#777", "stroke-width": 1.5, "stroke-dasharray": "3 4" }, svgB);
    if (Deff > 0) poly(svgB, TT, TT.map(t => pulseExact(PU.xobs, t, Deff)), bx, by, Y0, Y1, { stroke: COL.grey, "stroke-width": 3, "stroke-dasharray": "8 5" });
    poly(svgB, TT, exObs, bx, by, Y0, Y1, { stroke: COL.ink, "stroke-width": 2.5 });
    const thin = (rr: PulseRun) => {                  // at most ~1600 points
      const s = Math.max(1, Math.floor(rr.t.length / 1600));
      const ti: number[] = [], oi: number[] = [];
      for (let i = 0; i < rr.t.length; i += s) { ti.push(rr.t[i]); oi.push(rr.obs[i]); }
      return [ti, oi];
    };
    if (st.central) { const [ti, oi] = thin(runs.cen); poly(svgB, ti, oi, bx, by, Y0, Y1, { stroke: COL.cen, "stroke-width": 2 }); }
    { const [ti, oi] = thin(runs.up); poly(svgB, ti, oi, bx, by, Y0, Y1, { stroke: COL.up, "stroke-width": 2.5 }); }

    // ---------- readout ----------
    ro.replaceChildren();
    const upOK = C + 2 * r <= 1 + 1e-12;
    const cenOK = C * C <= 2 * r + 1e-12 && 2 * r <= 1 + 1e-12;
    readItem(ro, `<b>Δt</b> = ${num(dt)} s, <b>r</b> = DΔt/Δx² = ${num(r)}, P = ${num(P)}`);
    readItem(ro, `<b>upwind</b>: C + 2r = ${fix(C + 2 * r, 3)} ${upOK ? "≤ 1, stable" : "> 1, unstable"}`, upOK ? "ad-ok" : "ad-warn");
    readItem(ro, `<b>central</b>: needs C ≤ 2/P = ${num(2 / P)} ${cenOK ? "(met)" : "(violated)"}`, cenOK ? "ad-ok" : "ad-warn");
    readItem(ro, `<b>D<sub>num</sub></b>: upwind uΔx(1−C)/2 = ${num(Dn)}, central −CuΔx/2 = ${num(-C * PU.u * dx / 2)} µm²/s`);
    const [pu, tu] = peak(runs.up), [pc, tc] = peak(runs.cen);
    const pk = (v: number, t: number, bl: number) => Number.isFinite(bl) ? `blew up at t = ${fix(bl, 1)} s`
      : Math.abs(v) > 10 ? `max ${num(v, 2)} (unstable)` : `${fix(v, 4)} at ${fix(t, 2)} s`;
    const rng = (rr: PulseRun) => {
      const sn = rr.snaps[si];
      let lo = Infinity, hi = -Infinity;
      for (const v of sn) { if (v < lo) lo = v; if (v > hi) hi = v; }
      return Number.isFinite(lo) && Number.isFinite(hi) && rr.maxAbs[si] < BLOW ? `${num(lo, 3)} … ${num(hi, 3)}` : "blown up";
    };
    readItem(ro, `<b>range of c at t = ${fix(st.t, 1)} s</b>: upwind ${rng(runs.up)}; central ${rng(runs.cen)}`);
    readItem(ro, `<b>peak at the cells</b>: exact ${fix(exPeak, 3)} at ${fix(exTpk, 2)} s; upwind ${pk(pu, tu, runs.up.blownAt)}; central ${pk(pc, tc, runs.cen.blownAt)}`);
  }
  draw();
}

// ---------------------------------------------------------------- T-sensor tab

function blue(v: number): [number, number, number] {
  const t = Math.max(0, Math.min(1, v)) * (BLUES.length - 1);
  const i = Math.min(BLUES.length - 2, Math.floor(t)), f = t - i;
  const h = (s: string, k: number) => parseInt(s.slice(1 + 2 * k, 3 + 2 * k), 16);
  return [0, 1, 2].map(k => Math.round(h(BLUES[i], k) * (1 - f) + h(BLUES[i + 1], k) * f)) as [number, number, number];
}

function buildTsensor(box: HTMLElement, ds: DOMStringMap): void {
  const D0 = Number(ds.d) > 0 ? Number(ds.d) : 1000;
  const U0 = Number(ds.ubar) > 0 ? Number(ds.ubar) : 1000;
  const L0 = Number(ds.l) > 0 ? Number(ds.l) : 5;
  const st = { D: D0, ubar: U0, Lmm: L0 };

  const ctr = htmlEl("div", "ad-row", "", box);
  const dS = slider(ctr, "D (µm²/s)", -2, 3.5, 0.05, Math.log10(st.D), v => { st.D = 10 ** v; dBtnSync(); schedule(); }, v => num(10 ** v, 2));
  const g = htmlEl("div", "ad-group", "", ctr);
  const dBtns = presets(g, SPECIES.map(s => ({ label: `${s.name} (${s.D})`, apply: () => { st.D = s.D; dS.set(Math.log10(s.D)); schedule(); } })));
  const dBtnSync = () => dBtns.forEach((b, i) => b.classList.toggle("active", Math.abs(SPECIES[i].D - st.D) / st.D < 1e-9));
  const ctr2 = htmlEl("div", "ad-row", "", box);
  const uS = slider(ctr2, "mean speed ū (µm/s)", 100, 5000, 100, st.ubar, v => { st.ubar = v; schedule(); }, v => String(v));
  const lS = slider(ctr2, "length L (mm)", 0.5, 40, 0.5, st.Lmm, v => { st.Lmm = v; schedule(); }, v => fix(v, 1));
  button(ctr2, "reset", () => {
    st.D = D0; st.ubar = U0; st.Lmm = L0;
    dS.set(Math.log10(st.D)); uS.set(st.ubar); lS.set(st.Lmm); dBtnSync(); schedule();
  });
  dBtnSync();

  const plots = htmlEl("div", "ad-plots", "", box);
  const H = 340, WA = 700, WB = 420, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "channel from above, concentration map" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 4px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "outlet profile and H-filter split" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 4px)`;
  const ro = htmlEl("div", "qb-readout ad-read", "", box);
  const canvas = document.createElement("canvas");
  const NSTEP = 800, NCOL = 400;
  canvas.width = NCOL; canvas.height = TS.M;
  const schedule = scheduler(draw);

  function draw(): void {
    const { D, ubar } = st, L = st.Lmm * 1000, dx = L / NSTEP;
    dBtnSync();
    const u = TS_Y.map(y => pois(y, ubar));
    // march the current species, keep NCOL columns, contour levels and the mixing point
    let c = c0();
    const cols: number[][] = [c];
    const y10: number[] = [TS.w / 2], y90: number[] = [TS.w / 2], xs: number[] = [0];
    let xmix = mixed(c) ? 0 : NaN;
    for (let k = 1; k <= NSTEP; k++) {
      c = marchStep(c, u, D, dx);
      if (Number.isNaN(xmix) && mixed(c)) xmix = k * dx;
      if (k % (NSTEP / NCOL) === 0) cols.push(c);
      y10.push(levelY(c, 0.1)); y90.push(levelY(c, 0.9)); xs.push(k * dx);
    }
    const cOut = c, fOut = crossFrac(cOut, u);
    // mixing point beyond the view: keep marching with steps of 1 % of x
    let xm = xmix, cc = cOut, xx = L;
    while (Number.isNaN(xm) && xx < 1e9) {
      const h = Math.max(dx, 0.01 * xx);
      cc = marchStep(cc, u, D, h); xx += h;
      if (mixed(cc)) xm = xx;
    }
    // ---------- heat map ----------
    const ctx = canvas.getContext("2d");
    if (ctx) {
      const img = ctx.createImageData(NCOL, TS.M);
      for (let i = 0; i < NCOL; i++) {
        const col = cols[i + 1];
        for (let j = 0; j < TS.M; j++) {
          const [r, g2, b] = blue(col[j]), p = 4 * ((TS.M - 1 - j) * NCOL + i);
          img.data[p] = r; img.data[p + 1] = g2; img.data[p + 2] = b; img.data[p + 3] = 255;
        }
      }
      ctx.putImageData(img, 0, 0);
    }
    svgA.replaceChildren();
    const left = 92, right = WA - 16, top = 58, bot = H - 96;
    const ax = scale(0, L, left, right), ay = scale(0, TS.w, bot, top);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgA).textContent = "(a) channel from above: sample concentration c(x, y)";
    svgEl("image", { href: canvas.toDataURL(), x: left, y: top, width: right - left, height: bot - top, preserveAspectRatio: "none" }, svgA);
    svgEl("rect", { x: left, y: top, width: right - left, height: bot - top, fill: "none", stroke: "#222", "stroke-width": 2 }, svgA);
    // contours of the interdiffusion zone
    const thinIdx = xs.map((_, i) => i).filter(i => i % 2 === 0);
    for (const ys of [y10, y90]) {
      let seg: number[] = [];
      const flush = () => {
        if (seg.length > 1) poly(svgA, seg.map(i => xs[i]), seg.map(i => ys[i]), ax, ay, 0, TS.w,
          { stroke: COL.orange, "stroke-width": 2.5, "stroke-dasharray": "7 4" });
        seg = [];
      };
      for (const i of thinIdx) { if (ys[i] > 0.5 && ys[i] < TS.w - 0.5) seg.push(i); else flush(); }
      flush();
    }
    // inlets
    svgEl("text", { x: left - 8, y: ay(25) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = "sample";
    svgEl("text", { x: left - 8, y: ay(75) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = "buffer";
    // outlet split
    svgEl("line", { x1: right - 2, x2: right + 14, y1: ay(50), y2: ay(50), stroke: "#222", "stroke-width": 3 }, svgA);
    // mixing length
    if (Number.isFinite(xm) && xm <= L) {
      svgEl("line", { x1: ax(xm), x2: ax(xm), y1: top - 6, y2: bot + 6, stroke: COL.purple, "stroke-width": 3 }, svgA);
      svgEl("text", { x: Math.min(ax(xm), right - 150), y: top - 10, fill: COL.purple, "font-weight": 600 }, svgA).textContent =
        `mixed to 1 %: ${lenStr(xm)}`;
    } else {
      svgEl("text", { x: right, y: top - 10, "text-anchor": "end", fill: COL.purple, "font-weight": 600 }, svgA).textContent =
        `mixed to 1 % only at ${Number.isFinite(xm) ? lenStr(xm) : "> 1 km"}`;
    }
    // x axis in mm
    const step = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20].find(s => st.Lmm / s <= 6) ?? 10;
    for (let v = 0; v <= st.Lmm + 1e-9; v += step) {
      svgEl("line", { x1: ax(v * 1000), x2: ax(v * 1000), y1: bot, y2: bot + 5, class: "axis" }, svgA);
      svgEl("text", { x: ax(v * 1000), y: bot + 23, "text-anchor": "middle", class: "sm" }, svgA).textContent = String(Number(v.toFixed(2)));
    }
    svgEl("text", { x: (left + right) / 2, y: bot + 48, "text-anchor": "middle" }, svgA).textContent = "distance downstream x (mm)";
    for (const v of [0, 50, 100]) svgEl("text", { x: left - 8, y: ay(v) + 6 + (v === 0 ? -2 : v === 100 ? 10 : 0), "text-anchor": "end", class: "sm" }, svgA).textContent = v === 50 ? "" : String(v);
    svgEl("text", { x: 16, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(top + bot) / 2})` }, svgA).textContent = "y (µm)";
    // colour bar + contour legend
    const cbY = H - 26, cbX = left;
    for (let i = 0; i < 40; i++) {
      const [r, g2, b] = blue(i / 39);
      svgEl("rect", { x: cbX + i * 4, y: cbY - 14, width: 4.5, height: 16, fill: `rgb(${r},${g2},${b})` }, svgA);
    }
    svgEl("text", { x: cbX - 6, y: cbY, "text-anchor": "end", class: "sm" }, svgA).textContent = "c = 0";
    svgEl("text", { x: cbX + 166, y: cbY }, svgA).textContent = "1";
    svgEl("line", { x1: cbX + 210, x2: cbX + 245, y1: cbY - 6, y2: cbY - 6, stroke: COL.orange, "stroke-width": 3, "stroke-dasharray": "7 4" }, svgA);
    svgEl("text", { x: cbX + 252, y: cbY }, svgA).textContent = "c = 0.1 and 0.9";

    // ---------- (b) outlet profile and split ----------
    svgB.replaceChildren();
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgB).textContent = "(b) outlet: profile and split";
    const bl = 62, br = 170, bx = scale(0, 1, bl, br), by = ay;
    for (const v of [0, 0.5, 1]) {
      svgEl("line", { x1: bx(v), x2: bx(v), y1: top, y2: bot, class: v === 0 ? "axis" : "grid" }, svgB);
      svgEl("text", { x: bx(v), y: bot + 23, "text-anchor": "middle", class: "sm" }, svgB).textContent = v === 0.5 ? "0.5" : String(v);
    }
    svgEl("line", { x1: bl, x2: br, y1: by(50), y2: by(50), stroke: "#222", "stroke-width": 2, "stroke-dasharray": "5 4" }, svgB);
    svgEl("text", { x: (bl + br) / 2, y: bot + 48, "text-anchor": "middle" }, svgB).textContent = "c at x = L";
    svgEl("text", { x: bl - 6, y: by(75) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = "buf.";
    svgEl("text", { x: bl - 6, y: by(25) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = "smp.";
    const pts = TS_Y.map((y, j) => `${bx(cOut[j]).toFixed(1)},${by(y).toFixed(1)}`).join(" ");
    svgEl("polyline", { points: pts, fill: "none", stroke: COL.up, "stroke-width": 3 }, svgB);
    // fractions crossing for the three species and the current D (only D L / ū matters)
    const fx = 192, fw = WB - fx - 8;
    svgEl("text", { x: fx, y: top - 10, class: "sm" }, svgB).textContent = "into buffer outlet";
    const list = SPECIES.map(s => ({ ...s, f: NaN }));
    const custom = !SPECIES.some(s => Math.abs(s.D - D) / D < 1e-9);
    if (custom) list.push({ name: "this D", D, col: COL.purple, f: NaN });
    for (const s of list) {
      if (Math.abs(s.D - D) / D < 1e-9) { s.f = fOut; continue; }
      let cs = c0(); for (let k = 0; k < NSTEP; k++) cs = marchStep(cs, u, s.D, dx);
      s.f = crossFrac(cs, u);
    }
    const fxs = scale(0, 0.5, fx, fx + fw);
    const rowH = (bot - top) / list.length;
    list.forEach((s, k) => {
      const yy = top + k * rowH + 6;
      svgEl("text", { x: fx, y: yy + 16, class: "sm" }, svgB).textContent = s.name;
      svgEl("rect", { x: fx, y: yy + 24, width: fw, height: 20, fill: "#eee" }, svgB);
      svgEl("rect", { x: fx, y: yy + 24, width: Math.max(1, fxs(Math.min(0.5, s.f)) - fx), height: 20, fill: s.col }, svgB);
      const isCur = Math.abs(s.D - D) / D < 1e-9;
      svgEl("text", { x: fx + fw, y: yy + 16, "text-anchor": "end", "font-weight": isCur ? 700 : 400 }, svgB).textContent = `${fix(100 * s.f, 1)} %`;
    });
    svgEl("text", { x: fx + fw, y: bot + 23, "text-anchor": "end", class: "sm" }, svgB).textContent = "50 % = fully mixed";

    // ---------- readout ----------
    ro.replaceChildren();
    const Pe = ubar * TS.w / D;
    readItem(ro, `<b>Pe</b> = ūw/D = ${num(Pe)}`);
    readItem(ro, `<b>ūw²/D</b> = ${lenStr(ubar * TS.w * TS.w / D)}`);
    readItem(ro, `<b>mixed to 1 %</b> (all |c − ½| ≤ 0.01) at ${Number.isFinite(xm) ? lenStr(xm) : "> 1 km"}`);
    const zw = y10[NSTEP] - y90[NSTEP];
    readItem(ro, `<b>10–90 % zone at outlet</b>: ${zw >= TS.w - 1e-9 ? "whole channel" : num(zw, 2) + " µm wide"}`);
    readItem(ro, `<b>into buffer outlet</b>: ${fix(100 * fOut, 1)} %`);
    readItem(ro, `Dx/ū = ${num(D * L / ubar, 3)} µm²`);
  }
  draw();
}

// ---------------------------------------------------------------- mount

const TABS: Record<string, { label: string; build: (b: HTMLElement, ds: DOMStringMap) => void }> = {
  layer: { label: "Steady layer: central vs upwind", build: buildLayer },
  pulse: { label: "Moving pulse: time steps and CFL", build: buildPulse },
  tsensor: { label: "T-sensor / H-filter", build: buildTsensor },
};

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  let names = (ds.tabs ?? "layer pulse tsensor").trim().split(/[\s,]+/).filter(n => n in TABS);
  if (names.length === 0) names = ["layer", "pulse", "tsensor"];
  const first = ds.tab && names.includes(ds.tab) ? ds.tab : names[0];
  const boxes: Record<string, HTMLElement> = {};
  const built = new Set<string>();
  let tabBtns: HTMLButtonElement[] = [];
  const show = (t: string) => {
    if (!built.has(t)) { built.add(t); TABS[t].build(boxes[t], ds); }
    for (const n of names) boxes[n].style.display = n === t ? "flex" : "none";
    tabBtns.forEach((b, i) => b.classList.toggle("active", names[i] === t));
  };
  if (names.length > 1) {
    const row = htmlEl("div", "ad-row", "", host);
    tabBtns = presets(row, names.map(n => ({ label: TABS[n].label, apply: () => show(n) })));
  }
  for (const n of names) boxes[n] = htmlEl("div", "ad-tab", "", host);
  show(first);
}

autoMount(NAME, mount);

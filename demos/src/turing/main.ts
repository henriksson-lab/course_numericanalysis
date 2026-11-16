// Reaction-diffusion explorer for the reaction-diffusion unit. Two tabs:
//  (a) "Morphogen gradient": the unit's synthesis-diffusion-degradation model in 1D,
//      dc/dt = D c'' - k(x) c on 0 < x < 500 µm, influx j at x = 0, wall at x = L,
//      time in minutes. IMEX Euler on a 1 µm grid (Thomas), dt = 0.5 min, from c = 0.
//      The exact steady state (cosh; piecewise with a temperature step: k doubles in
//      the warm half), lambda = sqrt(D/k) marked, two thresholds read as a "French flag",
//      their positions now / at steady state / for the reference D = 100, k = 0.01, j = 1.
//  (b) "Turing patterns": Schnakenberg u_t = a - u + u^2 v + Lap u, v_t = b - u^2 v + d Lap v,
//      a = 0.1, b = 0.9, periodic 80 x 80 box, IMEX Euler with an FFT solve (sim.ts), noise
//      from the unit's lcg (seed 1 gives the unit's 47 spots at t = 200). Grid (dx = 2, 1, 0.5)
//      and time step (0.4 ... 0.05) can be changed; finished runs are logged and the previous
//      one is shown for comparison. Dispersion relation: exact (solid), the grid's IMEX step
//      (dashed), unstable band shaded.
// Auto-mounts on <div data-demo="turing">. Optional attributes:
//   data-tab="gradient|turing"  data-thresholds="0.3 0.1"  data-d="20"  data-dx="1"
//   data-dt="0.1"  data-seed="1"  data-tend="200"
import { svgEl, htmlEl, slider, presets, checkbox, button, scale, numList, autoMount } from "../lib/ui";
import {
  Turing, makeTuring, stepTuring, countSpots, growth, growthIMEX, jac,
  makeGrad, stepGrad, steadyExact,
} from "./sim";

const NAME = "turing";
// Okabe-Ito, as in make-reaction-diffusion-figs.R; viridis (hcl.colors) for the 2D maps
const COL = { blue: "#0072B2", red: "#D55E00", orange: "#E69F00", sky: "#56B4E9", green: "#009E73", purple: "#CC79A7", grey: "#999999", ink: "#222" };
const VIR = ["#4B0055", "#403173", "#00588B", "#007B97", "#009B95", "#00B785", "#53CC67", "#B2DC3C", "#FDE333"];
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .35em; }
[data-demo="${NAME}"] .tu-row { display: flex; flex-wrap: wrap; gap: .3em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .tu-group { display: flex; gap: .35em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .tu-row .qb-slider { min-width: 250px; flex: 1 1 250px; max-width: 380px; grid-template-columns: auto 1fr 4.2em; }
[data-demo="${NAME}"] .tu-box { display: flex; flex-direction: column; gap: .35em; width: 100%; }
[data-demo="${NAME}"] .tu-plots { display: flex; flex-wrap: wrap; gap: 10px; width: 100%; align-items: flex-start; }
[data-demo="${NAME}"] .tu-col { display: flex; flex-direction: column; gap: .3em; min-width: 0; }
[data-demo="${NAME}"] svg { height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] canvas { display: block; width: 100%; height: auto; image-rendering: pixelated; }
[data-demo="${NAME}"] .tu-cap { font-weight: 600; }
[data-demo="${NAME}"] .tu-read { width: 100%; box-sizing: border-box; display: flex; flex-wrap: wrap; gap: .1em 1.3em; }
[data-demo="${NAME}"] .tu-log { font-family: ui-monospace, monospace; font-size: .85em; white-space: pre; background: var(--qb-panel); border-radius: 6px; padding: .3em .45em; overflow: hidden; }
[data-demo="${NAME}"] .tu-big { font-size: 1.3em; font-weight: 700; }
[data-demo="${NAME}"] .tu-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .tu-ok { color: #00704f; font-weight: 600; }
[data-demo="${NAME}"] .tu-tabs .qb-btn { font-size: .95em; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; gap: .2em; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; padding: .2em .45em; }
.reveal [data-demo="${NAME}"] .qb-btn { font-size: .9em; padding: .1em .45em; }
.reveal [data-demo="${NAME}"] .tu-log { font-size: 15px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

function viridis(t: number): [number, number, number] {
  const x = Math.min(1, Math.max(0, t)) * (VIR.length - 1);
  const i = Math.min(VIR.length - 2, Math.floor(x)), u = x - i;
  const a = VIR[i], b = VIR[i + 1];
  const ch = (s: string, o: number) => parseInt(s.slice(o, o + 2), 16);
  return [0, 1, 2].map(m => Math.round(ch(a, 1 + 2 * m) * (1 - u) + ch(b, 1 + 2 * m) * u)) as [number, number, number];
}
const VLUT: [number, number, number][] = Array.from({ length: 256 }, (_, i) => viridis(i / 255));

const num = (s: string | undefined, d: number, lo: number, hi: number) => {
  const v = Number(s);
  return s !== undefined && s !== "" && Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : d;
};
const f1 = (v: number) => v.toFixed(1);
function fmtSig(v: number, d = 3): string {
  const a = Math.abs(v);
  if (a === 0) return "0";
  if (a >= 1e-3 && a < 1e5) return Number(v.toPrecision(d)).toString();
  return v.toExponential(d - 1);
}
const isShown = (el: HTMLElement) => el.offsetParent !== null;

// ======================================================================= mount
function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const tabs = htmlEl("div", "tu-row tu-tabs", "", host);
  const boxA = htmlEl("div", "tu-box", "", host);
  const boxB = htmlEl("div", "tu-box", "", host);
  const tb = presets(tabs, [
    { label: "(a) Morphogen gradient", apply: () => show("gradient") },
    { label: "(b) Turing patterns", apply: () => show("turing") },
  ]);
  const A = buildGradient(boxA, ds);
  const B = buildTuring(boxB, ds);
  function show(t: string): void {
    boxA.style.display = t === "gradient" ? "flex" : "none";
    boxB.style.display = t === "turing" ? "flex" : "none";
    tb.forEach((b, i) => b.classList.toggle("active", (i === 0) === (t === "gradient")));
    A.setActive(t === "gradient"); B.setActive(t === "turing");
  }
  show(ds.tab === "turing" ? "turing" : "gradient");
}

interface TabCtl { setActive(on: boolean): void }

// ======================================================================= (a) gradient
const L = 500, NG = 500, DT1 = 0.5;          // µm, 1 µm grid, min
const REF = { D: 100, k: 0.01, j: 1 };
type Warm = "none" | "post" | "ant" | "all";

function buildGradient(box: HTMLElement, ds: DOMStringMap): TabCtl {
  const thr = numList(ds.thresholds, [0.3, 0.1]).slice(0, 2).sort((p, q) => q - p);
  while (thr.length < 2) thr.push(0.1);
  const st = { D: REF.D, k: REF.k, j: REF.j, warm: "none" as Warm, log: false, playing: false, active: false, done: false, tDone: NaN };
  const G = makeGrad(L, NG);
  const kx = new Float64Array(NG + 1);
  const xs = Array.from({ length: NG + 1 }, (_, i) => i * L / NG);

  // ---- controls
  const r1 = htmlEl("div", "tu-row", "", box);
  const dS = slider(r1, "D (µm²/min)", 25, 400, 5, st.D, v => { st.D = v; changed(); }, v => v.toFixed(0));
  const kS = slider(r1, "k (1/min)", 0.0025, 0.04, 0.0005, st.k, v => { st.k = v; changed(); }, v => v.toFixed(4));
  const jS = slider(r1, "influx j", 0.25, 4, 0.05, st.j, v => { st.j = v; changed(); }, v => v.toFixed(2));
  const r2 = htmlEl("div", "tu-row", "", box);
  const wg = htmlEl("div", "tu-group", "", r2);
  htmlEl("span", "", "warm half (k × 2):", wg);
  const WARMS: Warm[] = ["none", "post", "ant", "all"];
  const wBtns = presets(wg, [
    { label: "none", apply: () => setWarm("none") },
    { label: "posterior", apply: () => setWarm("post") },
    { label: "anterior", apply: () => setWarm("ant") },
    { label: "whole embryo", apply: () => setWarm("all") },
  ]);
  const pg = htmlEl("div", "tu-group", "", r2);
  const playB = button(pg, "▶ run", () => { if (st.done) restart(); st.playing = !st.playing; kick(); draw(); });
  button(pg, "reset c = 0", () => { restart(); draw(); });
  button(pg, "jump to steady", () => { jumpSteady(); draw(); });
  button(pg, "reference", () => {
    st.D = REF.D; st.k = REF.k; st.j = REF.j; dS.set(st.D); kS.set(st.k); jS.set(st.j); setWarm("none");
  });
  checkbox(r2, "log scale", st.log, v => { st.log = v; draw(); });

  const W = 1120, H = 380;
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Morphogen profile c(x) with steady state and thresholds" }, box);
  svg.style.width = "100%";
  const ro = htmlEl("div", "qb-readout tu-read", "", box);

  const left = 78, right = W - 18, top = 40, bot = H - 124;
  const sx = scale(0, L, left, right);

  function kAt(x: number): number {
    const warmHere = st.warm === "all" || (st.warm === "post" && x >= L / 2) || (st.warm === "ant" && x < L / 2);
    return warmHere ? 2 * st.k : st.k;
  }
  function steadyFn(): (x: number) => number {
    return steadyExact(st.D, kAt(0), kAt(L), st.j, L);
  }
  const refFn = steadyExact(REF.D, REF.k, REF.k, REF.j, L);
  function crossExact(f: (x: number) => number, T: number): number {
    if (f(0) <= T) return -1;
    if (f(L) > T) return Infinity;
    let a = 0, b = L;
    for (let i = 0; i < 60; i++) { const m = (a + b) / 2; if (f(m) > T) a = m; else b = m; }
    return (a + b) / 2;
  }
  function crossGrid(c: Float64Array, T: number): number {
    if (c[0] <= T) return -1;
    for (let i = 1; i <= NG; i++) if (c[i] <= T) return xs[i - 1] + (c[i - 1] - T) / (c[i - 1] - c[i]) * (xs[i] - xs[i - 1]);
    return Infinity;
  }
  const xTxt = (x: number) => x < 0 ? "none (c < threshold)" : !Number.isFinite(x) ? "beyond L" : f1(x) + " µm";
  const refX = thr.map(T => crossExact(refFn, T));
  const tTxt = () => Number.isFinite(G.t) ? `${G.t.toFixed(0)} min` : "∞ (steady state)";

  function setWarm(w: Warm): void { st.warm = w; changed(); }
  function changed(): void {
    for (let i = 0; i <= NG; i++) kx[i] = kAt(xs[i]);
    st.done = false; st.tDone = NaN;
    if (!Number.isFinite(G.t)) G.t = 0;
    if (st.active) { st.playing = true; kick(); }
    draw();
  }
  function restart(): void { G.c.fill(0); G.t = 0; st.done = false; st.tDone = NaN; }
  function jumpSteady(): void { const f = steadyFn(); for (let i = 0; i <= NG; i++) G.c[i] = f(xs[i]); st.playing = false; st.done = true; st.tDone = NaN; G.t = NaN; }
  function nearSteady(): boolean {
    const f = steadyFn();
    for (let i = 0; i <= NG; i += 5) { const s = f(xs[i]); if (Math.abs(G.c[i] - s) > 0.01 * s) return false; }
    return true;
  }

  // ---- animation: about 100 min of model time per second
  let raf = 0, last = 0;
  function kick(): void { if (!raf) { last = performance.now(); raf = requestAnimationFrame(tick); } }
  function tick(now: number): void {
    raf = 0;
    if (!st.playing || !st.active) return;
    const el = Math.min(100, now - last); last = now;
    if (isShown(box)) {
      const n = Math.max(1, Math.round(0.1 * el / DT1));
      stepGrad(G, st.D, kx, st.j, DT1, n);
      if (nearSteady()) { st.done = true; st.playing = false; st.tDone = G.t; }
      else if (G.t >= 5000) st.playing = false;
      draw();
    }
    if (st.playing) raf = requestAnimationFrame(tick);
  }

  function draw(): void {
    wBtns.forEach((b, i) => b.classList.toggle("active", WARMS[i] === st.warm));
    playB.textContent = st.playing ? "❚❚ pause" : "▶ run";
    const f = steadyFn();
    const lam0 = Math.sqrt(st.D / kAt(0));
    const c0 = f(0);
    // y scale
    const tops = [1.2, 2.5, 5, 10, 20];
    const ymax = tops.find(v => v >= 1.08 * Math.max(c0, refFn(0))) ?? 40;
    const y0 = st.log ? Math.log10(1e-4) : 0, y1 = st.log ? Math.log10(ymax) : ymax;
    const syy = scale(y0, y1, bot, top);
    const sy = (c: number) => st.log ? syy(Math.log10(Math.max(c, 10 ** y0))) : syy(Math.min(c, ymax * 1.02));

    svg.replaceChildren();
    // warm region
    if (st.warm !== "none") {
      const xa = st.warm === "post" ? L / 2 : 0, xb = st.warm === "ant" ? L / 2 : L;
      svgEl("rect", { x: sx(xa), y: top, width: sx(xb) - sx(xa), height: bot - top, fill: COL.orange, "fill-opacity": 0.13 }, svg);
      const t = svgEl("text", { x: (sx(xa) + sx(xb)) / 2, y: top + 24, "text-anchor": "middle", "font-weight": 600 }, svg);
      t.style.fill = "#9a6a00"; t.textContent = "warm: k × 2";
    }
    // axes
    const yt = st.log ? [1e-4, 1e-3, 0.01, 0.1, 1, 10].filter(v => v <= ymax) : (() => { const stp = ymax === 1.2 ? 0.2 : ymax / 5; return Array.from({ length: Math.round(ymax / stp) + 1 }, (_, i) => Number((i * stp).toPrecision(3))); })();
    for (const v of yt) {
      svgEl("line", { x1: left, x2: right, y1: sy(v), y2: sy(v), class: "grid" }, svg);
      svgEl("text", { x: left - 6, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = String(v);
    }
    svgEl("line", { x1: left, x2: right, y1: bot, y2: bot, class: "axis" }, svg);
    svgEl("line", { x1: left, x2: left, y1: top, y2: bot, class: "axis" }, svg);
    for (let x = 0; x <= L; x += 50) {
      svgEl("line", { x1: sx(x), x2: sx(x), y1: bot, y2: bot + 5, class: "axis" }, svg);
      svgEl("text", { x: sx(x), y: bot + 22, "text-anchor": "middle", class: "sm" }, svg).textContent = String(x);
    }
    const yl = svgEl("text", { x: 18, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 18 ${(top + bot) / 2})` }, svg);
    yl.textContent = st.log ? "c (log scale)" : "concentration c";
    svgEl("text", { x: left, y: 22, class: "ttl" }, svg).textContent = `profile at t = ${tTxt()}`;
    // thresholds
    thr.forEach((T, i) => {
      svgEl("line", { x1: left, x2: right, y1: sy(T), y2: sy(T), stroke: i === 0 ? COL.blue : COL.red, "stroke-width": 1.6, "stroke-dasharray": "3 4" }, svg);
      const t = svgEl("text", { x: right - 4, y: sy(T) - 6, "text-anchor": "end" }, svg);
      t.style.fill = i === 0 ? COL.blue : COL.red; t.textContent = `threshold ${T}`;
    });
    // lambda marker: x = lambda (anterior region), level c(0)/e
    if (lam0 <= L) {
      svgEl("line", { x1: sx(lam0), x2: sx(lam0), y1: sy(c0 / Math.E), y2: bot, stroke: COL.green, "stroke-width": 2, "stroke-dasharray": "6 4" }, svg);
      svgEl("line", { x1: left, x2: sx(lam0), y1: sy(c0 / Math.E), y2: sy(c0 / Math.E), stroke: COL.green, "stroke-width": 1.4, "stroke-dasharray": "2 4" }, svg);
      svgEl("circle", { cx: sx(lam0), cy: sy(c0 / Math.E), r: 6, fill: COL.green }, svg);
      const t = svgEl("text", { x: sx(lam0) + 8, y: sy(c0 / Math.E) - 10, "font-weight": 600 }, svg);
      t.style.fill = COL.green; t.textContent = `λ = √(D/k) = ${f1(lam0)} µm  (c(0)/e)`;
    }
    // reference, steady, current
    const path = (fn: (i: number) => number, step = 5) => {
      let d = "";
      for (let i = 0; i <= NG; i += step) d += (i ? "L" : "M") + sx(xs[i]).toFixed(1) + "," + sy(fn(i)).toFixed(1);
      return d;
    };
    svgEl("path", { d: path(i => refFn(xs[i])), fill: "none", stroke: COL.grey, "stroke-width": 2.5, "stroke-opacity": 0.7 }, svg);
    svgEl("path", { d: path(i => f(xs[i])), fill: "none", stroke: COL.ink, "stroke-width": 2.5, "stroke-dasharray": "9 6" }, svg);
    svgEl("path", { d: path(i => G.c[i], 2), fill: "none", stroke: COL.blue, "stroke-width": 4 }, svg);
    // legend
    const lg = (txt: string, x: number, c: string, dash: string, w: number) => {
      svgEl("line", { x1: x, x2: x + 34, y1: 16, y2: 16, stroke: c, "stroke-width": w, "stroke-dasharray": dash }, svg);
      svgEl("text", { x: x + 42, y: 22 }, svg).textContent = txt;
    };
    lg("simulation c(x, t)", 360, COL.blue, "", 4);
    lg("exact steady state", 578, COL.ink, "9 6", 2.5);
    lg("reference (D 100, k 0.01, j 1)", 795, COL.grey, "", 2.5);
    svgEl("text", { x: (left + right) / 2, y: bot + 44, "text-anchor": "middle" }, svg).textContent = "position x (µm)   anterior (source) → posterior (wall)";

    // French flag: now (from the simulation), and at steady state (ticks)
    const fy = bot + 58, fh = 26;
    const xNow = thr.map(T => crossGrid(G.c, T));
    const xSt = thr.map(T => crossExact(f, T));
    const clampX = (x: number) => x < 0 ? 0 : Number.isFinite(x) ? x : L;
    const b1 = clampX(xNow[0]), b2 = Math.max(b1, clampX(xNow[1]));
    svgEl("text", { x: left - 6, y: fy + 20, "text-anchor": "end", class: "sm" }, svg).textContent = "fate";
    svgEl("rect", { x: sx(0), y: fy, width: sx(b1) - sx(0), height: fh, fill: COL.blue }, svg);
    svgEl("rect", { x: sx(b1), y: fy, width: sx(b2) - sx(b1), height: fh, fill: "#f2f2f2" }, svg);
    svgEl("rect", { x: sx(b2), y: fy, width: sx(L) - sx(b2), height: fh, fill: COL.red }, svg);
    svgEl("rect", { x: sx(0), y: fy, width: sx(L) - sx(0), height: fh, fill: "none", stroke: "#666" }, svg);
    xSt.forEach(x => {
      if (x < 0 || !Number.isFinite(x)) return;
      svgEl("path", { d: `M${sx(x)},${fy - 1}l-7,-11h14z`, fill: COL.ink }, svg);
    });
    refX.forEach(x => {
      if (x < 0 || !Number.isFinite(x)) return;
      svgEl("path", { d: `M${sx(x)},${fy + fh + 1}l-7,11h14z`, fill: COL.grey }, svg);
    });

    const fl = svgEl("text", { x: right, y: fy + fh + 32, "text-anchor": "end", class: "sm" }, svg);
    fl.textContent = "bar: now · ▼ steady state · ▲ reference";
    // readout
    ro.replaceChildren();
    const item = (k: string, v: string, cls = "") => { const s = htmlEl("span", cls, "", ro); htmlEl("b", "", k + " ", s); s.append(v); };
    item("t =", tTxt());
    item("λ = √(D/k) =", st.warm === "none" || st.warm === "all" ? `${f1(lam0)} µm` :
      `${f1(Math.sqrt(st.D / kAt(0)))} µm (anterior), ${f1(Math.sqrt(st.D / kAt(L)))} µm (posterior)`);
    item("lifetime 1/k =", st.warm === "none" ? `${(1 / st.k).toFixed(0)} min` : `${(1 / st.k).toFixed(0)} / ${(1 / (2 * st.k)).toFixed(0)} min (cold / warm)`);
    item("steady c(0) =", c0.toFixed(4));
    thr.forEach((T, i) => {
      const sh = xSt[i] >= 0 && Number.isFinite(xSt[i]) && refX[i] >= 0 ? ` (${xSt[i] - refX[i] >= 0 ? "+" : "−"}${f1(Math.abs(xSt[i] - refX[i]))} vs reference ${f1(refX[i])})` : "";
      item(`c = ${T}:`, `now ${xTxt(xNow[i])}, steady ${xTxt(xSt[i])}${sh}`);
    });
    if (st.done && Number.isFinite(st.tDone)) item("", `every point within 1 % of the steady state at t = ${st.tDone.toFixed(0)} min`, "tu-ok");
  }

  for (let i = 0; i <= NG; i++) kx[i] = kAt(xs[i]);
  draw();
  return {
    setActive(on: boolean) { st.active = on; if (!on) st.playing = false; draw(); },
  };
}

// ======================================================================= (b) Turing
const A0 = 0.1, B0 = 0.9, LBOX = 80;
const DXS = [2, 1, 0.5], DTS = [0.4, 0.2, 0.1, 0.05];

interface Snap { u: Float64Array; N: number; label: string }

function buildTuring(box: HTMLElement, ds: DOMStringMap): TabCtl {
  const start = {
    d: num(ds.d, 20, 1, 40),
    dx: DXS.includes(Number(ds.dx)) ? Number(ds.dx) : 1,
    dt: DTS.includes(Number(ds.dt)) ? Number(ds.dt) : 0.1,
    seed: Math.round(num(ds.seed, 1, 1, 1e6)),
    tend: num(ds.tend, 200, 10, 2000),
  };
  const st = { ...start, playing: false, active: false, stretch: false, endLogged: false };
  let S: Turing = newSim();
  let prev: Snap | null = null;
  const log: string[] = [];

  function newSim(): Turing { return makeTuring(Math.round(LBOX / st.dx), LBOX, st.dt, st.d, A0, B0, st.seed); }

  // ---- controls
  const r1 = htmlEl("div", "tu-row", "", box);
  const dS = slider(r1, "d = D_v / D_u", 1, 40, 0.5, st.d, v => { st.d = v; restart(true); }, v => v.toFixed(1));
  const pg = htmlEl("div", "tu-group", "", r1);
  const PRE: [string, number, number, number][] = [
    ["spots d = 20", 20, 1, 0.1], ["near threshold d = 10", 10, 1, 0.1], ["no Turing d = 5", 5, 1, 0.1], ["d = 5, Δt = 0.4", 5, 1, 0.4]];
  const pBtns = presets(pg, PRE.map(([label, d, dx, dt]) => ({ label, apply: () => setAll(d, dx, dt) })));
  const r2 = htmlEl("div", "tu-row", "", box);
  const xg = htmlEl("div", "tu-group", "", r2);
  htmlEl("span", "", "grid Δx:", xg);
  const xBtns = presets(xg, DXS.map(v => ({ label: String(v), apply: () => { st.dx = v; restart(true); } })));
  const tg = htmlEl("div", "tu-group", "", r2);
  htmlEl("span", "", "time step Δt:", tg);
  const tBtns = presets(tg, DTS.map(v => ({ label: String(v), apply: () => { st.dt = v; restart(true); } })));
  const cg = htmlEl("div", "tu-group", "", r2);
  const playB = button(cg, "▶ run", () => { st.playing = !st.playing; kick(); drawInfo(); });
  button(cg, "reset", () => { restart(false); });
  button(cg, "new noise", () => { st.seed++; restart(true); });
  checkbox(r2, "stretch colours", st.stretch, v => { st.stretch = v; drawField(); });

  // ---- layout: field | dispersion + readout | previous run + log
  const plots = htmlEl("div", "tu-plots", "", box);
  const c1 = htmlEl("div", "tu-col", "", plots); c1.style.width = "37%";
  const c2 = htmlEl("div", "tu-col", "", plots); c2.style.width = "40%";
  const c3 = htmlEl("div", "tu-col", "", plots); c3.style.width = "calc(23% - 20px)";
  const cap = htmlEl("div", "tu-cap", "", c1);
  const cv = htmlEl("canvas", "", "", c1);
  const CW = 470, CH = 410;
  cv.width = CW; cv.height = CH;
  cv.style.imageRendering = "auto";
  const off = document.createElement("canvas");

  const DW = 470, DH = 340;
  const svgD = svgEl("svg", { viewBox: `0 0 ${DW} ${DH}`, role: "img", "aria-label": "Dispersion relation: growth rate against wavenumber" }, c2);
  svgD.style.width = "100%";
  const ro = htmlEl("div", "qb-readout tu-read", "", c2);
  htmlEl("div", "tu-cap", "previous run", c3);
  const pcv = htmlEl("canvas", "", "", c3);
  pcv.width = 220; pcv.height = 220;
  const pcap = htmlEl("div", "qb-note", "–", c3);
  htmlEl("div", "tu-cap", "runs so far", c3);
  const logEl = htmlEl("div", "tu-log", "", c3);

  function setAll(d: number, dx: number, dt: number): void { st.d = d; st.dx = dx; st.dt = dt; dS.set(d); restart(true); }
  function restart(run: boolean): void {
    if (S.steps > 0) {
      prev = { u: S.u.slice(), N: S.N, label: `d = ${S.d}, Δx = ${S.dx}, Δt = ${S.dt}, t = ${S.t.toFixed(0)}: ${spotTxt(S)}` };
      drawPrev();
    }
    S = newSim(); st.endLogged = false;
    if (run && st.active) { st.playing = true; kick(); }
    drawAll();
  }

  function stats(s: Turing) {
    let mn = Infinity, mx = -Infinity, me = 0;
    for (const x of s.u) { if (x < mn) mn = x; if (x > mx) mx = x; me += x; }
    return { mn, mx, me: me / s.u.length };
  }
  function spotTxt(s: Turing): string {
    const q = stats(s);
    const n = countSpots(s.u, s.N);
    return q.mx - q.mn < 0.05 ? "uniform" : `${n} spot${n === 1 ? "" : "s"}`;
  }

  // ---- drawing the field
  function paint(ctx: CanvasRenderingContext2D, u: Float64Array, N: number, x0: number, y0: number, size: number, lo: number, hi: number): void {
    off.width = N; off.height = N;
    const octx = off.getContext("2d");
    if (!octx) return;
    const img = octx.createImageData(N, N);
    const sc = hi > lo ? 255 / (hi - lo) : 0;
    for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
      // as R's image(): first index along x, second along y (upwards)
      const v = u[i + N * j];
      const c = VLUT[Math.max(0, Math.min(255, Math.round((v - lo) * sc)))];
      const p = 4 * (i + N * (N - 1 - j));
      img.data[p] = c[0]; img.data[p + 1] = c[1]; img.data[p + 2] = c[2]; img.data[p + 3] = 255;
    }
    octx.putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(off, x0, y0, size, size);
  }
  function drawField(): void {
    const ctx = cv.getContext("2d");
    if (!ctx) return;
    const q = stats(S);
    const lo = st.stretch ? q.mn : 0, hi = st.stretch ? Math.max(q.mx, q.mn + 1e-12) : 3.3;
    ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, CW, CH);
    paint(ctx, S.u, S.N, 0, 5, 400, lo, hi);
    ctx.strokeStyle = "#444"; ctx.lineWidth = 1; ctx.strokeRect(0.5, 5.5, 399, 399);
    // colour bar
    for (let k = 0; k < 256; k++) {
      const c = VLUT[k];
      ctx.fillStyle = `rgb(${c[0]},${c[1]},${c[2]})`;
      ctx.fillRect(410, 5 + 400 - (k + 1) * 400 / 256, 16, 400 / 256 + 0.5);
    }
    ctx.strokeRect(410.5, 5.5, 16, 399);
    ctx.fillStyle = "#222"; ctx.font = "18px system-ui, sans-serif"; ctx.textAlign = "left";
    const lab = (v: number) => Math.abs(hi - lo) < 0.01 ? v.toPrecision(5) : v.toFixed(2);
    ctx.fillText(lab(hi), 430, 20);
    ctx.fillText(lab(lo), 430, 403);
    ctx.fillText("u", 432, 210);
    cap.textContent = `u at t = ${S.t.toFixed(1)}   (${S.N} × ${S.N}, Δx = ${S.dx}, Δt = ${S.dt})`;
  }
  function drawPrev(): void {
    const ctx = pcv.getContext("2d");
    if (!ctx || !prev) return;
    ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, 220, 220);
    paint(ctx, prev.u, prev.N, 0, 0, 220, 0, 3.3);
    pcap.textContent = prev.label;
  }

  // ---- dispersion relation and readout
  const qmax = 1.6;
  const dl = 58, dr = DW - 22, dtp = 80, dbt = DH - 50;
  const qx = scale(0, qmax, dl, dr);
  function dcrit(): number {
    const [J11, J12, J21, J22] = jac(A0, B0);
    const det = J11 * J22 - J12 * J21;
    // (d J11 + J22)^2 = 4 d det, larger root
    const qa = J11 * J11, qb = 2 * J11 * J22 - 4 * det, qc = J22 * J22;
    return (-qb + Math.sqrt(qb * qb - 4 * qa * qc)) / (2 * qa);
  }
  const DC = dcrit();
  function drawInfo(): void {
    playB.textContent = st.playing ? "❚❚ pause" : "▶ run";
    xBtns.forEach((b, i) => b.classList.toggle("active", DXS[i] === st.dx));
    tBtns.forEach((b, i) => b.classList.toggle("active", DTS[i] === st.dt));
    pBtns.forEach((b, i) => b.classList.toggle("active", PRE[i][1] === st.d && PRE[i][2] === st.dx && PRE[i][3] === st.dt));
    const q = stats(S);
    ro.replaceChildren();
    const item = (k: string, v: string, cls = "") => { const s = htmlEl("span", cls, "", ro); htmlEl("b", "", k + " ", s); s.append(v); };
    const sp = q.mx - q.mn < 0.05 ? null : countSpots(S.u, S.N);
    htmlEl("span", "tu-big", sp === null ? "spots: – (uniform)" : `spots: ${sp}`, ro);
    item("t =", S.t.toFixed(1));
    item("u from", `${q.mn.toFixed(3)} to ${q.mx.toFixed(3)}, mean ${q.me.toFixed(3)}`);
    // Euler factor of the reaction (q = 0 mode)
    const [p, qq, r, s] = jac(A0, B0);
    const tr = p + s, det = p * s - qq * r;
    const re = tr / 2, im = Math.sqrt(Math.max(0, det - tr * tr / 4));
    const ef = Math.hypot(1 + st.dt * re, st.dt * im);
    item("reaction eigenvalues", `${re.toFixed(1)} ± ${im.toFixed(3)}i;`);
    item("|1 + Δt λ| =", ef.toFixed(4) + (ef > 1 + 1e-9 ? "  > 1: Δt too large for the reaction" : "  ≤ 1"), ef > 1 + 1e-9 ? "tu-warn" : "");
  }
  function drawDisp(): void {
    svgD.replaceChildren();
    const ys: number[] = [], yi: number[] = [], qs: number[] = [];
    for (let i = 0; i <= 320; i++) {
      const q = i * qmax / 320; qs.push(q);
      ys.push(growth(q, st.d, A0, B0));
      yi.push(growthIMEX(q, st.d, A0, B0, st.dx, st.dt));
    }
    const gmax = Math.max(...ys, ...yi);
    const yhi = Math.max(0.3, Math.ceil(gmax * 10 + 0.5) / 10), ylo = -0.6;
    const gy = scale(ylo, yhi, dbt, dtp);
    const cl = (v: number) => gy(Math.max(ylo, Math.min(yhi, v)));
    svgEl("text", { x: 6, y: 22, class: "ttl" }, svgD).textContent = "dispersion relation";
    // band
    let qa = -1, qb = -1, qf = 0, gf = -Infinity;
    for (let i = 0; i <= 1600; i++) {          // band and fastest mode on a finer q grid
      const q = i / 1000, g = growth(q, st.d, A0, B0);
      if (g > 0) { if (qa < 0) qa = q; qb = q; }
      if (g > gf) { gf = g; qf = q; }
    }
    if (qa >= 0) {
      svgEl("rect", { x: qx(qa), y: dtp, width: qx(qb) - qx(qa), height: dbt - dtp, fill: COL.orange, "fill-opacity": 0.2 }, svgD);
    }
    for (let v = Math.ceil(ylo * 5) / 5; v <= yhi + 1e-9; v += 0.2) {
      const vv = Math.round(v * 10) / 10;
      svgEl("line", { x1: dl, x2: dr, y1: gy(vv), y2: gy(vv), class: vv === 0 ? "axis" : "grid", "stroke-width": vv === 0 ? 1.5 : 1 }, svgD);
      svgEl("text", { x: dl - 5, y: gy(vv) + 6, "text-anchor": "end", class: "sm" }, svgD).textContent = vv.toFixed(1);
    }
    for (let q = 0; q <= qmax + 1e-9; q += 0.4) {
      svgEl("text", { x: qx(q), y: dbt + 20, "text-anchor": "middle", class: "sm" }, svgD).textContent = q.toFixed(1);
    }
    svgEl("text", { x: (dl + dr) / 2, y: DH - 8, "text-anchor": "middle" }, svgD).textContent = "wavenumber q = 2π / wavelength";
    const yl = svgEl("text", { x: 14, y: (dtp + dbt) / 2, "text-anchor": "middle", transform: `rotate(-90 14 ${(dtp + dbt) / 2})` }, svgD);
    yl.textContent = "growth rate";
    // grid Nyquist
    const qn = Math.PI / st.dx;
    if (qn < qmax) {
      svgEl("line", { x1: qx(qn), x2: qx(qn), y1: dtp, y2: dbt, stroke: COL.purple, "stroke-width": 2, "stroke-dasharray": "2 4" }, svgD);
      const t = svgEl("text", { x: qx(qn) - 4, y: dbt - 8, "text-anchor": "end" }, svgD);
      t.style.fill = COL.purple; t.textContent = "π/Δx";
    }
    const pl = (v: number[], c: string, dash: string, w: number) => svgEl("polyline", {
      points: qs.map((q, i) => `${qx(q).toFixed(1)},${cl(v[i]).toFixed(1)}`).join(" "), fill: "none", stroke: c, "stroke-width": w, "stroke-dasharray": dash,
    }, svgD);
    pl(ys, COL.blue, "", 3.5);
    pl(yi, COL.red, "8 5", 2.5);
    if (gf > 0) {
      svgEl("circle", { cx: qx(qf), cy: gy(gf), r: 6, fill: COL.blue }, svgD);
      const t = svgEl("text", { x: qx(qf) + 8, y: gy(gf) - 8 }, svgD);
      t.textContent = `fastest: wavelength ${(2 * Math.PI / qf).toFixed(1)}`;
    }
    // legend + band text
    const lg = (txt: string, y: number, c: string, dash: string) => {
      svgEl("line", { x1: y, x2: y + 30, y1: 64, y2: 64, stroke: c, "stroke-width": 3, "stroke-dasharray": dash }, svgD);
      const t = svgEl("text", { x: y + 36, y: 70 }, svgD);
      t.style.fill = c; t.textContent = txt;
    };
    lg("exact", dl, COL.blue, "");
    lg("IMEX step on this grid", dl + 110, COL.red, "8 5");
    const bt = svgEl("text", { x: dl, y: 46, class: "sm" }, svgD);
    bt.textContent = qa >= 0 ? `unstable band q = ${qa.toFixed(3)}–${qb.toFixed(3)}   (needs d > d_c = ${DC.toFixed(2)})` : `no unstable band: d < d_c = ${DC.toFixed(2)}`;
  }
  function drawLog(): void { logEl.textContent = log.length ? log.slice(-7).join("\n") : "(a run is logged\n when it reaches\n t = " + st.tend + ")"; }
  function drawAll(): void {
    drawField(); drawDisp(); drawInfo(); drawLog();
  }

  // ---- animation: as many steps as fit into ~14 ms per frame
  let raf = 0;
  function kick(): void { if (!raf) raf = requestAnimationFrame(tick); }
  function tick(): void {
    raf = 0;
    if (!st.playing || !st.active) return;
    if (isShown(box)) {
      const t0 = performance.now();
      const target = st.endLogged ? Infinity : Math.round(st.tend / S.dt);
      do { stepTuring(S); } while (S.steps < target && performance.now() - t0 < 14);
      if (!st.endLogged && S.steps >= target) {
        st.endLogged = true; st.playing = false;
        log.push(`d=${S.d} Δx=${S.dx} Δt=${S.dt}${st.seed !== 1 ? " s" + st.seed : ""}: ${spotTxt(S)}`);
        drawLog();
      }
      drawField(); drawInfo();
    }
    if (st.playing) raf = requestAnimationFrame(tick);
  }

  drawAll();
  return {
    setActive(on: boolean) { st.active = on; if (!on) { st.playing = false; } drawInfo(); },
  };
}

autoMount(NAME, mount);

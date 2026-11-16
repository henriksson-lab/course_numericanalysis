// FTCS diffusion simulator for the diffusion-explicit unit: the flow-free chemotaxis
// chamber. Source reservoir c = 1 at x = 0; at x = L either a sink (c = 0) or a wall
// (no flux, ghost point c_{N+1} = c_{N-1}). Channel empty at t = 0.
//  (a) the profile c(x, t): FTCS on the grid (dots), the exact series solution (line,
//      optional) and the steady state (dashed), with the time in s and minutes;
//  (b) the distance from the steady state, max_i |c_i - c_inf(x_i)|, against time, both
//      on log axes, with L^2/D and the 5 % level marked. The exact distance curve is
//      dashed; earlier runs stay as faint traces, so doubling L shows a shift by x4.
// The animation runs at a fixed number of seconds per decade of simulated time, so
// early creeping and the late approach are both visible. Above r = 0.5 the zigzag grows.
// Auto-mounts on <div data-demo="ftcs">. Optional attributes:
//   data-d="1000" (µm²/s)  data-l="500" (µm)  data-n="50"  data-r="0.4"
//   data-bc="sink|wall"  data-exact="1|0"  data-autoplay="0|1"
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, autoMount } from "../lib/ui";

const NAME = "ftcs";
// Okabe-Ito, as in make-diffusion-explicit-figs.R
const COL = {
  ftcs: "#0072B2", exact: "#D55E00", green: "#009E73", sky: "#56B4E9", orange: "#E69F00",
  purple: "#CC79A7", grey: "#999999", ink: "#222222", warn: "#b00020",
};
const GHOST = [COL.sky, COL.orange, COL.purple];
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .35em; }
[data-demo="${NAME}"] .fd-row { display: flex; flex-wrap: wrap; gap: .3em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .fd-group { display: flex; gap: .35em; align-items: center; }
[data-demo="${NAME}"] .fd-row .qb-slider { flex: 1 1 400px; min-width: 300px; max-width: 560px; grid-template-columns: 7.2em 1fr 6.6em; }
[data-demo="${NAME}"] .fd-plots { display: flex; flex-wrap: wrap; gap: 8px; width: 100%; }
[data-demo="${NAME}"] .fd-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .fd-read { width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .fd-bad output { color: ${COL.warn}; font-weight: 600; }
[data-demo="${NAME}"] .qb-btn.fd-play { min-width: 5.5em; font-weight: 600; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 18px; padding: .25em .6em; }
.reveal [data-demo="${NAME}"] .fd-row .qb-slider { flex: 1 1 520px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

type BC = "sink" | "wall";

// ---- exact solution, dimensionless: xi = x/L in [0, 1], tau = t D / L^2 ----
function erfc(x: number): number {   // Abramowitz-Stegun 7.1.26, |error| < 1.5e-7 (drawing only)
  const s = x < 0 ? -1 : 1, z = Math.abs(x), t = 1 / (1 + 0.3275911 * z);
  const y = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
  const e = y * Math.exp(-z * z);
  return s > 0 ? e : 2 - e;
}
const steady = (bc: BC, xi: number): number => (bc === "sink" ? 1 - xi : 1);
/** c(xi, tau). For tau < 1e-3 the far end is out of reach: the half-line solution erfc. */
function exact(bc: BC, xi: number, tau: number): number {
  if (tau <= 0) return xi === 0 ? 1 : 0;
  if (tau < 1e-3) return erfc(xi / (2 * Math.sqrt(tau)));
  const nmax = Math.ceil(Math.sqrt(40 / (Math.PI * Math.PI * tau))) + 3;
  let s = 0;
  if (bc === "sink") {
    // c = 1 - xi - sum 2/(n pi) sin(n pi xi) exp(-n^2 pi^2 tau)
    for (let n = 1; n <= nmax; n++)
      s += 2 / (n * Math.PI) * Math.sin(n * Math.PI * xi) * Math.exp(-n * n * Math.PI * Math.PI * tau);
    return 1 - xi - s;
  }
  // wall at xi = 1: c = 1 - sum 4/(k pi) sin(k pi xi / 2) exp(-k^2 pi^2 tau / 4), k odd
  for (let k = 1; k <= 2 * nmax + 1; k += 2)
    s += 4 / (k * Math.PI) * Math.sin(k * Math.PI * xi / 2) * Math.exp(-k * k * Math.PI * Math.PI * tau / 4);
  return 1 - s;
}
/** Exact max_x |c - c_inf| as a function of tau, on a log grid (computed once per BC). */
const distCache: Partial<Record<BC, { lt: number[]; d: number[] }>> = {};
function exactDist(bc: BC): { lt: number[]; d: number[] } {
  const hit = distCache[bc];
  if (hit) return hit;
  const lt: number[] = [], d: number[] = [];
  for (let l = -7; l <= 0.75 + 1e-9; l += 0.025) {
    const tau = 10 ** l;
    const xs: number[] = [];
    for (let i = 0; i <= 400; i++) xs.push(i / 400);
    for (let k = 1; k <= 40; k++) xs.push(k * 0.25 * Math.sqrt(tau));   // resolve the front
    let m = 0;
    for (const xi of xs) if (xi <= 1) m = Math.max(m, Math.abs(exact(bc, xi, tau) - steady(bc, xi)));
    lt.push(l); d.push(m);
  }
  return (distCache[bc] = { lt, d });
}

// ---- formatting ----
function fmtSec(s: number): string {
  if (s < 10) return s.toFixed(2) + " s";
  if (s < 1000) return s.toFixed(1) + " s";
  return s.toFixed(0) + " s";
}
function fmtMin(s: number): string {
  const m = s / 60;
  if (m < 10) return m.toFixed(2) + " min";
  if (m < 100) return m.toFixed(1) + " min";
  return m.toFixed(0) + " min";
}
const fmtTime = (s: number): string => `${fmtMin(s)}  (${fmtSec(s)})`;
function sci(v: number, d = 2): string {
  if (!Number.isFinite(v)) return String(v);
  if (v !== 0 && (Math.abs(v) >= 1e4 || Math.abs(v) < 1e-3)) return v.toExponential(d - 1).replace("e+", "e");
  return String(+v.toPrecision(d + 1));
}
const sigD = (v: number): number => Number(v.toPrecision(3));
const sup: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const pow10 = (k: number): string => k === 0 ? "1" : "10" + String(k).split("").map(c => sup[c]).join("");

interface Trace { pts: [number, number][]; col: string; label: string; l2d: number }

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (v: string | undefined, def: number, lo: number, hi: number) => {
    const x = Number(v);
    return v !== undefined && Number.isFinite(x) ? Math.min(hi, Math.max(lo, x)) : def;
  };
  const START = {
    D: num(ds.d, 1000, 50, 3000), L: num(ds.l, 500, 100, 2000), N: Math.round(num(ds.n, 50, 4, 200)),
    r: num(ds.r, 0.4, 0.05, 0.7), bc: (ds.bc === "wall" ? "wall" : "sink") as BC, exact: ds.exact !== "0",
  };
  const st = { ...START };
  let secPerDecade = 2;

  // ---- simulation state ----
  let c = new Float64Array(0), cn = new Float64Array(0);
  let nstep = 0, dx = 0, dt = 0, t5 = -1, blown = false, finished = false;
  let trace: [number, number][] = [];
  let ghosts: Trace[] = [];
  let ghostNo = 0;

  const dist = (): number => {
    let m = 0;
    for (let i = 0; i <= st.N; i++) {
      const v = Math.abs(c[i] - steady(st.bc, i / st.N));
      if (!(v <= m)) m = v;               // NaN propagates
    }
    return m;
  };
  function restart(keepGhost: boolean): void {
    if (keepGhost && trace.length > 5 && nstep > 0) {
      ghosts.push({ pts: trace, col: GHOST[ghostNo++ % GHOST.length], label: lastLabel, l2d: lastL2D });
      if (ghosts.length > 3) ghosts.shift();
    }
    dx = st.L / st.N; dt = st.r * dx * dx / st.D;
    c = new Float64Array(st.N + 1); cn = new Float64Array(st.N + 1);
    c[0] = 1;
    nstep = 0; t5 = -1; blown = false; finished = false;
    trace = [];
    lastL2D = st.L * st.L / st.D; lastLabel = label();
  }
  let lastL2D = 0, lastLabel = "";
  const label = () => `D ${sigD(st.D)} · L ${st.L} · N ${st.N} · r ${st.r}${st.bc === "wall" ? " · wall" : ""}`;
  const tauEnd = () => (st.bc === "sink" ? 1.0 : 3.5);

  /** One FTCS step: interior from old values, wall update with the ghost point. */
  function step(): void {
    const N = st.N, r = st.r;
    for (let i = 1; i < N; i++) cn[i] = c[i] + r * (c[i - 1] - 2 * c[i] + c[i + 1]);
    cn[0] = 1;
    cn[N] = st.bc === "sink" ? 0 : c[N] + 2 * r * (c[N - 1] - c[N]);
    const tmp = c; c = cn; cn = tmp;
    nstep++;
    const d = dist();
    if (t5 < 0 && d <= 0.05) t5 = nstep * dt;
    if (!Number.isFinite(d) || d > 1e300) blown = true;
  }
  const record = () => { const t = nstep * dt; if (t > 0) trace.push([t, dist()]); };

  // ---- controls ----
  const row1 = htmlEl("div", "fd-row", "", host);
  const g1 = htmlEl("div", "fd-group", "", row1);
  htmlEl("span", "", "molecule:", g1);
  const molBtns = presets(g1, [
    { label: "small molecule", apply: () => { st.D = 1000; dS.set(3); restart(true); draw(); } },
    { label: "protein", apply: () => { st.D = 100; dS.set(2); restart(true); draw(); } },
  ]);
  const g2 = htmlEl("div", "fd-group", "", row1);
  htmlEl("span", "", "at x = L:", g2);
  const bcBtns = presets(g2, [
    { label: "sink c = 0", apply: () => { st.bc = "sink"; restart(true); draw(); } },
    { label: "wall (no flux)", apply: () => { st.bc = "wall"; restart(true); draw(); } },
  ]);
  checkbox(row1, "exact solution", st.exact, v => { st.exact = v; draw(); });

  const row2 = htmlEl("div", "fd-row", "", host);
  const dS = slider(row2, "D (µm²/s)", Math.log10(50), Math.log10(3000), 0.01, Math.log10(st.D),
    v => { st.D = sigD(10 ** v); restart(true); draw(); }, v => String(sigD(10 ** v)));
  const lS = slider(row2, "L (µm)", 100, 2000, 50, st.L, v => { st.L = v; restart(true); draw(); }, v => v + " µm");
  const row3 = htmlEl("div", "fd-row", "", host);
  const nS = slider(row3, "N intervals", 4, 200, 1, st.N, v => { st.N = v; restart(true); draw(); }, v => String(v));
  const rS = slider(row3, "r = DΔt/Δx²", 0.05, 0.7, 0.01, st.r, v => { st.r = v; restart(true); draw(); },
    v => v.toFixed(2) + (v > 0.5 + 1e-9 ? " !" : ""));
  // tick mark at r = 0.5 (shown by the browser under the range input)
  const dl = htmlEl("datalist", "", "", row3);
  dl.id = `fd-rmarks-${Math.random().toString(36).slice(2, 8)}`;
  for (const v of [0.25, 0.5]) { const o = htmlEl("option", "", "", dl); o.value = String(v); }
  rS.input.setAttribute("list", dl.id);
  const rLab = rS.input.parentElement as HTMLElement;

  const row4 = htmlEl("div", "fd-row", "", host);
  const g4 = htmlEl("div", "fd-group", "", row4);
  const playBtn = button(g4, "▶ run", () => { if (playing) pause(); else play(); }, "fd-play");
  button(g4, "1 step", () => { pause(); if (!blown) { step(); record(); } draw(); });
  button(g4, "restart run", () => { restart(false); draw(); });
  const g5 = htmlEl("div", "fd-group", "", row4);
  htmlEl("span", "", "r =", g5);
  const rBtns = presets(g5, [0.25, 0.4, 0.5, 0.52, 0.6].map(v => ({
    label: String(v), apply: () => { st.r = v; rS.set(v); restart(true); draw(); },
  })));
  const g6 = htmlEl("div", "fd-group", "", row4);
  htmlEl("span", "", "speed:", g6);
  const spBtns = presets(g6, [
    { label: "slow", apply: () => { secPerDecade = 4; } },
    { label: "normal", apply: () => { secPerDecade = 2; } },
    { label: "fast", apply: () => { secPerDecade = 0.7; } },
  ]);
  spBtns[1].classList.add("active");
  button(row4, "clear traces", () => { ghosts = []; draw(); });
  button(row4, "reset all", () => {
    pause();
    Object.assign(st, START);
    dS.set(Math.log10(st.D)); lS.set(st.L); nS.set(st.N); rS.set(st.r);
    secPerDecade = 2; spBtns.forEach((b, i) => b.classList.toggle("active", i === 1));
    ghosts = []; ghostNo = 0;
    restart(false); draw();
  });

  // ---- plots ----
  const plots = htmlEl("div", "fd-plots", "", host);
  const H = 330, WA = 570, WB = 570, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Concentration profile in the channel" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 4px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Distance from steady state against time, log axes" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 4px)`;
  const ro = htmlEl("div", "qb-readout fd-read", "", host);

  const top = 64, bot = H - 48;
  // (a): x in µm, c in [-0.2, 1.2]
  const AX0 = 92, AX1 = WA - 46;
  const ay = scale(-0.15, 1.12, bot, top);
  // (b): log10 t (s) in [-0.5, 5], log10 distance in [-5, 1]
  const BX0 = 70, BX1 = WB - 14, LT0 = -0.5, LT1 = 5, LD0 = -5, LD1 = 1;
  const bx = scale(LT0, LT1, BX0, BX1), by = scale(LD0, LD1, bot, top);
  const cy = (v: number) => Math.max(LD0 - 0.3, Math.min(LD1 + 0.15, v));
  const uid = Math.random().toString(36).slice(2, 8);

  function drawProfile(): void {
    const svg = svgA, L = st.L, N = st.N;
    const ax = scale(0, L, AX0, AX1);
    svg.replaceChildren();
    const t = nstep * dt;
    const ttl = svgEl("text", { x: 8, y: 24, class: "ttl" }, svg);
    ttl.textContent = `(a) c(x, t)   t = ${fmtTime(t)}`;
    for (const v of [0, 0.5, 1]) {
      svgEl("line", { x1: AX0, x2: AX1, y1: ay(v), y2: ay(v), class: v === 0 ? "axis" : "grid" }, svg);
      svgEl("text", { x: AX0 - 30, y: ay(v) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = String(v);
    }
    const stepX = [20, 25, 50, 100, 200, 250, 500].find(s => L / s <= 5.01) ?? 500;
    for (let x = 0; x <= L + 1e-9; x += stepX) {
      svgEl("line", { x1: ax(x), x2: ax(x), y1: bot, y2: bot + 5, class: "axis" }, svg);
      svgEl("text", { x: ax(x), y: bot + 22, "text-anchor": "middle", class: "sm" }, svg).textContent = String(x);
    }
    svgEl("text", { x: (AX0 + AX1) / 2, y: H - 6, "text-anchor": "middle" }, svg).textContent = "position x (µm)";
    const ym = (top + bot) / 2;
    svgEl("text", { x: 16, y: ym, "text-anchor": "middle", transform: `rotate(-90 16 ${ym})` }, svg).textContent = "concentration c";
    // reservoirs / wall
    svgEl("rect", { x: AX0 - 24, y: ay(1), width: 20, height: ay(0) - ay(1), fill: "#0072B233", stroke: COL.ftcs, "stroke-width": 2 }, svg);
    const vlab = (x: number, txt: string) => {
      const yy = (ay(0) + ay(1)) / 2;
      svgEl("text", { x, y: yy + 6, "text-anchor": "middle", transform: `rotate(-90 ${x} ${yy})`, class: "sm" }, svg).textContent = txt;
    };
    vlab(AX0 - 14, "source");
    if (st.bc === "sink") {
      svgEl("rect", { x: AX1 + 4, y: ay(1), width: 20, height: ay(0) - ay(1), fill: "#f2f2f2", stroke: COL.grey, "stroke-width": 2 }, svg);
      vlab(AX1 + 14, "sink");
    } else {
      svgEl("rect", { x: AX1 + 4, y: ay(1.1), width: 12, height: ay(-0.1) - ay(1.1), fill: "#555" }, svg);
      svgEl("text", { x: AX1 + 10, y: ay(1.1) - 6, "text-anchor": "middle", class: "sm" }, svg).textContent = "wall";
    }
    // steady state
    const ss = st.bc === "sink" ? [[0, 1], [L, 0]] : [[0, 1], [L, 1]];
    svgEl("line", { x1: ax(ss[0][0]), y1: ay(ss[0][1]), x2: ax(ss[1][0]), y2: ay(ss[1][1]),
      stroke: COL.grey, "stroke-width": 2.5, "stroke-dasharray": "8 6" }, svg);
    // exact
    const tau = t * st.D / (L * L);
    if (st.exact) {
      const pts: string[] = [];
      for (let k = 0; k <= 250; k++) {
        const xi = k / 250;
        pts.push(`${ax(xi * L).toFixed(1)},${ay(exact(st.bc, xi, tau)).toFixed(1)}`);
      }
      svgEl("polyline", { points: pts.join(" "), fill: "none", stroke: COL.exact, "stroke-width": 2.5 }, svg);
    }
    // FTCS, clipped to the plot
    const clipId = "fd-clip-" + uid;
    const cp = svgEl("clipPath", { id: clipId }, svg);
    svgEl("rect", { x: AX0 - 8, y: top - 8, width: AX1 - AX0 + 16, height: bot - top + 16 }, cp);
    const g = svgEl("g", { "clip-path": `url(#${clipId})` }, svg);
    const clampC = (v: number) => (Number.isFinite(v) ? Math.max(-0.14, Math.min(1.11, v)) : (v > 0 ? 1.11 : -0.14));
    const pts = Array.from(c, (v, i) => `${ax(i * dx).toFixed(1)},${ay(clampC(v)).toFixed(1)}`).join(" ");
    svgEl("polyline", { points: pts, fill: "none", stroke: COL.ftcs, "stroke-width": 1.8, opacity: 0.7 }, g);
    const rad = N <= 60 ? 4.5 : N <= 120 ? 3 : 2;
    c.forEach((v, i) => {
      const out = !(v >= -1e-9 && v <= 1 + 1e-9);
      svgEl("circle", { cx: ax(i * dx), cy: ay(clampC(v)), r: rad, fill: out ? COL.warn : COL.ftcs }, g);
    });
    // legend
    legend(svg, AX0 - 60, [["● FTCS", COL.ftcs], ...(st.exact ? [["— exact series", COL.exact]] : []),
      [st.bc === "sink" ? "- - steady 1 − x/L" : "- - steady c = 1", "#666"]] as [string, string][]);
    // instability warning
    let lo = Infinity, hi = -Infinity;
    for (const v of c) { if (v < lo || Number.isNaN(v)) lo = v; if (v > hi || Number.isNaN(v)) hi = v; }
    if (blown || lo < -1e-9 || hi > 1 + 1e-9) {
      const w = svgEl("text", { x: AX0 + 8, y: ay(-0.12), "font-weight": 600 }, svg);
      w.style.fill = COL.warn;
      w.textContent = blown ? "overflow: values reached ±∞ — unstable"
        : `outside [0, 1]: c from ${sci(lo)} to ${sci(hi)}`;
    }
  }

  /** One-line legend under the panel title. */
  function legend(svg: SVGSVGElement, x: number, items: [string, string][]): void {
    const e = svgEl("text", { x, y: 50, "font-weight": 600 }, svg);
    items.forEach(([txt, col], i) => {
      const t = svgEl("tspan", { dx: i ? 18 : 0 }, e);
      t.style.fill = col; t.textContent = txt;
    });
  }

  function drawDist(): void {
    const svg = svgB;
    svg.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svg).textContent = "(b) distance from steady state";
    for (let k = LD0; k <= LD1; k++) {
      svgEl("line", { x1: BX0, x2: BX1, y1: by(k), y2: by(k), class: k === LD0 ? "axis" : "grid" }, svg);
      svgEl("text", { x: BX0 - 6, y: by(k) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = pow10(k);
    }
    const ticks: [number, string][] = [[1, "1 s"], [10, "10 s"], [60, "1 min"], [600, "10 min"], [3600, "1 h"], [36000, "10 h"]];
    for (const [s, lab] of ticks) {
      const X = bx(Math.log10(s));
      svgEl("line", { x1: X, x2: X, y1: top, y2: bot, class: "grid" }, svg);
      svgEl("text", { x: X, y: bot + 22, "text-anchor": "middle", class: "sm" }, svg).textContent = lab;
    }
    svgEl("text", { x: (BX0 + BX1) / 2, y: H - 6, "text-anchor": "middle" }, svg).textContent = "time t (log scale)";
    const ym = (top + bot) / 2;
    svgEl("text", { x: 16, y: ym, "text-anchor": "middle", transform: `rotate(-90 16 ${ym})` }, svg).textContent = "max |c − c∞|";
    // 5 % level
    svgEl("line", { x1: BX0, x2: BX1, y1: by(Math.log10(0.05)), y2: by(Math.log10(0.05)), stroke: "#666", "stroke-width": 1.5, "stroke-dasharray": "3 4" }, svg);
    svgEl("text", { x: BX0 + 6, y: by(Math.log10(0.05)) - 6, class: "sm" }, svg).textContent = "5 %";

    const clipId = "fd-clipb-" + uid;
    const cp = svgEl("clipPath", { id: clipId }, svg);
    svgEl("rect", { x: BX0, y: top - 6, width: BX1 - BX0, height: bot - top + 12 }, cp);
    const g = svgEl("g", { "clip-path": `url(#${clipId})` }, svg);
    const path = (p: [number, number][], stroke: string, w: number, dash = "", op = 1) => {
      const s = p.filter(([t]) => t > 0).map(([t, d]) => `${bx(Math.log10(t)).toFixed(1)},${by(cy(Math.log10(Math.max(d, 1e-30)))).toFixed(1)}`);
      if (s.length > 1) svgEl("polyline", { points: s.join(" "), fill: "none", stroke, "stroke-width": w, "stroke-dasharray": dash, opacity: op }, g);
    };
    // ghosts
    let gy = bot - 10;
    for (const gh of ghosts) {
      path(gh.pts, gh.col, 3, "", 0.8);
      const X = bx(Math.log10(gh.l2d));
      svgEl("line", { x1: X, x2: X, y1: by(LD1), y2: by(LD1) + 18, stroke: gh.col, "stroke-width": 3 }, g);
      const e = svgEl("text", { x: BX0 + 8, y: gy, "font-weight": 600 }, svg);
      e.style.fill = gh.col; e.textContent = gh.label; gy -= 22;
    }
    // exact distance for the current D, L
    const l2d = st.L * st.L / st.D;
    if (st.exact) {
      const ed = exactDist(st.bc);
      path(ed.lt.map((l, i) => [10 ** l * l2d, ed.d[i]] as [number, number]), COL.exact, 2.2, "7 5");
    }
    // L^2/D marker
    const X = bx(Math.log10(l2d));
    svgEl("line", { x1: X, x2: X, y1: top, y2: bot, stroke: COL.green, "stroke-width": 2.5 }, g);
    if (X > BX0 && X < BX1) {
      const lx = X > BX1 - 90 ? X - 6 : X + 6;
      const e = svgEl("text", { x: lx, y: top + 18, "text-anchor": X > BX1 - 90 ? "end" : "start", "font-weight": 600 }, svg);
      e.style.fill = COL.green; e.textContent = `L²/D = ${fmtSec(l2d)}`;
    }
    // current trace
    path(trace, COL.ftcs, 3.5);
    const last = trace[trace.length - 1];
    if (last) {
      const ld = Math.log10(Math.max(last[1], 1e-30));
      svgEl("circle", { cx: bx(Math.max(LT0, Math.min(LT1, Math.log10(last[0])))), cy: by(cy(ld)), r: 6, fill: COL.ftcs, stroke: "#fff", "stroke-width": 1.5 }, svg);
      if (ld > LD1) {
        const e = svgEl("text", { x: BX0 + 6, y: top + 44, "font-weight": 600 }, svg);
        e.style.fill = COL.warn; e.textContent = `↑ off scale: ${sci(last[1])}`;
      }
    }
    if (t5 > 0) {
      const X5 = bx(Math.log10(t5));
      svgEl("circle", { cx: X5, cy: by(Math.log10(0.05)), r: 6, fill: "none", stroke: COL.ink, "stroke-width": 2 }, g);
    }
    legend(svg, BX0 - 40, [["— FTCS, this run", COL.ftcs], ...(st.exact ? [["- - exact", COL.exact]] : []),
      ["│ L²/D", COL.green]] as [string, string][]);
  }

  function draw(): void {
    molBtns.forEach((b, i) => b.classList.toggle("active", st.D === (i === 0 ? 1000 : 100)));
    bcBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.bc === "sink")));
    rBtns.forEach((b, i) => b.classList.toggle("active", Math.abs([0.25, 0.4, 0.5, 0.52, 0.6][i] - st.r) < 1e-9));
    rLab.classList.toggle("fd-bad", st.r > 0.5 + 1e-9);
    playBtn.textContent = playing ? "❚❚ pause" : "▶ run";
    drawProfile();
    drawDist();
    const t = nstep * dt, l2d = st.L * st.L / st.D;
    const f = 1 - 4 * st.r;
    const stab = st.r <= 0.5 + 1e-9 ? `stable (zigzag factor 1 − 4r = ${f.toFixed(2)})`
      : `UNSTABLE: zigzag factor 1 − 4r = ${f.toFixed(2)}, grows ${(Math.abs(f)).toFixed(2)}× per step`;
    ro.textContent =
      `Δx = ${sci(dx, 3)} µm   Δt = ${sci(dt, 3)} s   steps = ${nstep}   tD/L² = ${(t / l2d).toFixed(4)}   L²/D = ${fmtSec(l2d)}   ${stab}\n` +
      (t5 > 0 ? `within 5 % of steady state after ${fmtTime(t5)} = ${(t5 / l2d).toFixed(4)} L²/D  (${Math.round(t5 / dt)} steps)`
        : blown ? "the run blew up: values overflowed to ±∞ (stopped)"
        : "not yet within 5 % of the steady state");
  }

  // ---- animation: a fixed number of real seconds per decade of simulated time ----
  let playing = false, raf = 0, lastTs = 0;
  function play(): void {
    if (finished || blown) restart(false);
    playing = true; lastTs = 0;
    if (!raf) raf = window.requestAnimationFrame(tick);
    draw();
  }
  function pause(): void {
    playing = false;
    if (raf) { window.cancelAnimationFrame(raf); raf = 0; }
    draw();
  }
  function tick(ts: number): void {
    raf = 0;
    if (!playing) return;
    const el = lastTs ? Math.min(0.1, (ts - lastTs) / 1000) : 1 / 60;
    lastTs = ts;
    const t = nstep * dt;
    const target = Math.max(t + dt, t * 10 ** (el / secPerDecade));
    const maxSteps = Math.max(1, Math.floor(4e6 / st.N));
    let k = 0;
    while (nstep * dt < target - 1e-12 * dt && k < maxSteps && !blown) { step(); k++; }
    record();
    const tau = nstep * dt * st.D / (st.L * st.L);
    if (blown || (tau >= tauEnd() && st.r <= 0.5 + 1e-9) || Math.log10(nstep * dt) > LT1) { finished = true; playing = false; }
    draw();
    if (playing) raf = window.requestAnimationFrame(tick);
  }

  restart(false);
  draw();
  if (ds.autoplay === "1") play();
}

autoMount(NAME, mount);

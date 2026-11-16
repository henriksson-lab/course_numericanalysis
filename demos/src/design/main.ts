// Sampling-design and uncertainty explorer for the identifiability unit.
// The oral one-compartment PK model (D = 100 mg, true ka = 1, ke = 0.2 /h, V = 20 L)
// sampled at draggable times. Each sample keeps its own noise value, taken from the
// course data (noise SD 0.1 mg/L), so the presets "all 12" and "from 3 h" reproduce the
// lecture's fits exactly; with the noise switched off the data are the true curve and
// the SEs use sigma = 0.1 (the lecture's se_design).
//  (a) the data, the true and the fitted curve; drag a sample along the time axis;
//  (b) Delta S / s^2 over (ka, ke) with V re-fitted, and the Wald ellipse at 3.84 from
//      the (ka, ke) block of s^2 (J^T J)^{-1};
//  (c) the profile of ka (ke and V re-fitted for each fixed ka) and the Wald parabola;
//  below: singular values of the log-parameter Jacobian, SEs, correlations.
// Grey in (b) and (c): the full 12-sample design, for comparison.
// Auto-mounts on <div data-demo="design">. Optional attributes:
//   data-preset="full|late3|spread|early|late"  data-noise="1|0"  data-zoom="0|1"
import { svgEl, htmlEl, presets, checkbox, button, scale, svgPoint, draggable, autoMount } from "../lib/ui";
import {
  T_DATA, EPS_DATA, SIGMA, TRUTH, pk, analyse, profileKa, sseKaKe, Analysis,
} from "./math";

const NAME = "design";
// Okabe-Ito, as in make-identifiability-figs.R
const COL = { blue: "#0072B2", red: "#D55E00", orange: "#E69F00", green: "#009E73", grey: "#999999", ink: "#222" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; font-size: 18px; }
[data-demo="${NAME}"] .dz-row { display: flex; flex-wrap: wrap; gap: .4em 1.2em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .dz-main { width: 100%; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 19px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .dz-bottom { display: flex; flex-wrap: wrap; gap: .6em 1.2em; width: 100%; align-items: flex-start; }
[data-demo="${NAME}"] .dz-sv { width: 380px; max-width: 100%; height: auto; flex: 0 1 380px; }
[data-demo="${NAME}"] table.dz-t { border-collapse: collapse; font-family: ui-monospace, monospace; font-size: 18px; }
[data-demo="${NAME}"] table.dz-t th, [data-demo="${NAME}"] table.dz-t td { padding: .08em .5em; text-align: right; }
[data-demo="${NAME}"] table.dz-t th { font-family: system-ui, sans-serif; font-weight: 600; border-bottom: 1px solid #bbb; }
[data-demo="${NAME}"] table.dz-t td:first-child { text-align: left; font-family: system-ui, sans-serif; }
[data-demo="${NAME}"] .dz-info { display: flex; flex-direction: column; gap: .25em; font-size: 18px; flex: 1 1 260px; }
[data-demo="${NAME}"] .dz-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .dz-key { display: inline-block; width: 22px; height: 0; border-top: 3px solid; vertical-align: middle; margin-right: .3em; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 18px; }
.reveal [data-demo="${NAME}"] .qb-btn { font-size: 18px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- designs: each sample = (index of its noise value, time) ----
interface Sample { k: number; t: number }
const PRESETS: Record<string, { label: string; k: number[]; t: number[] }> = {
  full: { label: "all 12 samples", k: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], t: T_DATA },
  late3: { label: "from 3 h (6)", k: [6, 7, 8, 9, 10, 11], t: [3, 4, 6, 8, 10, 12] },
  spread: { label: "spread (6)", k: [1, 3, 5, 7, 9, 11], t: [0.5, 1, 2, 4, 8, 12] },
  early: { label: "first 3 h (6)", k: [0, 1, 3, 4, 5, 6], t: [0.25, 0.5, 1, 1.5, 2, 3] },
  late: { label: "all late (6)", k: [6, 7, 8, 9, 10, 11], t: [7, 8, 9, 10, 11, 12] },
};
const PRESET_KEYS = Object.keys(PRESETS);
const T_MIN = 0.1, T_MAX = 15;
const CRIT = 3.84;

const yOf = (s: Sample, noisy: boolean): number =>
  pk(s.t, TRUTH[0], TRUTH[1], TRUTH[2]) + (noisy ? SIGMA * EPS_DATA[s.k] : 0);

// profile grid of ka, log-spaced 0.05 ... 50
const KA0 = 0.05, KA1 = 50, NKA = 97;
const KA_GRID = Array.from({ length: NKA }, (_, i) => 10 ** (Math.log10(KA0) + i * Math.log10(KA1 / KA0) / (NKA - 1)));

interface Result {
  a: Analysis;
  Smin: number;
  grid: number[];                    // KA_GRID plus the estimates of ka and ke (the two wells)
  prof: number[];                    // Delta = (PL - Smin) / s2 on grid
  lo: number | null; hi: number | null; // profile interval (null = no bound in range)
  ts: number[]; ys: number[];
}

export function compute(samples: Sample[], noisy: boolean): Result {
  const ts = samples.map(s => s.t), ys = samples.map(s => yOf(s, noisy));
  const a = analyse(ts, ys, noisy);
  const grid = KA_GRID.concat(a.th.slice(0, 2).filter(v => v > KA0 && v < KA1)).sort((x, y) => x - y);
  const pl = grid.map(ka => profileKa(ts, ys, ka));
  const Smin = Math.min(a.S, ...pl);
  const delta = (v: number) => (v - Smin) / a.s2;
  const prof = pl.map(delta);
  // profile interval: walk out from the grid point nearest the estimate, refine by bisection
  let lo: number | null = null, hi: number | null = null;
  const lk = Math.log(a.th[0]);
  if (a.th[0] > KA0 && a.th[0] < KA1) {
    let i0 = 0;
    grid.forEach((k, i) => { if (Math.abs(Math.log(k) - lk) < Math.abs(Math.log(grid[i0]) - lk)) i0 = i; });
    const f = (l: number) => delta(profileKa(ts, ys, Math.exp(l))) - CRIT;
    const bisect = (a0: number, b0: number) => {      // f(a0) < 0 < f(b0)
      let x0 = a0, x1 = b0;
      for (let it = 0; it < 30; it++) { const m = (x0 + x1) / 2; if (f(m) < 0) x0 = m; else x1 = m; }
      return Math.exp((x0 + x1) / 2);
    };
    for (let i = i0 + 1; i < grid.length; i++) if (prof[i] > CRIT) { hi = bisect(Math.max(lk, Math.log(grid[i - 1])), Math.log(grid[i])); break; }
    for (let i = i0 - 1; i >= 0; i--) if (prof[i] > CRIT) { lo = bisect(Math.min(lk, Math.log(grid[i + 1])), Math.log(grid[i])); break; }
  }
  return { a, Smin, grid, prof, lo, hi, ts, ys };
}

// ---- formatting ----
function fmt(v: number, d = 3): string {
  if (!Number.isFinite(v)) return "∞";
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e5 || a < 1e-3)) { const [m, e] = v.toExponential(1).split("e"); return `${m}e${Number(e)}`; }
  return v.toPrecision(d);
}
const fmtSig = (v: number): string => Number.isFinite(v) && Math.abs(v) < 1e4 ? v.toPrecision(3) : "∞";

/** Contour lines of a grid z[i][j] (x index i, y index j) at level c, as an SVG path (marching squares). */
function contourPath(z: number[][], c: number, X: (i: number) => number, Y: (j: number) => number): string {
  const nx = z.length, ny = z[0].length;
  let d = "";
  const lerp = (a: number, b: number) => (c - a) / (b - a);
  for (let i = 0; i < nx - 1; i++) for (let j = 0; j < ny - 1; j++) {
    const v = [z[i][j], z[i + 1][j], z[i + 1][j + 1], z[i][j + 1]];
    const pts: [number, number][] = [];
    const edges: [number, number, number, number][] = [[0, 1, 0, 0], [1, 2, 1, 0], [3, 2, 0, 1], [0, 3, 0, 0]];
    // corners: 0 (i,j) 1 (i+1,j) 2 (i+1,j+1) 3 (i,j+1)
    const cx = [i, i + 1, i + 1, i], cy = [j, j, j + 1, j + 1];
    for (const [p, q] of edges) {
      if ((v[p] < c) !== (v[q] < c)) {
        const f = lerp(v[p], v[q]);
        pts.push([X(cx[p] + f * (cx[q] - cx[p])), Y(cy[p] + f * (cy[q] - cy[p]))]);
      }
    }
    for (let k = 0; k + 1 < pts.length; k += 2)
      d += `M${pts[k][0].toFixed(1)},${pts[k][1].toFixed(1)}L${pts[k + 1][0].toFixed(1)},${pts[k + 1][1].toFixed(1)}`;
  }
  return d;
}

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const uid = ++mountCount;
  const ds = host.dataset;
  const startPreset = ds.preset && ds.preset in PRESETS ? ds.preset : "full";
  const startNoise = ds.noise !== "0";
  const startZoom = ds.zoom === "1";
  const st = { noisy: startNoise, zoom: startZoom, samples: [] as Sample[], preset: startPreset as string | null };
  const setPreset = (key: string) => {
    const p = PRESETS[key];
    st.samples = p.k.map((k, i) => ({ k, t: p.t[i] }));
    st.preset = key;
  };
  setPreset(startPreset);

  // ---- controls ----
  const row = htmlEl("div", "dz-row", "", host);
  htmlEl("span", "", "design:", row);
  const pBtns = presets(row, PRESET_KEYS.map(key => ({ label: PRESETS[key].label, apply: () => { setPreset(key); rebuildHandles(); draw(); } })));
  const nBox = checkbox(row, "noise (the course data)", st.noisy, v => { st.noisy = v; draw(); });
  const zBox = checkbox(row, "zoom (b)", st.zoom, v => { st.zoom = v; draw(); });
  button(row, "reset", () => {
    st.noisy = startNoise; st.zoom = startZoom; nBox.checked = st.noisy; zBox.checked = st.zoom;
    setPreset(startPreset); rebuildHandles(); draw();
  });

  // ---- main SVG: three panels ----
  const W = 1150, H = 340, top = 38, bot = H - 52;
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, class: "dz-main", role: "img", "aria-label": "PK samples, confidence ellipse and profile of ka" }, host);
  const defs = svgEl("defs", {}, svg);
  const gStatic = svgEl("g", {}, svg);
  const gHandles = svgEl("g", {}, svg);
  // panel extents
  const A = { x0: 62, x1: 400 }, B = { x0: 500, x1: 770 }, C = { x0: 852, x1: 1140 };
  const clipB = `dz-clipB-${uid}`, clipC = `dz-clipC-${uid}`;
  const cpB = svgEl("clipPath", { id: clipB }, defs);
  svgEl("rect", { x: B.x0, y: top, width: B.x1 - B.x0, height: bot - top }, cpB);
  const cpC = svgEl("clipPath", { id: clipC }, defs);
  svgEl("rect", { x: C.x0, y: top, width: C.x1 - C.x0, height: bot - top }, cpC);

  const ax = scale(0, T_MAX, A.x0, A.x1), ay = scale(0, 4, bot, top);
  const cx = scale(Math.log(KA0), Math.log(KA1), C.x0, C.x1), cy = scale(0, 16, bot, top);

  // ---- bottom: singular values + tables ----
  const bottom = htmlEl("div", "dz-bottom", "", host);
  const SW = 380, SH = 190;
  const svgS = svgEl("svg", { viewBox: `0 0 ${SW} ${SH}`, class: "dz-sv", role: "img", "aria-label": "Singular values of the Jacobian" }, bottom);
  const tabBox = htmlEl("div", "", "", bottom);
  const info = htmlEl("div", "dz-info", "", bottom);

  // ---- handles (persistent, so a drag survives the redraws) ----
  let handles: { g: SVGGElement; dot: SVGCircleElement; hit: SVGCircleElement }[] = [];
  function rebuildHandles(): void {
    gHandles.replaceChildren();
    handles = st.samples.map((_, i) => {
      const g = svgEl("g", {}, gHandles);
      const hit = svgEl("circle", { r: 17, fill: "transparent" }, g);
      const dot = svgEl("circle", { r: 8.5, fill: COL.ink, stroke: "#fff", "stroke-width": 2 }, g);
      draggable(g as unknown as SVGElement, ev => svgPoint(svg, ev), x => {
        st.samples[i].t = Math.round(Math.min(T_MAX, Math.max(T_MIN, ax.invert(x))) * 20) / 20;
        st.preset = null;
        schedule();
      });
      return { g, dot, hit };
    });
  }
  let pending = false;
  const schedule = () => { if (pending) return; pending = true; requestAnimationFrame(() => { pending = false; draw(); }); };

  // ghost (full design) per noise mode
  const ghostCache: Record<string, Result> = {};
  const ghost = (noisy: boolean): Result => {
    const key = noisy ? "n" : "c";
    if (!ghostCache[key]) ghostCache[key] = compute(PRESETS.full.k.map((k, i) => ({ k, t: PRESETS.full.t[i] })), noisy);
    return ghostCache[key];
  };

  const txt = (x: number, y: number, s: string, attrs: Record<string, string | number> = {}, parent: Element = gStatic) => {
    const e = svgEl("text", { x, y, ...attrs }, parent); e.textContent = s; return e;
  };

  function ellipsePath(r: Result, X: (v: number) => number, Y: (v: number) => number): string {
    const a = r.a;
    if (a.singular) return "";
    const c11 = a.cov[0][0], c12 = a.cov[0][1], c22 = a.cov[1][1];
    const l11 = Math.sqrt(c11), l21 = c12 / l11, l22 = Math.sqrt(Math.max(0, c22 - l21 * l21));
    const k = Math.sqrt(CRIT);
    let d = "";
    for (let i = 0; i <= 120; i++) {
      const th = 2 * Math.PI * i / 120, u = Math.cos(th), v = Math.sin(th);
      const ka = a.th[0] + k * l11 * u, ke = a.th[1] + k * (l21 * u + l22 * v);
      const px = Math.max(-1e4, Math.min(1e4, X(ka))), py = Math.max(-1e4, Math.min(1e4, Y(ke)));
      d += (i ? "L" : "M") + px.toFixed(1) + "," + py.toFixed(1);
    }
    return d;
  }

  function draw(): void {
    pBtns.forEach((b, i) => b.classList.toggle("active", PRESET_KEYS[i] === st.preset));
    const r = compute(st.samples, st.noisy);
    const gr = ghost(st.noisy);
    const a = r.a;
    gStatic.replaceChildren();

    // ---------- (a) data and curves ----------
    txt(A.x0 - 54, 24, "(a) drag the samples along t", { class: "ttl" });
    for (let v = 0; v <= 4; v++) {
      svgEl("line", { x1: A.x0, x2: A.x1, y1: ay(v), y2: ay(v), class: v === 0 ? "axis" : "grid" }, gStatic);
      txt(A.x0 - 8, ay(v) + 6, String(v), { "text-anchor": "end", class: "sm" });
    }
    for (let t = 0; t <= T_MAX; t += 3) {
      svgEl("line", { x1: ax(t), x2: ax(t), y1: bot, y2: bot + 5, class: "axis" }, gStatic);
      txt(ax(t), bot + 22, String(t), { "text-anchor": "middle", class: "sm" });
    }
    txt((A.x0 + A.x1) / 2, H - 6, "time t (h)", { "text-anchor": "middle" });
    txt(16, (top + bot) / 2, "C (mg/L)", { "text-anchor": "middle", transform: `rotate(-90 16 ${(top + bot) / 2})` });
    const curve = (ka: number, ke: number, V: number) => Array.from({ length: 151 }, (_, i) => i * T_MAX / 150)
      .map(t => `${ax(t).toFixed(1)},${ay(Math.min(4.15, pk(t, ka, ke, V))).toFixed(1)}`).join(" ");
    svgEl("polyline", { points: curve(TRUTH[0], TRUTH[1], TRUTH[2]), fill: "none", stroke: COL.grey, "stroke-width": 2.5, "stroke-dasharray": "7 5" }, gStatic);
    if (a.th.every(Number.isFinite))
      svgEl("polyline", { points: curve(a.th[0], a.th[1], a.th[2]), fill: "none", stroke: COL.blue, "stroke-width": 3 }, gStatic);
    // rug of sample times
    for (const s of st.samples) svgEl("line", { x1: ax(s.t), x2: ax(s.t), y1: bot - 12, y2: bot, stroke: COL.red, "stroke-width": 3 }, gStatic);
    // legend
    const lx = ax(8.6);
    svgEl("line", { x1: lx, x2: lx + 24, y1: ay(3.75), y2: ay(3.75), stroke: COL.grey, "stroke-width": 2.5, "stroke-dasharray": "7 5" }, gStatic);
    txt(lx + 30, ay(3.75) + 6, "truth");
    svgEl("line", { x1: lx, x2: lx + 24, y1: ay(3.35), y2: ay(3.35), stroke: COL.blue, "stroke-width": 3 }, gStatic);
    txt(lx + 30, ay(3.35) + 6, "fit");
    // handles
    st.samples.forEach((s, i) => {
      const h = handles[i];
      const y = ay(Math.max(-0.2, Math.min(4.2, yOf(s, st.noisy))));
      h.g.setAttribute("transform", `translate(${ax(s.t).toFixed(1)},${y.toFixed(1)})`);
    });

    // ---------- (b) Delta S over (ka, ke), Wald ellipse ----------
    const kaR = st.zoom ? [0.5, 1.7] : [0, 3], keR = st.zoom ? [0.12, 0.28] : [0, 0.4];
    const bx = scale(kaR[0], kaR[1], B.x0, B.x1), by = scale(keR[0], keR[1], bot, top);
    txt(B.x0 - 52, 24, "(b) ΔS/s², V re-fitted", { class: "ttl" });
    const kaT = st.zoom ? [0.6, 1, 1.4] : [0, 1, 2, 3], keT = st.zoom ? [0.15, 0.2, 0.25] : [0, 0.1, 0.2, 0.3, 0.4];
    for (const v of keT) {
      svgEl("line", { x1: B.x0, x2: B.x1, y1: by(v), y2: by(v), class: "grid" }, gStatic);
      txt(B.x0 - 6, by(v) + 6, String(v), { "text-anchor": "end", class: "sm" });
    }
    for (const v of kaT) {
      svgEl("line", { x1: bx(v), x2: bx(v), y1: top, y2: bot, class: "grid" }, gStatic);
      txt(bx(v), bot + 22, String(v), { "text-anchor": "middle", class: "sm" });
    }
    svgEl("rect", { x: B.x0, y: top, width: B.x1 - B.x0, height: bot - top, fill: "none", class: "axis" }, gStatic);
    txt((B.x0 + B.x1) / 2, H - 6, "kₐ (1/h)", { "text-anchor": "middle" });
    txt(B.x0 - 46, (top + bot) / 2, "kₑ (1/h)", { "text-anchor": "middle", transform: `rotate(-90 ${B.x0 - 46} ${(top + bot) / 2})` });
    const gB = svgEl("g", { "clip-path": `url(#${clipB})` }, gStatic);
    const NX = 110, NY = 80;
    const gx = (i: number) => kaR[0] + (kaR[1] - kaR[0]) * i / (NX - 1), gy = (j: number) => keR[0] + (keR[1] - keR[0]) * j / (NY - 1);
    const z: number[][] = [];
    for (let i = 0; i < NX; i++) {
      const col: number[] = [];
      for (let j = 0; j < NY; j++) {
        const ka = Math.max(1e-4, gx(i)), ke = Math.max(1e-4, gy(j));
        col.push((sseKaKe(r.ts, r.ys, ka, ke) - r.Smin) / a.s2);
      }
      z.push(col);
    }
    const X = (i: number) => bx(gx(0) + (gx(1) - gx(0)) * i), Y = (j: number) => by(gy(0) + (gy(1) - gy(0)) * j);
    for (const [lev, w] of [[100, 1.2], [30, 1.2], [10, 1.2], [CRIT, 3]] as [number, number][])
      svgEl("path", { d: contourPath(z, lev, X, Y), fill: "none", stroke: COL.blue, "stroke-width": w, "stroke-opacity": lev === CRIT ? 1 : 0.55 }, gB);
    if (st.preset !== "full")
      svgEl("path", { d: ellipsePath(gr, bx, by), fill: "none", stroke: COL.grey, "stroke-width": 2 }, gB);
    const ep = ellipsePath(r, bx, by);
    if (ep) svgEl("path", { d: ep, fill: COL.red, "fill-opacity": 0.08, stroke: COL.red, "stroke-width": 3, "stroke-dasharray": "8 4" }, gB);
    svgEl("circle", { cx: bx(TRUTH[0]), cy: by(TRUTH[1]), r: 5, fill: "#fff", stroke: COL.ink, "stroke-width": 2 }, gB);
    if (a.th.every(Number.isFinite)) svgEl("circle", { cx: bx(a.th[0]), cy: by(a.th[1]), r: 5.5, fill: COL.ink }, gB);

    // ---------- (c) profile of ka ----------
    txt(C.x0 - 52, 24, "(c) profile of kₐ (ΔS/s²)", { class: "ttl" });
    for (const v of [0, 3.84, 8, 12, 16]) {
      svgEl("line", { x1: C.x0, x2: C.x1, y1: cy(v), y2: cy(v), class: v === 0 ? "axis" : "grid" }, gStatic);
      if (v !== 3.84) txt(C.x0 - 6, cy(v) + 6, String(v), { "text-anchor": "end", class: "sm" });
    }
    for (const v of [0.1, 1, 10]) {
      svgEl("line", { x1: cx(Math.log(v)), x2: cx(Math.log(v)), y1: top, y2: bot, class: "grid" }, gStatic);
      txt(cx(Math.log(v)), bot + 22, String(v), { "text-anchor": "middle", class: "sm" });
    }
    txt((C.x0 + C.x1) / 2, H - 6, "kₐ (1/h), log scale", { "text-anchor": "middle" });
    svgEl("line", { x1: C.x0, x2: C.x1, y1: cy(CRIT), y2: cy(CRIT), stroke: COL.grey, "stroke-width": 2, "stroke-dasharray": "6 4" }, gStatic);
    txt(C.x1 - 4, cy(CRIT) - 6, "3.84", { "text-anchor": "end", class: "sm" });
    const gC = svgEl("g", { "clip-path": `url(#${clipC})` }, gStatic);
    const profLine = (g: number[], p: number[]) => g.map((k, i) => `${cx(Math.log(k)).toFixed(1)},${cy(Math.max(-0.5, Math.min(30, p[i]))).toFixed(1)}`).join(" ");
    if (st.preset !== "full") svgEl("polyline", { points: profLine(gr.grid, gr.prof), fill: "none", stroke: COL.grey, "stroke-width": 2 }, gC);
    if (!a.singular) {
      const wald = KA_GRID.map(k => ((k - a.th[0]) / a.se[0]) ** 2);
      svgEl("polyline", { points: profLine(KA_GRID, wald), fill: "none", stroke: COL.red, "stroke-width": 2.5, "stroke-dasharray": "3 4" }, gC);
    }
    svgEl("polyline", { points: profLine(r.grid, r.prof), fill: "none", stroke: COL.blue, "stroke-width": 3.5 }, gC);
    if (r.lo !== null || r.hi !== null) {
      const x0 = r.lo !== null ? cx(Math.log(r.lo)) : C.x0, x1 = r.hi !== null ? cx(Math.log(r.hi)) : C.x1;
      svgEl("line", { x1: x0, x2: x1, y1: cy(0) - 4, y2: cy(0) - 4, stroke: COL.blue, "stroke-width": 6, "stroke-opacity": 0.5 }, gC);
    }

    // ---------- singular values ----------
    svgS.replaceChildren();
    const sx0 = 70, sx1 = SW - 10, sTop = 30, sBot = SH - 40;
    const sy = scale(-3, 2, sBot, sTop);
    const st2 = svgEl("text", { x: 4, y: 20, class: "ttl" }, svgS); st2.textContent = "(d) singular values of J (log θ)";
    for (let e = -3; e <= 2; e++) {
      svgEl("line", { x1: sx0, x2: sx1, y1: sy(e), y2: sy(e), class: e === -3 ? "axis" : "grid" }, svgS);
      const t = svgEl("text", { x: sx0 - 6, y: sy(e) + 6, "text-anchor": "end", class: "sm" }, svgS);
      t.textContent = e === 0 ? "1" : e === 1 ? "10" : e === 2 ? "100" : `1e${e}`;
    }
    const sNoise = Math.sqrt(a.s2);
    const bw = (sx1 - sx0) / 3;
    a.sv.forEach((v, i) => {
      const lv = Math.max(-3, Math.log10(Math.max(v, 1e-300)));
      const x = sx0 + bw * i + bw * 0.2;
      const below = v < sNoise;
      svgEl("rect", { x, y: sy(lv), width: bw * 0.6, height: Math.max(0, sBot - sy(lv)), fill: below ? COL.orange : COL.blue, "fill-opacity": 0.85 }, svgS);
      const t = svgEl("text", { x: x + bw * 0.3, y: Math.max(sTop + 14, sy(lv) - 6), "text-anchor": "middle" }, svgS);
      t.textContent = v < 1e-3 ? "≈0" : fmt(v, 3);
      const t2 = svgEl("text", { x: x + bw * 0.3, y: sBot + 20, "text-anchor": "middle", class: "sm" }, svgS);
      t2.textContent = `σ${"₁₂₃"[i]}`;
    });
    svgEl("line", { x1: sx0, x2: sx1, y1: sy(Math.log10(sNoise)), y2: sy(Math.log10(sNoise)), stroke: COL.red, "stroke-width": 2, "stroke-dasharray": "6 4" }, svgS);
    const tn = svgEl("text", { x: sx1, y: SH - 4, "text-anchor": "end", class: "sm" }, svgS);
    tn.textContent = `dashed: ${a.sEstimated ? "s" : "σ"} = ${sNoise.toFixed(3)}`;
    const tk = svgEl("text", { x: 4, y: SH - 4, class: "sm" }, svgS);
    tk.textContent = a.singular ? "κ = ∞" : `κ = ${fmt(a.sv[0] / a.sv[2], 3)}`;

    // ---------- tables ----------
    const names = ["kₐ", "kₑ", "V"];
    const truth = TRUTH;
    let h = `<table class="dz-t"><tr><th></th><th>truth</th><th>estimate</th><th>SE</th><th>SE/est.</th></tr>`;
    a.th.forEach((v, i) => {
      const rel = a.se[i] / v;
      h += `<tr><td>${names[i]}</td><td>${truth[i]}</td><td>${fmt(v, 4)}</td><td>${fmt(a.se[i], 3)}</td><td>${Number.isFinite(rel) && rel < 100 ? (100 * rel).toPrecision(2) + " %" : "∞"}</td></tr>`;
    });
    h += `</table><table class="dz-t" style="margin-top:.3em"><tr><th>corr.</th><th>kₐ,kₑ</th><th>kₐ,V</th><th>kₑ,V</th></tr>`;
    const cr = (i: number, j: number) => Number.isFinite(a.corr[i][j]) ? a.corr[i][j].toFixed(2) : "–";
    h += `<tr><td></td><td>${cr(0, 1)}</td><td>${cr(0, 2)}</td><td>${cr(1, 2)}</td></tr></table>`;
    tabBox.innerHTML = h;

    // ---------- info ----------
    info.replaceChildren();
    const line = (html: string, cls = "") => { const d = htmlEl("div", cls, "", info); d.innerHTML = html; };
    line(`<b>${st.samples.length} samples</b>, first at ${Math.min(...r.ts).toFixed(2)} h. SSE ${fmt(a.S, 4)}; ${a.sEstimated ? `s = √(SSE/(n−3)) = ${Math.sqrt(a.s2).toFixed(4)}` : "noise-free data: SEs use σ = 0.1"}`);
    const iv = (lo: number | null, hi: number | null) => `(${lo === null ? "&lt;" + KA0 : lo.toFixed(2)}, ${hi === null ? "∞" : hi.toFixed(2)})`;
    line(`<span class="dz-key" style="border-color:${COL.blue}"></span>95 % profile, ΔS/s² ≤ 3.84: kₐ ∈ ${a.th[0] > KA0 && a.th[0] < KA1 ? iv(r.lo, r.hi) : "–"}`);
    line(`<span class="dz-key" style="border-color:${COL.red};border-top-style:dotted"></span>95 % Wald, ±1.96 SE: kₐ ∈ ${a.singular ? "–" : `(${(a.th[0] - 1.96 * a.se[0]).toFixed(2)}, ${(a.th[0] + 1.96 * a.se[0]).toFixed(2)})`}`);
    if (!a.singular) {
      const v3 = a.V.map(row => row[2]);
      const sg = v3[0] < 0 ? -1 : 1;
      line(`weakest direction v₃ (log kₐ, log kₑ, log V) = (${v3.map(v => (sg * v).toFixed(2)).join(", ")}); SE along v₃ = ${fmtSig(Math.sqrt(a.s2) / a.sv[2])}`);
    }
    line(`${st.preset !== "full" ? `<span class="dz-key" style="border-color:${COL.grey}"></span>grey: all 12 samples. ` : ""}○ truth, ● estimate`);
    if (a.singular) line("J is singular at the fit: no standard errors. The fit runs off along a flat direction.", "dz-warn");
  }

  rebuildHandles();
  draw();
}

autoMount(NAME, mount);

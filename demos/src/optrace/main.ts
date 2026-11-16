// Optimiser race for the gradient-methods unit. A 2-D objective is drawn as a contour map;
// click it to choose a start. Steepest descent (Armijo backtracking), Newton (same line
// search, H + mu I if H is not positive definite), BFGS (B0 = I) and Nelder-Mead (R optim's
// rules, simplex drawn) all run from that start; "race" animates them iteration by
// iteration. The right panel shows f - f* against the iteration on a log axis (straight
// line = linear convergence) and the counts of f, gradient and Hessian calls.
// Presets: an elongated quadratic with a condition-number slider; the Rosenbrock banana;
// the lecture's 4PL SSE in (EC50, hill) with bottom and top fixed at the best fit, with
// EC50 in uM, in nM (bad scaling) and as log10 EC50.
// Auto-mounts on <div data-demo="optrace">. Optional attributes:
//   data-preset="quad|rosen|ec50|nm|log"  data-kappa="10"  data-start="x y"
//   data-axis="short|full"  (iteration axis of the convergence plot)
import { svgEl, htmlEl, slider, button, presets, scale, svgPoint, numList, autoMount } from "../lib/ui";
import {
  Problem, Trace, Vec, quadratic, ROSEN, FPL_UM, FPL_NM, FPL_LOG, polish, eig2,
  steepest, newton, bfgs, nelderMead,
} from "./math";

const NAME = "optrace";
// Okabe-Ito, as in make-gradient-methods-figs.R (fig04: SD red, Newton blue, BFGS orange)
const METHODS = [
  { id: "sd", label: "steepest descent", col: "#D55E00", run: steepest },
  { id: "newton", label: "Newton", col: "#0072B2", run: newton },
  { id: "bfgs", label: "BFGS", col: "#E69F00", run: bfgs },
  { id: "nm", label: "Nelder–Mead", col: "#009E73", run: nelderMead },
] as const;
const INK = "#222", GRID = "#e3e3e3", CONT = "#8c8c8c";

const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .or-row { display: flex; flex-wrap: wrap; gap: .35em 1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .or-group { display: flex; flex-wrap: wrap; gap: .35em; align-items: center; }
[data-demo="${NAME}"] .or-row .qb-slider { min-width: 260px; flex: 0 1 300px; grid-template-columns: auto 1fr 3.2em; }
[data-demo="${NAME}"] .or-main { display: flex; flex-wrap: wrap; gap: 12px; width: 100%; align-items: flex-start; }
[data-demo="${NAME}"] .or-map { flex: 1 1 400px; max-width: 560px; min-width: 0; }
[data-demo="${NAME}"] .or-right { flex: 1 1 440px; min-width: 0; display: flex; flex-direction: column; gap: .3em; }
[data-demo="${NAME}"] svg { width: 100%; height: auto; }
[data-demo="${NAME}"] .or-map svg { cursor: crosshair; touch-action: manipulation; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 19px; font-weight: 600; }
[data-demo="${NAME}"] table.or-tab { border-collapse: collapse; font-size: 17px; width: 100%; }
[data-demo="${NAME}"] .or-tab th, [data-demo="${NAME}"] .or-tab td { padding: .1em .45em; text-align: right; white-space: nowrap; border-bottom: 1px solid #e3e3e3; }
[data-demo="${NAME}"] .or-tab th { font-weight: 600; font-size: 16px; color: #444; }
[data-demo="${NAME}"] .or-tab td:first-child, [data-demo="${NAME}"] .or-tab th:first-child { text-align: left; }
[data-demo="${NAME}"] .or-tab td.st { text-align: left; font-size: 16px; }
[data-demo="${NAME}"] .or-sw { display: inline-block; width: 1.4em; height: 0; border-top: 4px solid; vertical-align: .3em; margin-right: .35em; }
[data-demo="${NAME}"] .or-warn { color: #b00020; }
[data-demo="${NAME}"] .or-info { font-size: 16px; color: #444; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 18px; margin: 0; gap: .25em; }
.reveal [data-demo="${NAME}"] .qb-btn { font-size: 16px; padding: .12em .5em; }
.reveal [data-demo="${NAME}"] .or-main { flex-wrap: nowrap; }
.reveal [data-demo="${NAME}"] .or-map { flex: 0 0 590px; max-width: 590px; }
.reveal [data-demo="${NAME}"] .or-tab th, .reveal [data-demo="${NAME}"] .or-tab td { padding: 0 .35em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- formatting ----
const fmtInt = (n: number): string => String(Math.round(n)).replace(/\B(?=(\d{3})+(?!\d))/g, " ");
function sci(v: number, d = 2): string {
  if (!Number.isFinite(v)) return String(v);
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e-3 && a < 1e5) return v.toPrecision(d + 1).replace(/\.?0+$/, m => (m.startsWith(".") ? "" : m));
  const [m, e] = v.toExponential(d - 1).split("e");
  return `${m}e${Number(e)}`;
}
const fmtCoord = (v: number): string => Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 10 ? v.toFixed(1) : v.toFixed(2);

function niceTicks(a: number, b: number, n = 5): number[] {
  const raw = (b - a) / n, p = 10 ** Math.floor(Math.log10(raw));
  const st = [1, 2, 2.5, 5, 10].map(m => m * p).find(s => s >= raw) ?? raw;
  const out: number[] = [];
  for (let v = Math.ceil(a / st - 1e-9) * st; v <= b + 1e-9 * st; v += st) out.push(Math.abs(v) < 1e-12 * st ? 0 : v);
  return out;
}
const tickLab = (v: number): string => {
  const s = Number(v.toPrecision(6));
  return String(s).replace("-", "−");
};

// ---- presets ----
type PresetId = "quad" | "rosen" | "ec50" | "nm" | "log";
const PRESET_IDS: PresetId[] = ["quad", "rosen", "ec50", "nm", "log"];
function problemFor(id: PresetId, kappa: number): Problem {
  switch (id) {
    case "quad": return quadratic(kappa);
    case "rosen": return ROSEN;
    case "ec50": return FPL_UM;
    case "nm": return FPL_NM;
    default: return FPL_LOG;
  }
}
const PRESET_LABEL: Record<PresetId, string> = {
  quad: "bowl (κ slider)", rosen: "Rosenbrock", ec50: "4PL, EC50 in µM", nm: "4PL, EC50 in nM", log: "4PL, log₁₀ EC50",
};

// ---- contour lines by marching squares ----
function contourPath(Z: Float64Array, nx: number, ny: number, lev: number,
  px: (i: number) => number, py: (j: number) => number): string {
  let d = "";
  const at = (i: number, j: number) => Z[j * nx + i];
  for (let j = 0; j < ny - 1; j++) for (let i = 0; i < nx - 1; i++) {
    const v = [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)];
    if (v.some(z => !Number.isFinite(z))) continue;
    const c = (v[0] > lev ? 1 : 0) | (v[1] > lev ? 2 : 0) | (v[2] > lev ? 4 : 0) | (v[3] > lev ? 8 : 0);
    if (c === 0 || c === 15) continue;
    // edge points: 0 bottom (0-1), 1 right (1-2), 2 top (3-2), 3 left (0-3)
    const e = (k: number): [number, number] => {
      const t = (a: number, b: number) => (lev - a) / (b - a);
      switch (k) {
        case 0: return [px(i + t(v[0], v[1])), py(j)];
        case 1: return [px(i + 1), py(j + t(v[1], v[2]))];
        case 2: return [px(i + t(v[3], v[2])), py(j + 1)];
        default: return [px(i), py(j + t(v[0], v[3]))];
      }
    };
    const segs: [number, number][] = {
      1: [[3, 0]], 2: [[0, 1]], 3: [[3, 1]], 4: [[1, 2]], 5: [[3, 2], [0, 1]], 6: [[0, 2]], 7: [[3, 2]],
      8: [[2, 3]], 9: [[0, 2]], 10: [[0, 3], [1, 2]], 11: [[1, 2]], 12: [[1, 3]], 13: [[0, 1]], 14: [[0, 3]],
    }[c] as [number, number][];
    for (const [a, b] of segs) {
      const p = e(a), q = e(b);
      d += `M${p[0].toFixed(1)} ${p[1].toFixed(1)}L${q[0].toFixed(1)} ${q[1].toFixed(1)}`;
    }
  }
  return d;
}

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const uid = ++mountCount;
  const ds = host.dataset;
  const pre0: PresetId = PRESET_IDS.includes(ds.preset as PresetId) ? ds.preset as PresetId : "quad";
  const kap0 = Number(ds.kappa) >= 1 && Number(ds.kappa) <= 1000 ? Number(ds.kappa) : 10;
  const st0 = numList(ds.start, []);
  const st = {
    preset: pre0, lk: Math.log10(kap0),
    start: (st0.length === 2 ? [st0[0], st0[1]] : null) as Vec | null,
    full: ds.axis === "full",
    k: 0, playing: false,
  };

  // ---- controls ----
  const row1 = htmlEl("div", "or-row", "", host);
  const pg = htmlEl("div", "or-group", "", row1);
  const pBtns = presets(pg, PRESET_IDS.map(id => ({
    label: PRESET_LABEL[id], apply: () => { st.preset = id; st.start = null; recompute(false); },
  })));
  const kS = slider(row1, "κ", 0, 3, 0.1, st.lk, v => { st.lk = v; recompute(false); }, v => fmtKappa(10 ** v));
  const row2 = htmlEl("div", "or-row", "", host);
  const rg = htmlEl("div", "or-group", "", row2);
  const playBtn = button(rg, "▶ race", () => { if (st.playing) stop(); else play(); });
  button(rg, "+1 iteration", () => { stop(); st.k = Math.min(st.k + 1, kEnd); drawDyn(); });
  button(rg, "show all", () => { stop(); st.k = kEnd; drawDyn(); });
  button(rg, "default start", () => { st.start = null; recompute(true); });
  const ag = htmlEl("div", "or-group", "", row2);
  htmlEl("span", "", "iteration axis:", ag);
  const aBtns = presets(ag, [
    { label: "short", apply: () => { st.full = false; drawDyn(); } },
    { label: "all", apply: () => { st.full = true; drawDyn(); } },
  ]);

  // ---- layout ----
  const main = htmlEl("div", "or-main", "", host);
  const mapBox = htmlEl("div", "or-map", "", main);
  const right = htmlEl("div", "or-right", "", main);
  const W = 540, H = 470;
  const L = 82, R = W - 26, T = 34, B = H - 66;
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Contour map with optimiser paths" }, mapBox);
  const defs = svgEl("defs", {}, svg);
  const clip = svgEl("clipPath", { id: `or-clip-${uid}` }, defs);
  svgEl("rect", { x: L, y: T, width: R - L, height: B - T }, clip);
  const gStatic = svgEl("g", {}, svg);
  const gDyn = svgEl("g", { "clip-path": `url(#or-clip-${uid})` }, svg);
  const gTop = svgEl("g", {}, svg);

  const CW = 600, CH = 276, cL = 84, cR = CW - 30, cT = 46, cB = CH - 56;
  const csvg = svgEl("svg", { viewBox: `0 0 ${CW} ${CH}`, role: "img", "aria-label": "f minus f* against iteration, log scale" }, right);
  const tab = htmlEl("table", "or-tab", "", right);
  const info = htmlEl("div", "or-info", "", right);

  // ---- state computed per landscape ----
  let prob: Problem = problemFor(st.preset, 10 ** st.lk);
  let fstar = 0, xstar: Vec = [0, 0], kappaStar = 1;
  let traces: Trace[] = [];
  let kEnd = 0, shortMax = 50, fullMax = 50, yBot = -16, yTop = 4;
  let sx = scale(0, 1, L, R), sy = scale(0, 1, B, T);
  let timer = 0;

  function fmtKappa(k: number): string { return k >= 10 ? k.toFixed(0) : k.toPrecision(2); }

  function stop(): void {
    st.playing = false;
    if (timer) { window.clearTimeout(timer); timer = 0; }
    playBtn.textContent = "▶ race";
  }
  function play(): void {
    if (st.k >= kEnd) st.k = 0;
    st.playing = true;
    playBtn.textContent = "❚❚ pause";
    const tick = () => {
      if (!host.isConnected) { stop(); return; }
      st.k = Math.min(kEnd, st.k + (st.k < 25 ? 1 : Math.ceil(st.k * 0.05)));
      drawDyn();
      if (st.k >= kEnd) { stop(); return; }
      timer = window.setTimeout(tick, st.k < 25 ? 220 : 90);
    };
    drawDyn();
    timer = window.setTimeout(tick, 300);
  }

  function recompute(race: boolean): void {
    stop();
    prob = problemFor(st.preset, 10 ** st.lk);
    const ps = polish(prob);
    xstar = ps.x; fstar = ps.f;
    const ev = eig2(prob.H(xstar));
    kappaStar = ev[1] / ev[0];
    const x0: Vec = st.start ?? prob.start;
    traces = METHODS.map(m => m.run(prob, x0));
    for (const t of traces) for (const f of t.fx) if (f < fstar) fstar = f;
    kEnd = Math.max(...traces.map(t => t.x.length - 1));
    const fast = Math.max(...traces.slice(1).map(t => t.x.length - 1));
    const nice = [10, 20, 30, 50, 100, 200, 300, 500, 1000, 2000, 5000, 10000, 20000];
    shortMax = nice.find(v => v >= fast * 1.08) ?? fast;
    fullMax = nice.find(v => v >= kEnd) ?? kEnd;
    const d0 = Math.max(...traces.map(t => t.fx[0] - fstar));
    yTop = Math.ceil(Math.log10(Math.max(d0, 1e-12)));
    yBot = fstar !== 0 ? Math.floor(Math.log10(1e-16 * Math.abs(fstar))) : -20;
    st.k = race ? 0 : kEnd;
    drawStatic();
    if (race) play(); else drawDyn();
  }

  // ---- static layer: heat map, contours, axes ----
  function drawStatic(): void {
    pBtns.forEach((b, i) => b.classList.toggle("active", PRESET_IDS[i] === st.preset));
    kS.input.closest("label")!.style.visibility = st.preset === "quad" ? "visible" : "hidden";
    const [x0, x1, y0, y1] = prob.dom;
    sx = scale(x0, x1, L, R); sy = scale(y0, y1, B, T);
    gStatic.replaceChildren(); gTop.replaceChildren();
    const nx = 181, ny = 151;
    const Z = new Float64Array(nx * ny);
    const ds: number[] = [];
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const v = prob.f([x0 + (x1 - x0) * i / (nx - 1), y0 + (y1 - y0) * j / (ny - 1)]);
      Z[j * nx + i] = v;
      if (Number.isFinite(v)) ds.push(Math.max(v - fstar, 0));
    }
    ds.sort((a, b) => a - b);
    const q = (p: number) => ds[Math.min(ds.length - 1, Math.floor(p * ds.length))];
    const lo = Math.log10(Math.max(q(0.01), 1e-12 * (1 + q(0.95)))), hi = Math.log10(Math.max(q(0.95), 1e-12));
    // heat map, drawn on a canvas and placed as an image
    const cv = document.createElement("canvas");
    cv.width = nx; cv.height = ny;
    const cx = cv.getContext("2d");
    if (cx) {
      const img = cx.createImageData(nx, ny);
      const top = Math.log10(Math.max(ds[ds.length - 1] ?? 1, 1e-12));
      for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
        const v = Z[j * nx + i], o = ((ny - 1 - j) * nx + i) * 4;
        let r = 226, g = 226, b = 226;
        if (Number.isFinite(v)) {
          const t = Math.min(1, Math.max(0, (Math.log10(Math.max(v - fstar, 1e-300)) - lo) / Math.max(top - lo, 1e-9)));
          r = 186 + (252 - 186) * t; g = 212 + (251 - 212) * t; b = 236 + (246 - 236) * t;
        }
        img.data[o] = r; img.data[o + 1] = g; img.data[o + 2] = b; img.data[o + 3] = 255;
      }
      cx.putImageData(img, 0, 0);
      svgEl("image", {
        href: cv.toDataURL(), x: L, y: T, width: R - L, height: B - T, preserveAspectRatio: "none",
        style: "image-rendering: auto",
      }, gStatic);
    }
    const px = (i: number) => L + (R - L) * i / (nx - 1), py = (j: number) => B - (B - T) * j / (ny - 1);
    const NL = 13;
    let dpath = "";
    for (let l = 0; l < NL; l++) dpath += contourPath(Z, nx, ny, fstar + 10 ** (lo + (hi - lo) * (l + 0.5) / NL), px, py);
    svgEl("path", { d: dpath, fill: "none", stroke: CONT, "stroke-width": 1.1 }, gStatic);
    svgEl("rect", { x: L, y: T, width: R - L, height: B - T, fill: "none", stroke: "#444" }, gStatic);
    for (const v of niceTicks(x0, x1, 5)) {
      svgEl("line", { x1: sx(v), x2: sx(v), y1: B, y2: B + 6, stroke: "#444" }, gStatic);
      svgEl("text", { x: sx(v), y: B + 25, "text-anchor": "middle" }, gStatic).textContent = tickLab(v);
    }
    for (const v of niceTicks(y0, y1, 5)) {
      svgEl("line", { x1: L - 6, x2: L, y1: sy(v), y2: sy(v), stroke: "#444" }, gStatic);
      svgEl("text", { x: L - 10, y: sy(v) + 6, "text-anchor": "end" }, gStatic).textContent = tickLab(v);
    }
    svgEl("text", { x: (L + R) / 2, y: H - 12, "text-anchor": "middle" }, gStatic).textContent = prob.xlab;
    svgEl("text", { x: 20, y: (T + B) / 2, "text-anchor": "middle", transform: `rotate(-90 20 ${(T + B) / 2})` }, gStatic)
      .textContent = prob.ylab;
    svgEl("text", { x: 6, y: 22, class: "ttl" }, gStatic).textContent = prob.formula;
    // the minimum
    const mx = sx(xstar[0]), my = sy(xstar[1]);
    svgEl("path", { d: `M${mx - 8} ${my - 8}L${mx + 8} ${my + 8}M${mx - 8} ${my + 8}L${mx + 8} ${my - 8}`, stroke: INK, "stroke-width": 3 }, gTop);
  }

  // ---- dynamic layer: paths, convergence plot, table ----
  const upto = (t: Trace) => Math.min(st.k, t.x.length - 1);

  function drawDyn(): void {
    aBtns.forEach((b, i) => b.classList.toggle("active", (i === 1) === st.full));
    gDyn.replaceChildren();
    const x0: Vec = traces[0].x[0];
    traces.forEach((t, m) => {
      const n = upto(t), col = METHODS[m].col;
      if (t.simplex) {
        const S = t.simplex[n];
        svgEl("polygon", {
          points: S.map(p => `${sx(p[0])},${sy(p[1])}`).join(" "),
          fill: col, "fill-opacity": 0.18, stroke: col, "stroke-width": 2,
        }, gDyn);
      }
      const pts = t.x.slice(0, n + 1).map(p => `${sx(p[0]).toFixed(1)},${sy(p[1]).toFixed(1)}`).join(" ");
      svgEl("polyline", {
        points: pts, fill: "none", stroke: col, "stroke-width": m === 0 ? 2 : 2.6,
        "stroke-dasharray": t.simplex ? "6 4" : "none", "stroke-linejoin": "round",
      }, gDyn);
      if (n <= 400 && !t.simplex) {
        for (let i = 1; i <= n; i++) svgEl("circle", { cx: sx(t.x[i][0]), cy: sy(t.x[i][1]), r: 3.2, fill: col }, gDyn);
      }
      const c = t.x[n];
      svgEl("circle", { cx: sx(c[0]), cy: sy(c[1]), r: 7, fill: col, stroke: "#fff", "stroke-width": 2 }, gDyn);
    });
    svgEl("rect", { x: sx(x0[0]) - 7, y: sy(x0[1]) - 7, width: 14, height: 14, fill: "#fff", stroke: INK, "stroke-width": 2.5 }, gDyn);
    drawConv();
    drawTable();
  }

  function drawConv(): void {
    csvg.replaceChildren();
    const xm = st.full ? fullMax : shortMax;
    const cx = scale(0, xm, cL, cR), cy = scale(yBot, yTop, cB, cT);
    svgEl("text", { x: 6, y: 22, class: "ttl" }, csvg).textContent = "f − f*  (log scale) against iteration";
    const dec = yTop - yBot, stepD = dec > 16 ? 4 : dec > 8 ? 2 : 1;
    for (let e = Math.ceil(yBot / stepD) * stepD; e <= yTop; e += stepD) {
      svgEl("line", { x1: cL, x2: cR, y1: cy(e), y2: cy(e), stroke: GRID }, csvg);
      svgEl("text", { x: cL - 8, y: cy(e) + 6, "text-anchor": "end" }, csvg).textContent = e === 0 ? "1" : `1e${e}`;
    }
    for (const v of niceTicks(0, xm, 5)) {
      svgEl("line", { x1: cx(v), x2: cx(v), y1: cB, y2: cB + 5, stroke: "#444" }, csvg);
      svgEl("text", { x: cx(v), y: cB + 23, "text-anchor": "middle" }, csvg).textContent = fmtInt(v);
    }
    svgEl("line", { x1: cL, x2: cR, y1: cB, y2: cB, stroke: "#444" }, csvg);
    svgEl("line", { x1: cL, x2: cL, y1: cT, y2: cB, stroke: "#444" }, csvg);
    svgEl("text", { x: (cL + cR) / 2, y: CH - 8, "text-anchor": "middle" }, csvg).textContent = "iteration k";
    if (st.k <= xm) svgEl("line", { x1: cx(st.k), x2: cx(st.k), y1: cT, y2: cB, stroke: "#999", "stroke-dasharray": "4 4" }, csvg);
    const cid = `or-cclip-${uid}`;
    const cd = svgEl("defs", {}, csvg);
    const cc = svgEl("clipPath", { id: cid }, cd);
    svgEl("rect", { x: cL, y: cT - 6, width: cR - cL + 6, height: cB - cT + 12 }, cc);
    const g = svgEl("g", { "clip-path": `url(#${cid})` }, csvg);
    const ly = (f: number) => cy(Math.max(yBot, Math.min(yTop + 0.3, Math.log10(Math.max(f - fstar, 1e-300)))));
    traces.forEach((t, m) => {
      const n = Math.min(upto(t), Math.ceil(xm) + 1), col = METHODS[m].col;
      const pts: string[] = [];
      for (let i = 0; i <= n; i++) pts.push(`${cx(i).toFixed(1)},${ly(t.fx[i]).toFixed(1)}`);
      svgEl("polyline", { points: pts.join(" "), fill: "none", stroke: col, "stroke-width": 2.6, "stroke-linejoin": "round" }, g);
      if (n <= 60) for (let i = 0; i <= n; i++) {
        const X = cx(i), Y = ly(t.fx[i]);
        if (m === 1) svgEl("circle", { cx: X, cy: Y, r: 4.5, fill: col }, g);
        else if (m === 2) svgEl("path", { d: `M${X} ${Y - 6}L${X + 5.5} ${Y + 4}L${X - 5.5} ${Y + 4}Z`, fill: col }, g);
        else if (m === 3) svgEl("rect", { x: X - 3.5, y: Y - 3.5, width: 7, height: 7, fill: col }, g);
      }
    });
  }

  function drawTable(): void {
    tab.replaceChildren();
    const hr = htmlEl("tr", "", "", tab);
    for (const h of ["method", "iter", "f calls", "∇f calls", "H calls", "f − f*", ""]) htmlEl("th", "", h, hr);
    traces.forEach((t, m) => {
      const n = upto(t), done = n === t.x.length - 1;
      const tr = htmlEl("tr", "", "", tab);
      const td0 = htmlEl("td", "", "", tr);
      const sw = htmlEl("span", "or-sw", "", td0);
      sw.style.borderTopColor = METHODS[m].col;
      if (t.simplex) sw.style.borderTopStyle = "dashed";
      td0.append(METHODS[m].label);
      htmlEl("td", "", fmtInt(n), tr);
      htmlEl("td", "", fmtInt(t.nf[n]), tr);
      htmlEl("td", "", m === 3 ? "–" : fmtInt(t.ng[n]), tr);
      htmlEl("td", "", m === 1 ? fmtInt(t.nh[n]) : "–", tr);
      htmlEl("td", "", sci(Math.max(t.fx[n] - fstar, 0), 2), tr);
      const s = !done ? "…" : t.status === "converged" ? "✓" : t.status;
      htmlEl("td", "st" + (done && t.status !== "converged" ? " or-warn" : ""), s, tr);
    });
    const x0 = traces[0].x[0];
    const nt = traces[1], ns = nt.shifted ? nt.shifted[upto(nt)] : 0;
    const shift = ns > 0 ? ` Newton used H + μI in ${ns} step${ns > 1 ? "s" : ""} (H not positive definite).` : "";
    const unit = prob.hessFD ? " Bottom and top fixed at the best fit; H = central differences of ∇f." : "";
    info.textContent = `start (${fmtCoord(x0[0])}, ${fmtCoord(x0[1])}) · κ(H) at the minimum = ${sci(kappaStar, 3)} · f* = ${Number(fstar.toPrecision(8))} · iteration ${fmtInt(st.k)} of ${fmtInt(kEnd)}. Click the map to choose a start. ✓ = converged.${shift}${unit}`;
  }

  // ---- click to set the start ----
  svg.addEventListener("click", ev => {
    const [px, py] = svgPoint(svg, ev);
    if (px < L || px > R || py < T || py > B) return;
    st.start = [sx.invert(px), sy.invert(py)];
    recompute(true);
  });

  if (st.start && !(st.start[0] >= prob.dom[0] && st.start[0] <= prob.dom[1])) st.start = null;
  recompute(false);
}

autoMount(NAME, mount);

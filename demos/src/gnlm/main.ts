// Gauss-Newton vs Levenberg-Marquardt step by step, for the nonlinear-least-squares unit.
// The lecture's oral one-compartment PK fit (Bateman function, D = 100 mg, 12 samples):
//  (a) the RSS landscape over (ka, ke) on log axes, V at its best value for each pair
//      ("profiled", as in fig02), with the fit (V = 20.2 L) and its flip-flop mirror
//      (V = 3.8 L) as stars, the line ka = ke, a draggable start point, and the paths of
//      Gauss-Newton (vermillion) and Levenberg-Marquardt (green) projected on (ka, ke);
//      rejected LM trial points are small crosses;
//  (b) the data and the current model curves of both methods;
//  (c) a readout per method: pass, (ka, ke, V), RSS, lambda, last step, status.
// The methods are the lecture's R code: forward-difference Jacobian, QR solve, LM with
// Marquardt scaling, lambda0 = 0.01, factor 10. They move all three parameters; the
// start V is set by a slider. Parameters either as (ka, ke, V) or as their logarithms.
// Auto-mounts on <div data-demo="gnlm">. Optional attributes:
//   data-start="3 0.05 50"  data-param="lin|log"  data-method="both|gn|lm"
//   data-lambda="0.01"  data-speed="350" (ms per pass when running)
import { svgEl, setAttrs, htmlEl, slider, presets, button, scale, svgPoint, draggable, numList, autoMount } from "../lib/ui";
import { TT, Y, P_HAT, P_FLIP, rssProfile, bateman, newRun, step, toNat, Run, Param, Status } from "./math";

const NAME = "gnlm";
// Okabe-Ito, as in make-nonlinear-least-squares-figs.R (fig04: GN vermillion, LM green;
// fig02: the fit and its mirror as stars)
const COL = { gn: "#D55E00", lm: "#009E73", hat: "#0072B2", flip: "#E69F00", ref: "#737373", ink: "#222", pt: "#444" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .gl-main { display: flex; flex-wrap: wrap; gap: 12px; width: 100%; align-items: flex-start; }
[data-demo="${NAME}"] .gl-land { flex: 0 1 470px; min-width: 280px; height: auto; }
[data-demo="${NAME}"] .gl-right { flex: 1 1 400px; min-width: 0; display: flex; flex-direction: column; gap: .35em; }
[data-demo="${NAME}"] .gl-row { display: flex; flex-wrap: wrap; gap: .35em 1em; align-items: center; }
[data-demo="${NAME}"] .gl-group { display: flex; gap: .35em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .gl-row .qb-slider { flex: 1 1 260px; max-width: 360px; grid-template-columns: auto 1fr 3.4em; }
[data-demo="${NAME}"] .gl-data { width: 100%; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 19px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .gl-btn-main { font-weight: 600; }
[data-demo="${NAME}"] table.gl-tab { border-collapse: collapse; font-variant-numeric: tabular-nums; font-size: .86em; width: 100%; }
[data-demo="${NAME}"] table.gl-tab td, [data-demo="${NAME}"] table.gl-tab th { padding: .05em .45em; text-align: left; white-space: nowrap; }
[data-demo="${NAME}"] table.gl-tab td:first-child { color: #555; }
[data-demo="${NAME}"] table.gl-tab tr:nth-child(even) td { background: var(--qb-panel); }
[data-demo="${NAME}"] .gl-acc { color: ${COL.lm}; font-weight: 600; }
[data-demo="${NAME}"] .gl-rej { color: #b00020; font-weight: 600; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; gap: .25em; }
.reveal [data-demo="${NAME}"] .gl-land { flex: 0 0 455px; }
.reveal [data-demo="${NAME}"] table.gl-tab { font-size: 17px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- formatting ----
function num(v: number, d = 4): string {
  if (!Number.isFinite(v)) return Number.isNaN(v) ? "NaN" : (v > 0 ? "Inf" : "−Inf");
  const a = Math.abs(v);
  const s = a !== 0 && (a < 1e-3 || a >= 1e5) ? v.toExponential(Math.max(1, d - 2)).replace("e+", "e") : v.toPrecision(d);
  return s.replace("-", "−");
}
const rssFmt = (v: number): string => Number.isFinite(v) && v < 1e5 ? v.toPrecision(7) : num(v, 3);
const lamFmt = (v: number): string => { const k = Math.round(Math.log10(v)); return Math.abs(v / 10 ** k - 1) < 1e-6 ? (k === 0 ? "1" : `1e${k}`.replace("-", "−")) : num(v, 3); };

// ---- precomputed landscape (log10 ka, log10 ke) ----
const L0 = -2, L1 = 1.5, NG = 141;
const gridV = Array.from({ length: NG }, (_, i) => L0 + (L1 - L0) * i / (NG - 1));
let LAND: number[][] | null = null;          // LAND[iy][ix] = profiled RSS
function land(): number[][] {
  if (!LAND) LAND = gridV.map(le => gridV.map(la => rssProfile(10 ** la, 10 ** le)[0]));
  return LAND;
}
const LEVELS = [0.07, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40];

/** Marching squares: line segments of the level set Z = c, in grid-index coordinates. */
function contourSegs(Z: number[][], c: number): [number, number, number, number][] {
  const segs: [number, number, number, number][] = [];
  const n = Z.length;
  for (let j = 0; j < n - 1; j++) for (let i = 0; i < n - 1; i++) {
    const v = [Z[j][i], Z[j][i + 1], Z[j + 1][i + 1], Z[j + 1][i]];
    const pts: [number, number][] = [];
    const corners: [number, number][] = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]];
    for (let e = 0; e < 4; e++) {
      const a = v[e], b = v[(e + 1) % 4];
      if ((a < c) !== (b < c)) {
        const t = (c - a) / (b - a);
        const [x0, y0] = corners[e], [x1, y1] = corners[(e + 1) % 4];
        pts.push([x0 + t * (x1 - x0), y0 + t * (y1 - y0)]);
      }
    }
    if (pts.length === 2) segs.push([pts[0][0], pts[0][1], pts[1][0], pts[1][1]]);
    else if (pts.length === 4) { segs.push([...pts[0], ...pts[1]]); segs.push([...pts[2], ...pts[3]]); }
  }
  return segs;
}

/** Light (low RSS) to dark blue (high RSS), like hcl.colors("Blues 3", rev = TRUE). */
function blue(u: number): [number, number, number] {
  const stops: [number, number, number][] = [[255, 255, 255], [222, 235, 247], [178, 210, 232], [130, 186, 219], [87, 158, 205]];
  const x = Math.max(0, Math.min(1, u)) * (stops.length - 1);
  const k = Math.min(stops.length - 2, Math.floor(x)), f = x - k;
  return [0, 1, 2].map(c => Math.round(stops[k][c] + f * (stops[k + 1][c] - stops[k][c]))) as [number, number, number];
}

let heatURL = "";
function heatmap(): string {
  if (heatURL) return heatURL;
  const Z = land();
  const cv = document.createElement("canvas");
  cv.width = NG; cv.height = NG;
  const ctx = cv.getContext("2d");
  if (!ctx) return "";
  const img = ctx.createImageData(NG, NG);
  const lo = Math.log10(0.06), hi = Math.log10(65);
  for (let j = 0; j < NG; j++) for (let i = 0; i < NG; i++) {
    const [r, g, b] = blue((Math.log10(Z[j][i]) - lo) / (hi - lo));
    const k = 4 * ((NG - 1 - j) * NG + i);       // row 0 of the image = top = largest ke
    img.data[k] = r; img.data[k + 1] = g; img.data[k + 2] = b; img.data[k + 3] = 255;
  }
  ctx.putImageData(img, 0, 0);
  heatURL = cv.toDataURL();
  return heatURL;
}

type Method = "both" | "gn" | "lm";

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const uid = ++mountCount;
  const s0 = numList(ds.start, [3, 0.05, 50]);
  const startDefault = s0.length === 3 && s0.every(v => v > 0) ? s0 : [3, 0.05, 50];
  const lambda0 = Number(ds.lambda) > 0 ? Number(ds.lambda) : 0.01;
  const speed = Number(ds.speed) > 0 ? Number(ds.speed) : 350;
  const st = {
    start: startDefault.slice(),
    par: (ds.param === "log" ? "log" : "lin") as Param,
    method: (ds.method === "gn" || ds.method === "lm" ? ds.method : "both") as Method,
  };
  const init = { par: st.par, method: st.method };
  let runs: { gn: Run; lm: Run } = fresh();
  let timer: number | null = null;

  function fresh() {
    return { gn: newRun("gn", st.par, st.start), lm: newRun("lm", st.par, st.start, lambda0) };
  }
  const active = (): Run[] => st.method === "both" ? [runs.gn, runs.lm] : [runs[st.method]];

  // ---------------- layout ----------------
  const main = htmlEl("div", "gl-main", "", host);
  const WL = 470, HL = 480;
  const svgL = svgEl("svg", { viewBox: `0 0 ${WL} ${HL}`, class: "gl-land", role: "img",
    "aria-label": "RSS landscape over ka and ke with the paths of Gauss-Newton and Levenberg-Marquardt" }, main);
  const right = htmlEl("div", "gl-right", "", main);

  // controls
  const r1 = htmlEl("div", "gl-row", "", right);
  const pg = htmlEl("div", "gl-group", "", r1);
  htmlEl("span", "", "start:", pg);
  const PRESETS: [string, number[]][] = [
    ["good", [0.5, 0.1, 10]], ["bad", [3, 0.05, 50]], ["mirror", [0.1, 0.5, 2]], ["far", [10, 4.5, 10]],
  ];
  const pBtns = presets(pg, PRESETS.map(([label, p]) => ({ label, apply: () => setStart(p, true) })));
  const r2 = htmlEl("div", "gl-row", "", right);
  const mg = htmlEl("div", "gl-group", "", r2);
  const mBtns = presets(mg, [
    { label: "GN + LM", apply: () => { st.method = "both"; draw(); } },
    { label: "GN", apply: () => { st.method = "gn"; draw(); } },
    { label: "LM", apply: () => { st.method = "lm"; draw(); } },
  ]);
  const qg = htmlEl("div", "gl-group", "", r2);
  const qBtns = presets(qg, [
    { label: "(ka, ke, V)", apply: () => { st.par = "lin"; restart(); } },
    { label: "log", apply: () => { st.par = "log"; restart(); } },
  ]);
  const r3 = htmlEl("div", "gl-row", "", right);
  const vS = slider(r3, "start V", 0, 2, 0.01, Math.log10(st.start[2]), v => {
    st.start[2] = 10 ** v; pBtns.forEach(b => b.classList.remove("active")); restart();
  }, v => { const V = 10 ** v; return V >= 10 ? V.toFixed(0) : V.toFixed(1); });
  const ag = htmlEl("div", "gl-group", "", r3);
  const stepB = button(ag, "Step", () => { stopRun(); doStep(); }, "gl-btn-main");
  const runB = button(ag, "Run", () => { if (timer !== null) stopRun(); else startRun(); }, "gl-btn-main");
  button(ag, "reset", () => {
    st.start = startDefault.slice(); st.par = init.par; st.method = init.method;
    vS.set(Math.log10(st.start[2])); pBtns.forEach(b => b.classList.remove("active")); restart();
  });
  void stepB;

  const WD = 620, HD = 250;
  const svgD = svgEl("svg", { viewBox: `0 0 ${WD} ${HD}`, class: "gl-data", role: "img",
    "aria-label": "Plasma concentration data with the current model curves" }, right);
  const tab = htmlEl("table", "gl-tab", "", right);

  // ---------------- landscape scales ----------------
  const lx0 = 62, lx1 = WL - 10, ly0 = 34, ly1 = HL - 54;
  const sx = scale(L0, L1, lx0, lx1), sy = scale(L0, L1, ly1, ly0);
  const gx = (i: number) => sx(L0 + (L1 - L0) * i / (NG - 1));
  const gy = (j: number) => sy(L0 + (L1 - L0) * j / (NG - 1));

  // static layer
  svgEl("text", { x: 6, y: 22, class: "ttl" }, svgL).textContent = "(a) RSS over (ka, ke), V at its best";
  const url = heatmap();
  if (url) svgEl("image", { href: url, x: lx0, y: ly0, width: lx1 - lx0, height: ly1 - ly0, preserveAspectRatio: "none" }, svgL);
  const Z = land();
  for (const c of LEVELS) {
    const d = contourSegs(Z, c).map(([a, b, e, f]) => `M${gx(a).toFixed(1)},${gy(b).toFixed(1)}L${gx(e).toFixed(1)},${gy(f).toFixed(1)}`).join("");
    svgEl("path", { d, fill: "none", stroke: "#2b3a4a", "stroke-width": c === 0.07 ? 1.6 : 1, "stroke-opacity": 0.55 }, svgL);
  }
  // ka = ke
  svgEl("line", { x1: sx(L0), y1: sy(L0), x2: sx(L1), y2: sy(L1), stroke: COL.ref, "stroke-width": 2, "stroke-dasharray": "7 5" }, svgL);
  const dl = svgEl("text", { x: sx(1.05), y: sy(1.05) + 26, "text-anchor": "middle", class: "sm",
    transform: `rotate(-45 ${sx(1.05)} ${sy(1.05) + 26})` }, svgL);
  dl.textContent = "ka = ke";
  svgEl("rect", { x: lx0, y: ly0, width: lx1 - lx0, height: ly1 - ly0, fill: "none", stroke: "#444" }, svgL);
  const ticks = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30];
  for (const t of ticks) {
    const l = Math.log10(t);
    svgEl("line", { x1: sx(l), x2: sx(l), y1: ly1, y2: ly1 + 5, class: "axis" }, svgL);
    svgEl("text", { x: sx(l), y: ly1 + 22, "text-anchor": "middle", class: "sm" }, svgL).textContent = String(t);
    svgEl("line", { x1: lx0 - 5, x2: lx0, y1: sy(l), y2: sy(l), class: "axis" }, svgL);
    svgEl("text", { x: lx0 - 8, y: sy(l) + 6, "text-anchor": "end", class: "sm" }, svgL).textContent = String(t);
  }
  svgEl("text", { x: (lx0 + lx1) / 2, y: HL - 8, "text-anchor": "middle" }, svgL).textContent = "absorption rate ka (per hour)";
  const yl = svgEl("text", { x: 16, y: (ly0 + ly1) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(ly0 + ly1) / 2})` }, svgL);
  yl.textContent = "elimination rate ke (per hour)";
  const star = (p: number[], c: string, lab: string, dx: number, dy: number, anchor: string) => {
    const x = sx(Math.log10(p[0])), y = sy(Math.log10(p[1]));
    const pts = Array.from({ length: 10 }, (_, k) => {
      const r = k % 2 ? 5 : 12, a = -Math.PI / 2 + k * Math.PI / 5;
      return `${(x + r * Math.cos(a)).toFixed(1)},${(y + r * Math.sin(a)).toFixed(1)}`;
    }).join(" ");
    svgEl("polygon", { points: pts, fill: c, stroke: "#fff", "stroke-width": 1.5 }, svgL);
    const t = svgEl("text", { x: x + dx, y: y + dy, "text-anchor": anchor, "font-weight": 700, stroke: "#fff", "stroke-width": 4, "paint-order": "stroke" }, svgL);
    t.style.fill = c; t.textContent = lab;
  };
  star(P_HAT, COL.hat, "fit, V = 20.2 L", -12, 30, "end");
  star(P_FLIP, COL.flip, "mirror, V = 3.8 L", 14, -16, "start");
  const dyn = svgEl("g", {}, svgL);
  const handle = svgEl("circle", { r: 12, class: "handle gl-handle", stroke: COL.ink, fill: "#fff", "stroke-width": 3 }, svgL);
  const handleLab = svgEl("text", { "font-weight": 700, stroke: "#fff", "stroke-width": 4, "paint-order": "stroke", "pointer-events": "none" }, svgL);
  handleLab.textContent = "start";

  // click on the map moves the start; the start marker is draggable
  const toWorld = (ev: { clientX: number; clientY: number }): [number, number] => {
    const [px, py] = svgPoint(svgL, ev);
    const cl = (v: number) => Math.max(L0 + 0.01, Math.min(L1 - 0.01, v));
    return [cl(sx.invert(px)), cl(sy.invert(py))];
  };
  const moveStart = (la: number, le: number) => {
    st.start[0] = Number((10 ** la).toPrecision(3)); st.start[1] = Number((10 ** le).toPrecision(3));
    pBtns.forEach(b => b.classList.remove("active"));
    restart();
  };
  draggable(handle, toWorld, (a, b) => moveStart(a, b));
  svgL.addEventListener("pointerdown", ev => {
    if ((ev.target as Element).classList.contains("gl-handle")) return;
    const [px, py] = svgPoint(svgL, ev);
    if (px < lx0 || px > lx1 || py < ly0 || py > ly1) return;
    const [a, b] = toWorld(ev);
    moveStart(a, b);
  });

  // ---------------- data panel scales ----------------
  const dx0 = 58, dx1 = WD - 10, dy0 = 30, dy1 = HD - 44;
  const tx = scale(0, 12.5, dx0, dx1), cy = scale(0, 4, dy1, dy0);
  const tg = Array.from({ length: 251 }, (_, i) => 0.01 + 12.49 * i / 250);

  // ---------------- actions ----------------
  function setStart(p: number[], fromPreset = false): void {
    st.start = p.slice();
    vS.set(Math.log10(p[2]));
    if (!fromPreset) pBtns.forEach(b => b.classList.remove("active"));
    restart();
  }
  function restart(): void { stopRun(); runs = fresh(); draw(); }
  function doStep(): boolean {
    let any = false;
    for (const r of active()) if (r.status === "running") { step(r); any = true; }
    draw();
    return any && active().some(r => r.status === "running");
  }
  function startRun(): void {
    if (!active().some(r => r.status === "running")) return;
    runB.textContent = "Pause";
    timer = window.setInterval(() => { if (!doStep()) stopRun(); }, speed);
  }
  function stopRun(): void {
    if (timer !== null) { window.clearInterval(timer); timer = null; }
    runB.textContent = "Run";
  }

  // ---------------- drawing ----------------
  const statusText: Record<Status, string> = {
    running: "running", converged: "converged", diverged: "diverged", singular: "singular gradient", maxit: "max. iterations",
  };

  function drawPath(r: Run, c: string): void {
    const pts = r.path.map(x => toNat(x, r.par));
    const map = (p: number[]): { x: number; y: number; out: boolean } => {
      // a non-positive value cannot be drawn on a log axis: put it at the low edge
      const lg = (v: number) => v > 0 ? Math.log10(v) : -Infinity;
      const la = lg(p[0]), le = lg(p[1]);
      const cla = Math.max(L0, Math.min(L1, la)), cle = Math.max(L0, Math.min(L1, le));
      return { x: sx(cla), y: sy(cle), out: cla !== la || cle !== le };
    };
    const m = pts.map(map);
    const g = svgEl("g", {}, dyn);
    // rejected LM trials
    for (const x of r.rejected) {
      const q = map(toNat(x, r.par));
      if (q.out) continue;
      svgEl("path", { d: `M${q.x - 6},${q.y - 6}L${q.x + 6},${q.y + 6}M${q.x - 6},${q.y + 6}L${q.x + 6},${q.y - 6}`, stroke: c, "stroke-width": 2.2, "stroke-opacity": 0.8 }, g);
    }
    if (m.length > 1) {
      svgEl("polyline", { points: m.map(q => `${q.x.toFixed(1)},${q.y.toFixed(1)}`).join(" "), fill: "none", stroke: "#fff", "stroke-width": 6, "stroke-opacity": 0.7, "stroke-linejoin": "round" }, g);
      svgEl("polyline", { points: m.map(q => `${q.x.toFixed(1)},${q.y.toFixed(1)}`).join(" "), fill: "none", stroke: c, "stroke-width": 3, "stroke-linejoin": "round" }, g);
    }
    m.forEach((q, k) => {
      if (k === 0) return;
      if (q.out) svgEl("rect", { x: q.x - 7, y: q.y - 7, width: 14, height: 14, fill: "#fff", stroke: c, "stroke-width": 3 }, g);
      else svgEl("circle", { cx: q.x, cy: q.y, r: k === m.length - 1 ? 7 : 4.5, fill: c, stroke: "#fff", "stroke-width": 1.5 }, g);
    });
    const last = m[m.length - 1];
    if (last.out) {
      const p = pts[pts.length - 1];
      const t = svgEl("text", { x: Math.min(Math.max(last.x, lx0 + 175), lx1 - 175), y: Math.min(Math.max(last.y + (last.y < (ly0 + ly1) / 2 ? 30 : -16), ly0 + 20), ly1 - 8), "text-anchor": "middle",
        "font-weight": 700, stroke: "#fff", "stroke-width": 4, "paint-order": "stroke" }, g);
      t.style.fill = c;
      t.textContent = `off the map: ka = ${num(p[0], 3)}, ke = ${num(p[1], 3)}`;
    }
  }

  function draw(): void {
    mBtns.forEach((b, i) => b.classList.toggle("active", ["both", "gn", "lm"][i] === st.method));
    qBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.par === "lin")));
    const showGN = st.method !== "lm", showLM = st.method !== "gn";

    // (a) paths and start
    dyn.replaceChildren();
    if (showLM) drawPath(runs.lm, COL.lm);
    if (showGN) drawPath(runs.gn, COL.gn);
    const hx = sx(Math.log10(st.start[0])), hy = sy(Math.log10(st.start[1]));
    setAttrs(handle, { cx: hx, cy: hy });
    setAttrs(handleLab, { x: hx + 16, y: hy + 6 });

    // (b) data and curves
    svgD.replaceChildren();
    svgEl("text", { x: 6, y: 20, class: "ttl" }, svgD).textContent = "(b) data and current model curves";
    for (let v = 0; v <= 4; v++) {
      svgEl("line", { x1: dx0, x2: dx1, y1: cy(v), y2: cy(v), class: v === 0 ? "axis" : "grid" }, svgD);
      svgEl("text", { x: dx0 - 6, y: cy(v) + 6, "text-anchor": "end", class: "sm" }, svgD).textContent = String(v);
    }
    for (let t = 0; t <= 12; t += 2) {
      svgEl("line", { x1: tx(t), x2: tx(t), y1: dy1, y2: dy1 + 5, class: "axis" }, svgD);
      svgEl("text", { x: tx(t), y: dy1 + 22, "text-anchor": "middle", class: "sm" }, svgD).textContent = String(t);
    }
    svgEl("text", { x: (dx0 + dx1) / 2, y: HD - 4, "text-anchor": "middle" }, svgD).textContent = "time after the dose (h)";
    const yl2 = svgEl("text", { x: 14, y: (dy0 + dy1) / 2, "text-anchor": "middle", transform: `rotate(-90 14 ${(dy0 + dy1) / 2})`, class: "sm" }, svgD);
    yl2.textContent = "mg/L";
    const dclip = `gl-dclip-${uid}`;
    const dd = svgEl("defs", {}, svgD);
    const dc = svgEl("clipPath", { id: dclip }, dd);
    svgEl("rect", { x: dx0, y: dy0 - 6, width: dx1 - dx0, height: dy1 - dy0 + 6 }, dc);
    const cg = svgEl("g", { "clip-path": `url(#${dclip})` }, svgD);
    const offScale: string[] = [];
    const curve = (p: number[], c: string, w: number, dash = "", who = "") => {
      const ys = tg.map(t => bateman(t, p[0], p[1], p[2]));
      if (who && !ys.every(v => Number.isFinite(v) && v > -1 && v < 6)) offScale.push(who);
      if (!ys.every(Number.isFinite)) return;
      const pts = tg.map((t, i) => `${tx(t).toFixed(1)},${cy(Math.max(-1, Math.min(6, ys[i]))).toFixed(1)}`).join(" ");
      svgEl("polyline", { points: pts, fill: "none", stroke: c, "stroke-width": w, "stroke-dasharray": dash }, cg);
    };
    curve(st.start, COL.ref, 2, "6 5");
    if (showLM) curve(toNat(runs.lm.x, st.par), COL.lm, 3.5, "", "LM");
    if (showGN) curve(toNat(runs.gn.x, st.par), COL.gn, 3, showLM ? "10 5" : "", "GN");
    TT.forEach((t, i) => svgEl("circle", { cx: tx(t), cy: cy(Y[i]), r: 5.5, fill: COL.pt, stroke: "#fff", "stroke-width": 1.5 }, svgD));
    let ly = dy0 + 6;
    const lg2 = (c: string, txt: string, dash = "", w = 3) => {
      const x = dx1 - 190;
      svgEl("line", { x1: x, x2: x + 28, y1: ly, y2: ly, stroke: c, "stroke-width": w, "stroke-dasharray": dash }, svgD);
      const t = svgEl("text", { x: x + 36, y: ly + 6, "font-weight": 600 }, svgD);
      t.style.fill = c; t.textContent = txt;
      ly += 24;
    };
    lg2(COL.ref, "start", "6 5", 2);
    if (showGN) lg2(COL.gn, "Gauss–Newton", showLM ? "10 5" : "");
    if (showLM) lg2(COL.lm, "Lev.–Marquardt", "", 3.5);
    if (offScale.length) {
      const t = svgEl("text", { x: dx1 - 4, y: ly + 8, "text-anchor": "end", "font-weight": 700 }, svgD);
      t.style.fill = "#b00020"; t.textContent = `${offScale.join(", ")} curve off the scale`;
    }

    // (c) readout table
    tab.replaceChildren();
    const cols = [showGN ? runs.gn : null, showLM ? runs.lm : null].filter((r): r is Run => r !== null);
    const head = htmlEl("tr", "", "", tab);
    htmlEl("th", "", st.par === "log" ? "in log-parameters" : "in (ka, ke, V)", head);
    for (const r of cols) {
      const th = htmlEl("th", "", r.method === "gn" ? "Gauss–Newton" : "Levenberg–Marquardt", head);
      th.style.color = r.method === "gn" ? COL.gn : COL.lm;
    }
    const row = (label: string, f: (r: Run) => string | [string, string]) => {
      const tr = htmlEl("tr", "", "", tab);
      htmlEl("td", "", label, tr);
      for (const r of cols) {
        const v = f(r);
        if (typeof v === "string") htmlEl("td", "", v, tr);
        else htmlEl("td", v[1], v[0], tr);
      }
    };
    row("pass", r => r.method === "gn" ? String(r.passes) : `${r.passes} (${r.accepted} acc., ${r.nRejected} rej.)`);
    row("ka, ke, V", r => toNat(r.x, r.par).map(v => num(v, 4)).join(", "));
    row("RSS", r => rssFmt(r.rss));
    row("λ", r => r.method === "gn" ? "0 (none)" : lamFmt(r.lambda));
    row("last step", r => {
      if (r.last === "") return "–";
      if (r.last === "taken") return r.passes > 1 && r.rss > rssOf(r, r.path.length - 2) ? ["taken, RSS went up", "gl-rej"] : "taken";
      return r.last === "accepted" ? [`accepted, λ ÷ 10`, "gl-acc"] : [`rejected (RSS ${num(r.lastTrialRss, 3)}), λ × 10`, "gl-rej"];
    });
    row("status", r => {
      const s = statusText[r.status];
      return r.status === "converged" ? [s, "gl-acc"] : r.status === "running" ? s : [s, "gl-rej"];
    });
  }
  const rssCache = new WeakMap<number[], number>();
  function rssOf(r: Run, k: number): number {
    const x = r.path[k];
    let v = rssCache.get(x);
    if (v === undefined) {
      const p = toNat(x, r.par);
      v = TT.reduce((s, t, i) => s + (Y[i] - bateman(t, p[0], p[1], p[2])) ** 2, 0);
      rssCache.set(x, v);
    }
    return v;
  }

  // initial preset highlight
  PRESETS.forEach(([, p], i) => pBtns[i].classList.toggle("active", p.every((v, j) => Math.abs(v - st.start[j]) < 1e-12)));
  draw();
}

autoMount(NAME, mount);

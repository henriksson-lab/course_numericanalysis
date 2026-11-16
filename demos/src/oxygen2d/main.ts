// Oxygen-in-a-chamber explorer for the steady-state-2d unit. The lecture's chamber
// (4000 x 2000 µm, lid 100 µm above a cell layer, D = 2000 µm^2/s, c0 = 200 µM on the open
// edges, 50 amol/s per cell, so R0 = 0.5 µM/s per 1e5 cells/cm^2, Km = 1 µM), solved live by
// nonlinear SOR on the five-point grid (solver.ts), warm-started from the previous field.
//  (a) top-view heat map (viridis, as in the unit's figures), the hypoxic contour (vermillion),
//      the c = 0 contour (white dashed, grey fill below 0), the outline of the cells (orange,
//      dotted) and the cross-section line; edges drawn blue when open, grey when sealed.
//      The mouse moves the cross-section, or paints / erases cells.
//  (b) the cross-section c(x) or c(y), with the threshold and c = 0, the previous grid's
//      curve for comparison; or a grid-check table (h = 200 ... 12.5 µm).
// Auto-mounts on <div data-demo="oxygen2d">. Optional attributes:
//   data-uptake="mm|constant|linear"  data-density="2" (1e5 cells/cm^2)  data-h="50"
//   data-threshold="10"  data-cells="all|spheroid|patches|none"  data-edges="lrbt" (open)
//   data-cut="x|y"  data-panel="section|check"
import { svgEl, htmlEl, slider, presets, checkbox, button, svgPoint, autoMount } from "../lib/ui";
import {
  Grid, Edges, Uptake, Problem, makeGrid, setEdges, solve, stats, nodeCount,
  LX, LY, C0, HFINE, R_PER_1E5,
} from "./solver";

const NAME = "oxygen2d";
// Okabe-Ito, as in make-steady-state-2d-figs.R; viridis (hcl.colors) for the map
const COL = { line: "#0072B2", hyp: "#D55E00", cells: "#CC79A7", open: "#56B4E9", sealed: "#555", ghost: "#999", neg: "#a8a8a8", ink: "#222" };
const VIR = ["#4B0055", "#403173", "#00588B", "#007B97", "#009B95", "#00B785", "#53CC67", "#B2DC3C", "#FDE333"];
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .35em; }
[data-demo="${NAME}"] .ox-row { display: flex; flex-wrap: wrap; gap: .3em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .ox-group { display: flex; gap: .35em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .ox-row .qb-slider { min-width: 250px; flex: 1 1 250px; max-width: 400px; grid-template-columns: auto 1fr 3em; }
[data-demo="${NAME}"] .ox-plots { display: flex; flex-wrap: wrap; gap: 8px; width: 100%; }
[data-demo="${NAME}"] .ox-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] svg .mono { font-family: ui-monospace, monospace; font-size: 18px; }
[data-demo="${NAME}"] .ox-read { width: 100%; box-sizing: border-box; display: flex; flex-wrap: wrap; gap: .1em 1.3em; }
[data-demo="${NAME}"] .ox-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .qb-check { gap: .2em; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 18px; gap: .2em; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 16px; padding: .2em .45em; }
.reveal [data-demo="${NAME}"] .qb-btn { font-size: .85em; padding: .12em .5em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

const HS = [200, 100, 50, 25, 12.5];
const [FX, FY] = nodeCount(HFINE);

function viridis(t: number): [number, number, number] {
  const x = Math.min(1, Math.max(0, t)) * (VIR.length - 1);
  const i = Math.min(VIR.length - 2, Math.floor(x)), u = x - i;
  const a = VIR[i], b = VIR[i + 1];
  const ch = (s: string, o: number) => parseInt(s.slice(o, o + 2), 16);
  return [0, 1, 2].map(m => Math.round(ch(a, 1 + 2 * m) * (1 - u) + ch(b, 1 + 2 * m) * u)) as [number, number, number];
}

/** Marching squares: line segments of the level set Z = lev, as an SVG path. */
function contour(Z: ArrayLike<number>, NX: number, NY: number, lev: number,
  px: (i: number) => number, py: (j: number) => number): string {
  let d = "";
  const f = (v: number) => v.toFixed(1);
  for (let j = 0; j < NY - 1; j++) for (let i = 0; i < NX - 1; i++) {
    const k = i + NX * j;
    const v = [Z[k], Z[k + 1], Z[k + 1 + NX], Z[k + NX]];   // corners: (i,j) (i+1,j) (i+1,j+1) (i,j+1)
    const cx = [i, i + 1, i + 1, i], cy = [j, j, j + 1, j + 1];
    const pts: [number, number][] = [];
    for (let e = 0; e < 4; e++) {
      const a = v[e], b = v[(e + 1) % 4];
      if ((a < lev) !== (b < lev)) {
        const t = (lev - a) / (b - a), e2 = (e + 1) % 4;
        pts.push([px(cx[e] + t * (cx[e2] - cx[e])), py(cy[e] + t * (cy[e2] - cy[e]))]);
      }
    }
    if (pts.length === 2) d += `M${f(pts[0][0])},${f(pts[0][1])}L${f(pts[1][0])},${f(pts[1][1])}`;
    else if (pts.length === 4) {
      const mid = (v[0] + v[1] + v[2] + v[3]) / 4;
      const [p, q] = (mid < lev) === (v[0] < lev) ? [[0, 3], [1, 2]] : [[0, 1], [2, 3]];
      d += `M${f(pts[p[0]][0])},${f(pts[p[0]][1])}L${f(pts[p[1]][0])},${f(pts[p[1]][1])}`;
      d += `M${f(pts[q[0]][0])},${f(pts[q[0]][1])}L${f(pts[q[1]][0])},${f(pts[q[1]][1])}`;
    }
  }
  return d;
}

type CellPreset = "all" | "spheroid" | "patches" | "none";
function fillMask(mask: Uint8Array, kind: CellPreset): void {
  for (let j = 0; j < FY; j++) for (let i = 0; i < FX; i++) {
    const x = i * HFINE, y = j * HFINE;
    let v = 0;
    if (kind === "all") v = 1;
    else if (kind === "spheroid") v = Math.hypot(x - 2000, y - 1000) <= 700 ? 1 : 0;
    else if (kind === "patches") v = (Math.hypot(x - 1100, y - 1000) <= 550 || Math.hypot(x - 2900, y - 1000) <= 550) ? 1 : 0;
    mask[i + FX * j] = v;
  }
}

const fmtH = (h: number) => (h === 12.5 ? "12.5" : String(h));
const fmtInt = (n: number) => String(n).replace(/\B(?=(\d{3})+(?!\d))/g, " ");
function fmtC(v: number): string {
  const a = Math.abs(v);
  return a >= 10 ? v.toFixed(2) : a >= 0.1 ? v.toFixed(3) : a >= 0.001 ? v.toFixed(4) : v.toExponential(1);
}

let uidCount = 0;

interface Row { h: number; g: Grid; sweeps: number; done: boolean }

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (s: string | undefined, d: number, lo: number, hi: number) => {
    const v = Number(s);
    return s !== undefined && Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : d;
  };
  const start = {
    uptake: (["constant", "linear", "mm"].includes(ds.uptake ?? "") ? ds.uptake : "mm") as Uptake,
    rho: num(ds.density, 2, 0, 4),
    h: HS.includes(Number(ds.h)) ? Number(ds.h) : 50,
    thr: num(ds.threshold, 10, 1, 50),
    cells: (["all", "spheroid", "patches", "none"].includes(ds.cells ?? "") ? ds.cells : "all") as CellPreset,
    edges: ((s: string) => ({ l: s.includes("l"), r: s.includes("r"), b: s.includes("b"), t: s.includes("t") }))(
      /[lrbt]/.test(ds.edges ?? "") ? ds.edges ?? "lrbt" : "lrbt"),
    cut: ds.cut === "y" ? "y" : "x",
    panel: ds.panel === "check" ? "check" : "section",
  };
  const st = {
    uptake: start.uptake, rho: start.rho, h: start.h, thr: start.thr,
    edges: { ...start.edges } as Edges,
    cut: start.cut as "x" | "y", cutPos: start.cut === "x" ? 1000 : 2000,
    mode: "cut" as "cut" | "paint" | "erase",
    panel: start.panel as "section" | "check",
  };
  const mask = new Uint8Array(FX * FY);
  fillMask(mask, start.cells);
  const prob = (): Problem => ({ uptake: st.uptake, R0: R_PER_1E5 * st.rho, edges: st.edges });

  let g = makeGrid(st.h, st.edges, mask);
  let sweeps = 0, ms = 0, done = false, warm = false;
  let ghost: { h: number; xs: number[]; ys: number[] } | null = null;
  let rows: Row[] = [];
  let maskPath = "";

  // ------------------------------------------------------------------ controls
  const r1 = htmlEl("div", "ox-row", "", host);
  const r2 = htmlEl("div", "ox-row", "", host);
  const r3 = htmlEl("div", "ox-row", "", host);
  // row 1: uptake law, cell density, reset
  const ug = htmlEl("div", "ox-group", "", r1);
  htmlEl("span", "", "uptake:", ug);
  const UP: Uptake[] = ["constant", "linear", "mm"];
  const uBtns = presets(ug, [
    { label: "constant", apply: () => setUptake("constant") },
    { label: "linear kc", apply: () => setUptake("linear") },
    { label: "Michaelis–Menten", apply: () => setUptake("mm") },
  ]);
  const rhoS = slider(r1, "cells (10⁵/cm²)", 0, 4, 0.05, st.rho, v => { st.rho = v; changed(false); }, v => v.toFixed(2));
  // row 2: grid, threshold, open edges
  const hg = htmlEl("div", "ox-group", "", r2);
  htmlEl("span", "", "grid h (µm):", hg);
  const hBtns = presets(hg, HS.map(h => ({ label: fmtH(h), apply: () => setH(h) })));
  const thrS = slider(r2, "hypoxic < (µM)", 1, 50, 1, st.thr, v => { st.thr = v; draw(); }, v => v.toFixed(0));
  const eg = htmlEl("div", "ox-group", "", r2);
  htmlEl("span", "", "open:", eg);
  const eKeys: (keyof Edges)[] = ["l", "r", "b", "t"];
  const eLab = { l: "left", r: "right", b: "bottom", t: "top" };
  const eBoxes = eKeys.map(key => checkbox(eg, eLab[key], st.edges[key], v => {
    st.edges[key] = v;
    if (!eKeys.some(q => st.edges[q])) { st.edges[key] = true; eBoxes[eKeys.indexOf(key)].checked = true; note = "At least one edge must stay open: without supply there is no steady state with oxygen."; draw(); return; }
    setEdges(g, st.edges); changed(false);
  }));
  // row 3: cells, mouse mode, cut direction, right panel
  const cg = htmlEl("div", "ox-group", "", r3);
  htmlEl("span", "", "cells:", cg);
  presets(cg, [
    { label: "everywhere", apply: () => setCells("all") },
    { label: "spheroid", apply: () => setCells("spheroid") },
    { label: "two patches", apply: () => setCells("patches") },
    { label: "none", apply: () => setCells("none") },
  ]).forEach(b => b.addEventListener("click", () => b.classList.remove("active")));
  const mg = htmlEl("div", "ox-group", "", r3);
  htmlEl("span", "", "drag:", mg);
  const mBtns = presets(mg, [
    { label: "move cut", apply: () => { st.mode = "cut"; draw(); } },
    { label: "paint", apply: () => { st.mode = "paint"; draw(); } },
    { label: "erase", apply: () => { st.mode = "erase"; draw(); } },
  ]);
  const xg = htmlEl("div", "ox-group", "", r3);
  const cutBtns = presets(xg, [
    { label: "cut x", apply: () => { if (st.cut !== "x") { st.cut = "x"; st.cutPos = 1000; } ghost = null; draw(); } },
    { label: "cut y", apply: () => { if (st.cut !== "y") { st.cut = "y"; st.cutPos = 2000; } ghost = null; draw(); } },
  ]);
  const pg = htmlEl("div", "ox-group", "", r3);
  const pBtns = presets(pg, [
    { label: "section", apply: () => { st.panel = "section"; draw(); } },
    { label: "grid check", apply: () => { st.panel = "check"; kick(); draw(); } },
  ]);
  button(r1, "reset", () => {
    st.uptake = start.uptake; st.rho = start.rho; st.thr = start.thr; st.edges = { ...start.edges };
    st.cut = start.cut as "x" | "y"; st.cutPos = st.cut === "x" ? 1000 : 2000; st.mode = "cut"; st.panel = start.panel as "section" | "check";
    rhoS.set(st.rho); thrS.set(st.thr); eKeys.forEach((q, i) => (eBoxes[i].checked = st.edges[q]));
    fillMask(mask, start.cells); updateMaskPath();
    ghost = null; st.h = start.h; g = makeGrid(st.h, st.edges, mask); warm = false;
    changed(true);
  });

  // ------------------------------------------------------------------ plots
  const plots = htmlEl("div", "ox-plots", "", host);
  const WA = 660, WB = 470, H = 362, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Oxygen in the chamber, top view" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 4px)`;
  svgA.style.touchAction = "none";
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Cross-section of the oxygen field" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 4px)`;
  const ro = htmlEl("div", "qb-readout ox-read", "", host);
  let note = "";

  // map geometry: 0.125 px per µm
  const ML = 80, MT = 44, MW = 500, MH = 250;
  const mx = (x: number) => ML + x * MW / LX, my = (y: number) => MT + MH - y * MH / LY;
  const mxInv = (px: number) => (px - ML) * LX / MW, myInv = (py: number) => (MT + MH - py) * LY / MH;

  // static parts of (a)
  const defs = svgEl("defs", {}, svgA);
  const clip = svgEl("clipPath", { id: `ox-clip-${++uidCount}` }, defs);
  svgEl("rect", { x: ML, y: MT, width: MW, height: MH }, clip);
  svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent = "(a) O₂ under the lid, top view";
  const img = svgEl("image", { x: 0, y: 0, width: 1, height: 1, preserveAspectRatio: "none", "clip-path": `url(#${clip.id})` }, svgA);
  const negPath = svgEl("path", { fill: "none", stroke: "#fff", "stroke-width": 2.5, "stroke-dasharray": "7 5" }, svgA);
  const hypPath = svgEl("path", { fill: "none", stroke: COL.hyp, "stroke-width": 3.5 }, svgA);
  const cellPath = svgEl("path", { fill: "none", stroke: COL.cells, "stroke-width": 3, "stroke-dasharray": "4 4" }, svgA);
  const cutLine = svgEl("line", { stroke: COL.ink, "stroke-width": 2, "stroke-dasharray": "10 5" }, svgA);
  const edgeLines = eKeys.map(() => svgEl("line", { "stroke-width": 7, "stroke-linecap": "square" }, svgA));
  const minMark = svgEl("circle", { r: 6, fill: "none", stroke: "#fff", "stroke-width": 2.5 }, svgA);
  for (let x = 0; x <= LX; x += 1000) {
    svgEl("line", { x1: mx(x), x2: mx(x), y1: MT + MH + 4, y2: MT + MH + 9, class: "axis" }, svgA);
    svgEl("text", { x: mx(x), y: MT + MH + 28, "text-anchor": "middle", class: "sm" }, svgA).textContent = String(x);
  }
  for (let y = 0; y <= LY; y += 1000) {
    svgEl("text", { x: ML - 10, y: my(y) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = String(y);
  }
  svgEl("text", { x: ML + MW / 2, y: MT + MH + 52, "text-anchor": "middle" }, svgA).textContent = "x (µm)";
  const ylab = svgEl("text", { x: 16, y: MT + MH / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${MT + MH / 2})` }, svgA);
  ylab.textContent = "y (µm)";
  // colour bar
  const CBX = ML + MW + 18;
  for (let s = 0; s < 50; s++) {
    const [r, gg, b] = viridis((s + 0.5) / 50);
    svgEl("rect", { x: CBX, y: MT + MH - (s + 1) * MH / 50, width: 16, height: MH / 50 + 0.6, fill: `rgb(${r},${gg},${b})` }, svgA);
  }
  for (const v of [0, 100, 200]) {
    svgEl("text", { x: CBX + 20, y: MT + MH - v * MH / 200 + 6, class: "sm" }, svgA).textContent = String(v);
  }
  svgEl("text", { x: CBX - 4, y: MT - 10, class: "sm" }, svgA).textContent = "µM";
  svgEl("rect", { x: CBX, y: MT + MH + 12, width: 16, height: 16, fill: COL.neg }, svgA);
  svgEl("text", { x: CBX + 20, y: MT + MH + 26, class: "sm" }, svgA).textContent = "< 0";

  // ------------------------------------------------------------------ mouse
  let painting = false;
  const act = (ev: PointerEvent) => {
    const [px, py] = svgPoint(svgA, ev);
    const x = mxInv(px), y = myInv(py);
    if (st.mode === "cut") {
      if (st.cut === "x") st.cutPos = Math.min(LY, Math.max(0, y));
      else st.cutPos = Math.min(LX, Math.max(0, x));
      draw();
      return;
    }
    const v = st.mode === "paint" ? 1 : 0, rad = 250;
    let any = false;
    const i0 = Math.max(0, Math.floor((x - rad) / HFINE)), i1 = Math.min(FX - 1, Math.ceil((x + rad) / HFINE));
    const j0 = Math.max(0, Math.floor((y - rad) / HFINE)), j1 = Math.min(FY - 1, Math.ceil((y + rad) / HFINE));
    for (let j = j0; j <= j1; j++) for (let i = i0; i <= i1; i++) {
      if (Math.hypot(i * HFINE - x, j * HFINE - y) <= rad && mask[i + FX * j] !== v) { mask[i + FX * j] = v; any = true; }
    }
    if (any) { updateMaskPath(); refreshCells(); changed(false); }
  };
  svgA.addEventListener("pointerdown", ev => {
    const [px, py] = svgPoint(svgA, ev);
    if (px < ML - 10 || px > ML + MW + 10 || py < MT - 10 || py > MT + MH + 10) return;
    ev.preventDefault();
    painting = true;
    svgA.setPointerCapture(ev.pointerId);
    act(ev);
  });
  svgA.addEventListener("pointermove", ev => { if (painting) act(ev); });
  const stop = () => { painting = false; };
  svgA.addEventListener("pointerup", stop);
  svgA.addEventListener("pointercancel", stop);

  // ------------------------------------------------------------------ state changes
  function setUptake(u: Uptake): void { st.uptake = u; changed(false); }
  function setH(h: number): void {
    if (h === st.h) return;
    ghost = section(g);
    st.h = h;
    g = makeGrid(h, st.edges, mask, g);   // warm start: interpolate the old field
    sweeps = 0; ms = 0; done = false; warm = true;
    kick(); draw();
  }
  function setCells(kind: CellPreset): void { fillMask(mask, kind); updateMaskPath(); refreshCells(); changed(false); }
  function refreshCells(): void {
    const m = Math.round(g.h / HFINE);
    for (let j = 0; j < g.NY; j++) for (let i = 0; i < g.NX; i++) g.cells[i + g.NX * j] = mask[i * m + FX * j * m];
  }
  function changed(cold: boolean): void {
    note = "";
    ghost = null;
    rows = [];
    sweeps = 0; ms = 0; done = false; warm = !cold;
    kick(); draw();
  }
  function updateMaskPath(): void {
    // outline of the cells on a 25 µm grid
    const s = 2, NX = (FX - 1) / s + 1, NY = (FY - 1) / s + 1;
    const Z = new Float64Array(NX * NY);
    for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) Z[i + NX * j] = mask[i * s + FX * j * s];
    let all = true;
    for (let k = 0; k < Z.length; k++) if (!Z[k]) { all = false; break; }
    maskPath = all ? "" : contour(Z, NX, NY, 0.5, i => mx(i * s * HFINE), j => my(j * s * HFINE));
  }
  updateMaskPath();

  // ------------------------------------------------------------------ solve loop
  let raf = 0;
  function kick(): void { if (!raf) raf = requestAnimationFrame(tick); }
  function tick(): void {
    raf = 0;
    const p = prob();
    if (!done) {
      const t0 = performance.now();
      const r = solve(g, p, 40);
      ms += performance.now() - t0; sweeps += r.sweeps; done = r.done;
      draw();
      if (!done) { kick(); return; }
    }
    if (st.panel === "check") {
      // grid check: h = 200 ... 12.5, each started from the previous (coarser) field
      let row = rows.find(rw => !rw.done);
      if (!row && rows.length < HS.length) {
        const h = HS[rows.length];
        const prev = rows.length ? rows[rows.length - 1].g : undefined;
        row = { h, g: makeGrid(h, st.edges, mask, prev), sweeps: 0, done: false };
        rows.push(row);
      }
      if (row) {
        const r = solve(row.g, p, 40);
        row.sweeps += r.sweeps; row.done = r.done;
        draw();
        kick();
      }
    }
  }

  // ------------------------------------------------------------------ drawing
  function section(gr: Grid): { h: number; xs: number[]; ys: number[] } {
    const xs: number[] = [], ys: number[] = [];
    if (st.cut === "x") {
      const j = Math.round(st.cutPos / gr.h);
      for (let i = 0; i < gr.NX; i++) { xs.push(i * gr.h); ys.push(gr.c[i + gr.NX * j]); }
    } else {
      const i = Math.round(st.cutPos / gr.h);
      for (let j = 0; j < gr.NY; j++) { xs.push(j * gr.h); ys.push(gr.c[i + gr.NX * j]); }
    }
    return { h: gr.h, xs, ys };
  }

  function drawMap(): void {
    const { NX, NY, c, h } = g;
    const cv = document.createElement("canvas");
    cv.width = NX; cv.height = NY;
    const ctx = cv.getContext("2d");
    if (ctx) {
      const im = ctx.createImageData(NX, NY);
      for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
        const v = c[i + NX * j], o = ((NY - 1 - j) * NX + i) * 4;
        let rgb: [number, number, number] = [168, 168, 168];
        if (v >= 0) rgb = viridis(v / C0);
        im.data[o] = rgb[0]; im.data[o + 1] = rgb[1]; im.data[o + 2] = rgb[2]; im.data[o + 3] = 255;
      }
      ctx.putImageData(im, 0, 0);
      // one pixel per node, centred on it
      img.setAttribute("href", cv.toDataURL());
      img.setAttribute("x", String(mx(-h / 2))); img.setAttribute("y", String(my(LY + h / 2)));
      img.setAttribute("width", String(mx(LX + h / 2) - mx(-h / 2))); img.setAttribute("height", String(my(-h / 2) - my(LY + h / 2)));
      img.style.imageRendering = "pixelated";
    }
    const px = (i: number) => mx(i * h), py = (j: number) => my(j * h);
    hypPath.setAttribute("d", contour(c, NX, NY, st.thr, px, py));
    negPath.setAttribute("d", contour(c, NX, NY, 0, px, py));
    cellPath.setAttribute("d", maskPath);
    const ends: [number, number, number, number][] = [
      [0, 0, 0, LY], [LX, 0, LX, LY], [0, 0, LX, 0], [0, LY, LX, LY],
    ];
    eKeys.forEach((q, i) => {
      const [x1, y1, x2, y2] = ends[i];
      edgeLines[i].setAttribute("x1", String(mx(x1))); edgeLines[i].setAttribute("y1", String(my(y1)));
      edgeLines[i].setAttribute("x2", String(mx(x2))); edgeLines[i].setAttribute("y2", String(my(y2)));
      edgeLines[i].setAttribute("stroke", st.edges[q] ? COL.open : COL.sealed);
    });
    if (st.cut === "x") {
      const y = Math.round(st.cutPos / h) * h;
      Object.entries({ x1: mx(0), x2: mx(LX), y1: my(y), y2: my(y) }).forEach(([a, v]) => cutLine.setAttribute(a, String(v)));
    } else {
      const x = Math.round(st.cutPos / h) * h;
      Object.entries({ x1: mx(x), x2: mx(x), y1: my(0), y2: my(LY) }).forEach(([a, v]) => cutLine.setAttribute(a, String(v)));
    }
  }

  function drawSection(s0: { min: number }): void {
    svgB.replaceChildren();
    const L = 74, R = WB - 26, T = 44, B = 294;
    if (st.panel === "check") { drawCheck(); return; }
    const sec = section(g);
    const len = st.cut === "x" ? LX : LY;
    const lo = Math.min(0, Math.floor(Math.min(s0.min, ...(ghost ? ghost.ys : [0])) / 50) * 50);
    const sx = (x: number) => L + (R - L) * x / len, sy = (v: number) => B - (B - T) * (v - lo) / (210 - lo);
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgB).textContent =
      st.cut === "x" ? `(b) along x at y = ${Math.round(st.cutPos / g.h) * g.h} µm` : `(b) along y at x = ${Math.round(st.cutPos / g.h) * g.h} µm`;
    if (lo < 0) {
      svgEl("rect", { x: L, y: sy(0), width: R - L, height: B - sy(0), fill: "#ececec" }, svgB);
      svgEl("text", { x: (L + R) / 2, y: (sy(0) + B) / 2 + 6, "text-anchor": "middle", class: "sm" }, svgB).textContent = "impossible: c < 0";
    }
    const step = 210 - lo > 300 ? 100 : 50;
    for (let v = lo; v <= 200 + 1e-9; v += step) {
      svgEl("line", { x1: L, x2: R, y1: sy(v), y2: sy(v), class: v === 0 ? "axis" : "grid" }, svgB);
      svgEl("text", { x: L - 8, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = String(v);
    }
    const tstep = len > 2000 ? 1000 : 500;
    for (let x = 0; x <= len; x += tstep) {
      svgEl("line", { x1: sx(x), x2: sx(x), y1: B, y2: B + 5, class: "axis" }, svgB);
      svgEl("text", { x: sx(x), y: B + 24, "text-anchor": "middle", class: "sm" }, svgB).textContent = String(x);
    }
    svgEl("text", { x: (L + R) / 2, y: B + 44, "text-anchor": "middle" }, svgB).textContent = st.cut === "x" ? "x (µm)" : "y (µm)";
    const yl = svgEl("text", { x: 16, y: (T + B) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(T + B) / 2})` }, svgB);
    yl.textContent = "O₂ (µM)";
    // where the cells are along the cut
    for (let k = 0; k < sec.xs.length; k++) {
      const i = st.cut === "x" ? k : Math.round(st.cutPos / g.h), j = st.cut === "x" ? Math.round(st.cutPos / g.h) : k;
      if (g.cells[i + g.NX * j] > 0) {
        const x0 = Math.max(0, sec.xs[k] - g.h / 2), x1 = Math.min(len, sec.xs[k] + g.h / 2);
        svgEl("rect", { x: sx(x0), y: B - 7, width: sx(x1) - sx(x0) + 0.5, height: 7, fill: COL.cells }, svgB);
      }
    }
    // threshold
    svgEl("line", { x1: L, x2: R, y1: sy(st.thr), y2: sy(st.thr), stroke: COL.hyp, "stroke-width": 2.5, "stroke-dasharray": "8 5" }, svgB);
    svgEl("text", { x: R - 4, y: sy(st.thr) - 6, "text-anchor": "end", fill: COL.hyp }, svgB).textContent = `hypoxic < ${st.thr} µM`;
    const poly = (xs: number[], ys: number[]) => xs.map((x, k) => `${sx(x).toFixed(1)},${sy(ys[k]).toFixed(1)}`).join(" ");
    if (ghost && ghost.xs.length) {
      svgEl("polyline", { points: poly(ghost.xs, ghost.ys), fill: "none", stroke: COL.ghost, "stroke-width": 3, "stroke-dasharray": "6 4" }, svgB);
      svgEl("text", { x: L + 8, y: T + 46, fill: "#666" }, svgB).textContent = `- - previous grid, h = ${fmtH(ghost.h)} µm`;
    }
    svgEl("polyline", { points: poly(sec.xs, sec.ys), fill: "none", stroke: COL.line, "stroke-width": 3 }, svgB);
    if (sec.xs.length <= 41) for (let k = 0; k < sec.xs.length; k++) {
      svgEl("circle", { cx: sx(sec.xs[k]), cy: sy(sec.ys[k]), r: 4, fill: COL.line }, svgB);
    }
    svgEl("text", { x: L + 8, y: T + 20, fill: COL.line }, svgB).textContent = `— h = ${fmtH(g.h)} µm`;
  }

  function drawCheck(): void {
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgB).textContent = "(b) grid check: halve h";
    const xr = [12, 190, 270, 380, 458];
    ["h (µm)", "unknowns", "sweeps", "min c", "hypoxic"].forEach((t, i) => {
      svgEl("text", { x: xr[i], y: 66, "text-anchor": i ? "end" : "start", "font-weight": 600 }, svgB).textContent = t;
    });
    svgEl("line", { x1: 8, x2: WB - 8, y1: 76, y2: 76, class: "axis" }, svgB);
    HS.forEach((h, k) => {
      const y = 104 + k * 30;
      const row = rows[k];
      const cur = h === g.h;
      if (cur) svgEl("rect", { x: 6, y: y - 21, width: WB - 12, height: 28, fill: "#fde7d8" }, svgB);
      const cells: string[] = [fmtH(h)];
      if (row) {
        const s = stats(row.g, st.thr);
        cells.push(fmtInt(s.nUnk), row.done ? String(row.sweeps) : "…", fmtC(s.min), (100 * s.hyp).toFixed(1) + " %");
      } else cells.push("", "", "", "");
      cells.forEach((t, i) => {
        svgEl("text", { x: xr[i], y, "text-anchor": i ? "end" : "start", class: "mono" }, svgB).textContent = t;
      });
    });
    const allDone = rows.length === HS.length && rows.every(r => r.done);
    const t = svgEl("text", { x: 12, y: 104 + 5 * 30 + 10, class: "sm" }, svgB);
    t.textContent = allDone ? "each grid started from the coarser one" : "solving …";
    svgEl("text", { x: 12, y: 104 + 5 * 30 + 36, class: "sm" }, svgB).textContent = "min c: second order; area: first order";
  }

  function draw(): void {
    uBtns.forEach((b, i) => b.classList.toggle("active", UP[i] === st.uptake));
    hBtns.forEach((b, i) => b.classList.toggle("active", HS[i] === st.h));
    pBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.panel === "section")));
    mBtns.forEach((b, i) => b.classList.toggle("active", ["cut", "paint", "erase"][i] === st.mode));
    cutBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.cut === "x")));
    const s = stats(g, st.thr);
    drawMap();
    const kx = s.kmin % g.NX, ky = Math.floor(s.kmin / g.NX);
    minMark.setAttribute("cx", String(mx(kx * g.h))); minMark.setAttribute("cy", String(my(ky * g.h)));
    drawSection(s);

    const nCellsAny = g.cells.some((v, k) => v > 0 && !g.fixed[k]);
    const parts: string[] = [];
    parts.push(`<span>min c = <b>${fmtC(s.min)} µM</b> at (${kx * g.h}, ${ky * g.h}) µm</span>`);
    parts.push(`<span>hypoxic (&lt; ${st.thr} µM): <b>${(100 * s.hyp).toFixed(1)} %</b> of the chamber${nCellsAny ? `, ${(100 * s.hypCells).toFixed(1)} % of the cells` : ""}</span>`);
    parts.push(`<span>grid ${g.NX - 2} × ${g.NY - 2} interior, ${fmtInt(s.nUnk)} unknowns · SOR ω = ${g.omega.toFixed(3)} · ${done ? "" : "solving … "}${sweeps} sweeps, ${ms.toFixed(0)} ms (${warm ? "warm start" : "from c = 200"})</span>`);
    if (st.uptake === "constant" && s.min < 0) {
      parts.push(`<span class="ox-warn">c &lt; 0 on ${(100 * s.neg).toFixed(1)} % of the chamber: impossible. Constant uptake keeps consuming oxygen that is not there.</span>`);
    } else if (st.uptake === "linear") {
      parts.push(`<span>linear uptake k = R₀/c₀ = ${(R_PER_1E5 * st.rho / C0).toExponential(1)} /s: right at the edge, ${(C0 / 10 / (1 + 1 / 10)).toFixed(0)}× too low at 10 µM</span>`);
    }
    if (note) parts.push(`<span class="ox-warn">${note}</span>`);
    ro.innerHTML = parts.join("");
  }

  changed(true);
}

autoMount(NAME, mount);

// Spectral unmixing explorer for the pivoting-conditioning unit. Two dyes with Gaussian
// emission spectra (width 30 nm), detectors at the two peaks, overlap c, so the mixing
// matrix is A = [1 c; c 1] and kappa_2 = (1 + c) / (1 - c) (c = 0.95 gives the lecture's 39).
//  (a) the two spectra; the peak spacing follows from c = exp(-d^2 / (2 * 30^2));
//  (b) detector space: readings b of simulated cells that all carry 1 unit of each dye,
//      b = A (1, 1) + noise; the dashed circle is 2 sd of the noise. One cell ("your cell",
//      vermillion) can be dragged: its reading is b0 + delta b;
//  (c) dye space: the two equations of your cell as lines (solid: b0, dashed: b0 + delta b),
//      their intersection (the unmixed amounts), and the unmixed cloud A^-1 b with its
//      predicted 2-sd ellipse (half-axes 2 sd / (1 + c) and 2 sd / (1 - c)).
// Fixed seed, so every run shows the same cells.
// Auto-mounts on <div data-demo="unmixing">. Optional attributes:
//   data-c="0.95"  data-noise="0.01"  data-db="0.01 -0.01"  data-cells="400"  data-seed="2"
import { svgEl, setAttrs, htmlEl, slider, presets, button, checkbox, scale, svgPoint, draggable, rng, numList, autoMount } from "../lib/ui";

const NAME = "unmixing";
// Okabe-Ito, as in make-pivoting-conditioning-figs.R
const COL = { d1: "#0072B2", d2: "#E69F00", bad: "#D55E00", grn: "#009E73", grey: "#8c8c8c", cell: "#555", ink: "#222" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .45em; }
[data-demo="${NAME}"] .um-row { display: flex; flex-wrap: wrap; gap: .4em 1.3em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .um-group { display: flex; gap: .4em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .um-row .qb-slider { min-width: 300px; flex: 1 1 300px; max-width: 470px; grid-template-columns: auto 1fr 4.2em; }
[data-demo="${NAME}"] .um-plots { display: flex; flex-wrap: wrap; gap: 8px; width: 100%; }
[data-demo="${NAME}"] .um-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .um-read { display: flex; flex-wrap: wrap; gap: .15em 1.4em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .um-read span b { font-weight: 600; }
[data-demo="${NAME}"] .um-kappa { font-weight: 700; color: ${COL.bad}; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 20px; gap: .3em; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 18px; padding: .3em .5em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

const WIDTH = 30;                                  // spectrum width (nm), as in the lecture
const spacing = (c: number): number => WIDTH * Math.sqrt(-2 * Math.log(c));
const kappa = (c: number): number => (1 + c) / (1 - c);

/** Solve [1 c; c 1] x = b by elimination (no swap needed: |c| < 1 is already the smaller entry). */
function solve2(c: number, b1: number, b2: number): [number, number] {
  const m = c;                      // multiplier a21 / a11
  const u22 = 1 - m * c;
  const y2 = b2 - m * b1;
  const x2 = y2 / u22;
  return [b1 - c * x2, x2];
}

const f2 = (v: number) => v.toFixed(2);
const f3 = (v: number) => v.toFixed(3);
const sgn = (v: number, d = 3) => (v < -5e-13 ? "−" : "+") + Math.abs(v).toFixed(d);
const pct = (v: number) => (v >= 0.1 ? (100 * v).toFixed(1) : (100 * v).toFixed(2)) + " %";

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (s: string | undefined, d: number) => (s !== undefined && s !== "" && Number.isFinite(Number(s)) ? Number(s) : d);
  const C0 = Math.min(0.99, Math.max(0.1, num(ds.c, 0.95)));
  const S0 = Math.min(0.02, Math.max(0, num(ds.noise, 0.01)));
  const DB0 = numList(ds.db, [0, 0]);
  const NCELL = Math.max(10, Math.min(2000, Math.round(num(ds.cells, 400))));
  const r = rng(num(ds.seed, 2));
  const Z = Array.from({ length: NCELL }, () => [r.norm(), r.norm()] as [number, number]);
  const uid = ++mountCount;

  const st = { c: C0, s: S0, db: [DB0[0] ?? 0, DB0[1] ?? 0] as [number, number], cloud: true, lines: true };
  const BW = 0.05;      // detector panel: half-width around b0
  const XW = 0.6;       // dye panel: half-width around (1, 1)
  const clampDb = (v: number) => Math.max(-BW * 0.96, Math.min(BW * 0.96, v));

  // ---- controls ----
  const row1 = htmlEl("div", "um-row", "", host);
  const cS = slider(row1, "overlap c", 0.1, 0.99, 0.01, st.c, v => { st.c = v; markPresets(); draw(); }, v => v.toFixed(2));
  const cg = htmlEl("div", "um-group", "", row1);
  const cVals = [0.5, 0.9, 0.95, 0.99];
  const cBtns = presets(cg, cVals.map(v => ({ label: "c = " + v, apply: () => { st.c = v; cS.set(v); draw(); } })));
  const nS = slider(row1, "noise sd", 0, 0.02, 0.001, st.s, v => { st.s = v; draw(); }, v => v.toFixed(3));
  const row2 = htmlEl("div", "um-row", "", host);
  const dg = htmlEl("div", "um-group", "", row2);
  htmlEl("span", "", "your cell, Δb =", dg);
  const dbVals: [number, number][] = [[0, 0], [0.01, -0.01], [0.01, 0.01], [0.01, 0]];
  const dbBtns = presets(dg, dbVals.map(([a, b]) => ({
    label: a === 0 && b === 0 ? "0" : `(${sgn(a, 2)}, ${sgn(b, 2)})`,
    apply: () => { st.db = [a, b]; draw(); },
  })));
  htmlEl("span", "qb-note", "or drag the red point in (b)", dg);
  const cb = checkbox(row2, "cells", st.cloud, v => { st.cloud = v; draw(); });
  const lb = checkbox(row2, "lines", st.lines, v => { st.lines = v; draw(); });
  button(row2, "reset", () => {
    st.c = C0; st.s = S0; st.db = [DB0[0] ?? 0, DB0[1] ?? 0]; st.cloud = true; st.lines = true;
    cS.set(st.c); nS.set(st.s); cb.checked = true; lb.checked = true; markPresets(); draw();
  });
  function markPresets(): void {
    cBtns.forEach((b, i) => b.classList.toggle("active", Math.abs(cVals[i] - st.c) < 1e-9));
  }

  // ---- panels ----
  const plots = htmlEl("div", "um-plots", "", host);
  const H = 370, WA = 330, WB = 350, WC = 440, WT = WA + WB + WC;
  const mk = (w: number, label: string) => {
    const s = svgEl("svg", { viewBox: `0 0 ${w} ${H}`, role: "img", "aria-label": label }, plots);
    s.style.width = `calc(${(100 * w / WT).toFixed(2)}% - 6px)`;
    return s;
  };
  const svgA = mk(WA, "Emission spectra of the two dyes");
  const svgB = mk(WB, "Detector readings of the cells");
  const layB = svgEl("g", {}, svgB);     // redrawn; the drag handle stays on top, never detached
  const svgC = mk(WC, "Unmixed dye amounts and the two equations as lines");
  const ro = htmlEl("div", "qb-readout um-read", "", host);

  const top = 36, bot = H - 50;
  // (a)
  const ax = scale(430, 640, 58, WA - 10), ay = scale(0, 1.62, bot, top);
  // (b) square plot area, centred on b0
  const bL = 62, bR = WB - 14, bSide = Math.min(bR - bL, bot - top);
  const bxs = scale(-BW, BW, bL, bL + bSide), bys = scale(-BW, BW, top + bSide, top);
  // (c) square plot area, centred on (1, 1)
  const cL = 62, cSide = Math.min(WC - 14 - cL, bot - top);
  const cxs = scale(1 - XW, 1 + XW, cL, cL + cSide), cys = scale(1 - XW, 1 + XW, top + cSide, top);

  // drag handle for your cell (kept across redraws)
  const handle = svgEl("circle", { r: 11, fill: COL.bad, stroke: "#fff", "stroke-width": 2.5 }, svgB);
  draggable(handle, ev => {
    const [px, py] = svgPoint(svgB, ev);
    return [clampDb(bxs.invert(px)), clampDb(bys.invert(py))];
  }, (dx, dy) => {
    st.db = [Math.round(dx * 1000) / 1000, Math.round(dy * 1000) / 1000];
    draw();
  });

  const axisLabels = (svg: Element, x0: number, x1: number, y0: number, y1: number, xl: string, yl: string) => {
    svgEl("text", { x: (x0 + x1) / 2, y: H - 8, "text-anchor": "middle" }, svg).textContent = xl;
    const ym = (y0 + y1) / 2;
    svgEl("text", { x: 16, y: ym, "text-anchor": "middle", transform: `rotate(-90 16 ${ym})` }, svg).textContent = yl;
  };

  function draw(): void {
    const c = st.c, s = st.s, k = kappa(c), d = spacing(c);
    const b0 = 1 + c;                                  // A (1, 1), both detectors
    const [db1, db2] = st.db;
    const b1 = b0 + db1, b2 = b0 + db2;
    const [x1, x2] = solve2(c, b1, b2);
    dbBtns.forEach((b, i) => b.classList.toggle("active", Math.abs(dbVals[i][0] - db1) < 1e-9 && Math.abs(dbVals[i][1] - db2) < 1e-9));

    // ---------- (a) spectra ----------
    svgA.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent = "(a) emission spectra";
    svgEl("line", { x1: 58, x2: WA - 10, y1: bot, y2: bot, class: "axis" }, svgA);
    for (let w = 450; w <= 600; w += 50) {
      svgEl("line", { x1: ax(w), x2: ax(w), y1: bot, y2: bot + 5, class: "axis" }, svgA);
      svgEl("text", { x: ax(w), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgA).textContent = String(w);
    }
    for (const v of [0, 0.5, 1]) {
      if (v > 0) svgEl("line", { x1: 58, x2: WA - 10, y1: ay(v), y2: ay(v), class: "grid" }, svgA);
      svgEl("text", { x: 52, y: ay(v) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = String(v);
    }
    axisLabels(svgA, 58, WA, top, bot, "wavelength (nm)", "emission");
    const p1 = 500, p2 = 500 + d;
    const spec = (p: number, stroke: string, dash: string) => {
      const pts: string[] = [];
      for (let w = 430; w <= 640; w += 1.5) pts.push(`${ax(w)},${ay(Math.exp(-((w - p) ** 2) / (2 * WIDTH * WIDTH)))}`);
      svgEl("polygon", { points: `${ax(430)},${ay(0)} ${pts.join(" ")} ${ax(640)},${ay(0)}`, fill: stroke, "fill-opacity": 0.18, stroke: "none" }, svgA);
      svgEl("polyline", { points: pts.join(" "), fill: "none", stroke, "stroke-width": 3, "stroke-dasharray": dash }, svgA);
    };
    spec(p1, COL.d1, "");
    if (p2 < 640) spec(p2, COL.d2, "9 5");
    for (const [p, lab, anc] of [[p1, "det 1", "end"], [p2, "det 2", "start"]] as [number, string, string][]) {
      if (p > 638) continue;
      svgEl("line", { x1: ax(p), x2: ax(p), y1: ay(1.18), y2: bot, stroke: COL.grey, "stroke-width": 2, "stroke-dasharray": "3 4" }, svgA);
      svgEl("text", { x: ax(p) + (anc === "end" ? -4 : 4), y: ay(1.08), "text-anchor": anc, class: "sm" }, svgA).textContent = lab;
    }
    // c = dye 2's signal in detector 1
    svgEl("circle", { cx: ax(p1), cy: ay(c), r: 6, fill: COL.d2, stroke: COL.ink, "stroke-width": 1.5 }, svgA);
    const ct = svgEl("text", { x: ax(p1) - 9, y: ay(c) + 6, "text-anchor": "end", "font-weight": 600 }, svgA);
    Object.assign(ct.style, { paintOrder: "stroke", stroke: "#fff", strokeWidth: "4px" });
    ct.textContent = "c";
    svgEl("text", { x: 64, y: ay(1.3), class: "sm" }, svgA).textContent = `c = ${c.toFixed(2)}, peaks ${d.toFixed(0)} nm apart`;
    const lg = (txt: string, x: number, col: string) => {
      const t = svgEl("text", { x, y: ay(1.47), "font-weight": 600 }, svgA);
      t.style.fill = col; t.textContent = txt;
    };
    lg("— dye 1", 64, COL.d1);
    lg("- - dye 2", 160, COL.d2);

    // ---------- (b) detector space ----------
    layB.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, layB).textContent = "(b) detector readings b";
    const gridB = [-0.04, -0.02, 0, 0.02, 0.04];
    for (const v of gridB) {
      svgEl("line", { x1: bxs(v), x2: bxs(v), y1: top, y2: top + bSide, class: v === 0 ? "axis" : "grid", "stroke-opacity": v === 0 ? 0.35 : 1 }, layB);
      svgEl("line", { x1: bL, x2: bL + bSide, y1: bys(v), y2: bys(v), class: v === 0 ? "axis" : "grid", "stroke-opacity": v === 0 ? 0.35 : 1 }, layB);
    }
    for (const v of [-0.04, 0, 0.04]) {
      svgEl("text", { x: bxs(v), y: top + bSide + 22, "text-anchor": "middle", class: "sm" }, layB).textContent = f2(b0 + v);
      svgEl("text", { x: bL - 6, y: bys(v) + 6, "text-anchor": "end", class: "sm" }, layB).textContent = f2(b0 + v);
    }
    svgEl("rect", { x: bL, y: top, width: bSide, height: bSide, fill: "none", stroke: "#bbb" }, layB);
    axisLabels(layB, bL, bL + bSide, top, top + bSide, "detector 1", "detector 2");
    const clipB = svgEl("clipPath", { id: `um-cb-${uid}` }, layB);
    svgEl("rect", { x: bL, y: top, width: bSide, height: bSide }, clipB);
    const gB = svgEl("g", { "clip-path": `url(#um-cb-${uid})` }, layB);
    if (st.cloud && s > 0) {
      for (const [z1, z2] of Z) svgEl("circle", { cx: bxs(s * z1), cy: bys(s * z2), r: 3, fill: COL.cell, "fill-opacity": 0.5 }, gB);
      svgEl("ellipse", { cx: bxs(0), cy: bys(0), rx: bxs(2 * s) - bxs(0), ry: bxs(2 * s) - bxs(0), fill: "none", stroke: COL.grn, "stroke-width": 2.5, "stroke-dasharray": "8 5" }, gB);
    }
    svgEl("circle", { cx: bxs(0), cy: bys(0), r: 6, fill: COL.ink }, gB);
    if (db1 !== 0 || db2 !== 0) svgEl("line", { x1: bxs(0), y1: bys(0), x2: bxs(db1), y2: bys(db2), stroke: COL.bad, "stroke-width": 2.5 }, gB);
    setAttrs(handle, { cx: bxs(db1), cy: bys(db2) });
    const tb = svgEl("text", { x: bL + 6, y: top + bSide - 8, class: "sm" }, layB);
    tb.textContent = "● exact b = A·(1, 1)";

    // ---------- (c) dye space ----------
    svgC.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgC).textContent = "(c) unmixed dye amounts x = A⁻¹b";
    for (let v = 0.5; v <= 1.5 + 1e-9; v += 0.25) {
      const main = Math.abs(v - 1) < 1e-9;
      svgEl("line", { x1: cxs(v), x2: cxs(v), y1: top, y2: top + cSide, class: main ? "axis" : "grid", "stroke-opacity": main ? 0.35 : 1 }, svgC);
      svgEl("line", { x1: cL, x2: cL + cSide, y1: cys(v), y2: cys(v), class: main ? "axis" : "grid", "stroke-opacity": main ? 0.35 : 1 }, svgC);
      if (Math.abs(v * 2 - Math.round(v * 2)) < 1e-9) {
        svgEl("text", { x: cxs(v), y: top + cSide + 22, "text-anchor": "middle", class: "sm" }, svgC).textContent = String(v);
        svgEl("text", { x: cL - 6, y: cys(v) + 6, "text-anchor": "end", class: "sm" }, svgC).textContent = String(v);
      }
    }
    svgEl("rect", { x: cL, y: top, width: cSide, height: cSide, fill: "none", stroke: "#bbb" }, svgC);
    axisLabels(svgC, cL, cL + cSide, top, top + cSide, "dye 1", "dye 2");
    const clipC = svgEl("clipPath", { id: `um-cc-${uid}` }, svgC);
    svgEl("rect", { x: cL, y: top, width: cSide, height: cSide }, clipC);
    const gC = svgEl("g", { "clip-path": `url(#um-cc-${uid})` }, svgC);
    const gL = svgEl("g", {}, gC);         // lines underneath the cells
    let outside = 0;
    const X1: number[] = [];
    if (s > 0) for (const [z1, z2] of Z) {
      const [u1, u2] = solve2(c, b0 + s * z1, b0 + s * z2);
      X1.push(u1);
      if (Math.abs(u1 - 1) > XW || Math.abs(u2 - 1) > XW) outside++;
      if (st.cloud) svgEl("circle", { cx: cxs(u1), cy: cys(u2), r: 3, fill: COL.cell, "fill-opacity": 0.5 }, gC);
    }
    if (st.cloud && s > 0) {
      // predicted 2-sd ellipse: half-axes 2s/(1-c) along (1,-1), 2s/(1+c) along (1,1)
      const px = cxs(1) - cxs(0), a1 = 2 * s / (1 - c) * px, a2 = 2 * s / (1 + c) * px;
      svgEl("ellipse", { cx: cxs(1), cy: cys(1), rx: a1, ry: a2, transform: `rotate(45 ${cxs(1)} ${cys(1)})`, fill: "none", stroke: COL.grn, "stroke-width": 2.5, "stroke-dasharray": "8 5" }, gC);
    }
    // the two equations as lines: x1 + c x2 = b1 (dye-1 detector), c x1 + x2 = b2
    const lineEq = (a: number, bb: number, rhs: number, stroke: string, dash: string, w: number) => {
      if (!st.lines) return;
      // a*x1 + bb*x2 = rhs; draw across the window
      const lo = 1 - 2 * XW, hi = 1 + 2 * XW;
      let pA: [number, number], pB: [number, number];
      if (Math.abs(bb) >= Math.abs(a)) pA = [lo, (rhs - a * lo) / bb], pB = [hi, (rhs - a * hi) / bb];
      else pA = [(rhs - bb * lo) / a, lo], pB = [(rhs - bb * hi) / a, hi];
      svgEl("line", { x1: cxs(pA[0]), y1: cys(pA[1]), x2: cxs(pB[0]), y2: cys(pB[1]), stroke, "stroke-width": w, "stroke-dasharray": dash, "stroke-opacity": dash ? 1 : 0.6 }, gL);
    };
    lineEq(1, c, b0, COL.d1, "", 3);
    lineEq(c, 1, b0, COL.d2, "", 3);
    const moved = db1 !== 0 || db2 !== 0;
    if (moved) {
      lineEq(1, c, b1, COL.d1, "7 5", 2.5);
      lineEq(c, 1, b2, COL.d2, "7 5", 2.5);
      svgEl("line", { x1: cxs(1), y1: cys(1), x2: cxs(x1), y2: cys(x2), stroke: COL.bad, "stroke-width": 2.5 }, gC);
    }
    svgEl("circle", { cx: cxs(1), cy: cys(1), r: 6, fill: COL.ink }, gC);
    if (moved) svgEl("circle", { cx: cxs(x1), cy: cys(x2), r: 9, fill: COL.bad, stroke: "#fff", "stroke-width": 2.5 }, gC);
    const lgC = (txt: string, y: number, col: string) => {
      const t = svgEl("text", { x: cL + cSide + 8, y, "font-weight": 600 }, svgC);
      t.style.fill = col; t.textContent = txt;
    };
    if (st.lines && WC - cL - cSide > 60) {
      lgC("det 1", top + 40, COL.d1); lgC("line", top + 60, COL.d1);
      lgC("det 2", top + 100, COL.d2); lgC("line", top + 120, COL.d2);
    }
    if (outside > 0 && st.cloud) {
      const t = svgEl("text", { x: cL + 8, y: top + 22, class: "sm" }, svgC);
      Object.assign(t.style, { paintOrder: "stroke", stroke: "#fff", strokeWidth: "4px" });
      t.textContent = `${outside} of ${NCELL} cells outside`;
    }

    // ---------- readout ----------
    ro.replaceChildren();
    const item = (kk: string, v: string, cls = "") => {
      const sp = htmlEl("span", cls, "", ro);
      if (kk) htmlEl("b", "", kk + " ", sp);
      sp.append(v);
    };
    item("A =", `[1 ${f2(c)}; ${f2(c)} 1]`);
    item("σ =", `${f3(1 + c)}, ${f3(1 - c)}`);
    item("κ = (1+c)/(1−c) =", k >= 100 ? k.toFixed(0) : k.toFixed(1), "um-kappa");
    item("your cell: b =", `(${f3(b1)}, ${f3(b2)})`);
    item("x̂ =", `(${f3(x1)}, ${f3(x2)})`);
    if (moved) {
      const rb = Math.hypot(db1, db2) / (Math.SQRT2 * b0);
      const rx = Math.hypot(x1 - 1, x2 - 1) / Math.SQRT2;
      item("‖Δb‖/‖b‖ =", pct(rb));
      item("‖Δx‖/‖x‖ =", pct(rx));
      item("amplification", `${(rx / rb).toFixed(1)} (≤ κ)`);
    }
    if (s > 0) {
      const mean = X1.reduce((a, v) => a + v, 0) / X1.length;
      const sd = Math.sqrt(X1.reduce((a, v) => a + (v - mean) ** 2, 0) / (X1.length - 1));
      const pred = s * Math.sqrt(1 + c * c) / (1 - c * c);
      item(`cells: sd of unmixed dye 1 =`, `${f3(sd)}  (predicted ${f3(pred)}, ×${(pred / s).toFixed(1)} the noise)`);
    }
  }
  markPresets();
  draw();
}

autoMount(NAME, mount);

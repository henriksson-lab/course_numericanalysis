// Calibration-curve fitter for the least-squares-qr unit. The lecture's two standard
// curves (BCA protein assay, A562 against BSA in ug/mL; qPCR, Cq against log10 copies),
// fitted by a polynomial of degree d in one of several units for x (raw, rescaled,
// centred and scaled).
//  (a) the standards (drag them up and down; readings snap to the data's resolution),
//      the fitted curve by Householder QR, the normal-equations curve (dashed), the
//      residuals as vertical segments, and a residual strip below.
//  (b) for every degree: the relative error of the coefficients from the normal equations
//      and from QR, both computed here in doubles, against the exact least-squares
//      solution (computed in rational arithmetic with BigInt), with the lines kappa^2 eps
//      and kappa eps, kappa = kappa(A) for the chosen units. Degrees where
//      kappa(A^T A) = kappa^2 > 1/eps are shaded: there R's solve() stops.
// Auto-mounts on <div data-demo="lsqfit">. Optional attributes:
//   data-data="bca|qpcr"  data-units="raw|scaled|centred"  data-degree="2"
import { svgEl, htmlEl, slider, presets, button, scale, svgPoint, autoMount } from "../lib/ui";
import { EPS, vander, solveNormal, solveQR, singularValues, exactLS, ratToNum } from "./math";

const NAME = "lsqfit";
// Okabe-Ito, as in make-least-squares-qr-figs.R (QR blue, normal equations vermillion)
const COL = { qr: "#0072B2", ne: "#D55E00", res: "#009E73", pt: "#000000", ref: "#737373", shade: "#D55E00" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .45em; }
[data-demo="${NAME}"] .lq-row { display: flex; flex-wrap: wrap; gap: .4em 1.2em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .lq-group { display: flex; gap: .4em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .lq-row .qb-slider { min-width: 240px; flex: 0 1 300px; grid-template-columns: auto 1fr 1.6em; }
[data-demo="${NAME}"] .lq-plots { display: flex; flex-wrap: wrap; gap: 10px; width: 100%; }
[data-demo="${NAME}"] .lq-plots svg { min-width: 0; height: auto; width: calc(50% - 5px); }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .lq-read { display: flex; flex-wrap: wrap; gap: .15em 1.4em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .lq-read b { font-weight: 600; }
[data-demo="${NAME}"] .lq-warn { color: #b00020; font-weight: 600; }
@media (max-width: 760px) { [data-demo="${NAME}"] .lq-plots svg { width: 100%; } }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; gap: .3em; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; padding: .3em .5em; }
.reveal [data-demo="${NAME}"] .lq-plots svg { width: calc(50% - 5px); }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- the two data sets of the unit (x integer; y stored as integer multiples of 1/ys) ----
interface Units { key: string; label: string; m: number; h: number; desc: string }
interface DataSet {
  key: string; label: string; x: number[]; y0: number[]; ys: number; dmax: number;
  xlab: string; ylab: string; xlim: [number, number]; xticks: number[]; ylim: [number, number]; ystep: number;
  rlim: number; units: Units[]; ydec: number;
}
const BCA: DataSet = {
  key: "bca", label: "BCA protein", ys: 1000, dmax: 8, ydec: 3,
  x: [0, 25, 125, 250, 500, 750, 1000, 1500, 2000].flatMap(c => [c, c]),
  y0: [0.087, 0.103, 0.150, 0.109, 0.226, 0.229, 0.358, 0.344, 0.601, 0.569, 0.780, 0.789, 0.951, 0.942, 1.300, 1.238, 1.545, 1.532],
  xlab: "BSA (µg/mL)", ylab: "A562", xlim: [-60, 2060], xticks: [0, 500, 1000, 1500, 2000],
  ylim: [0, 1.8], ystep: 0.4, rlim: 0.12,
  units: [
    { key: "raw", label: "µg/mL", m: 0, h: 1, desc: "x in µg/mL (0 … 2000)" },
    { key: "scaled", label: "mg/mL", m: 0, h: 1000, desc: "x in mg/mL (0 … 2)" },
    { key: "centred", label: "centred, [−1, 1]", m: 1000, h: 1000, desc: "t = (x − 1000) / 1000" },
  ],
};
const QPCR: DataSet = {
  key: "qpcr", label: "qPCR", ys: 100, dmax: 6, ydec: 2,
  x: [1, 2, 3, 4, 5, 6, 7].flatMap(v => [v, v, v]),
  y0: [34.52, 34.69, 34.48, 31.62, 31.37, 31.14, 28.05, 28.10, 28.07, 24.54, 24.90, 24.68, 21.13, 20.81, 21.47, 17.89, 17.90, 18.09, 14.71, 14.67, 14.73],
  xlab: "log₁₀ copies", ylab: "Cq", xlim: [0.7, 7.3], xticks: [1, 2, 3, 4, 5, 6, 7],
  ylim: [12, 38], ystep: 4, rlim: 0.6,
  units: [
    { key: "raw", label: "log₁₀ N", m: 0, h: 1, desc: "x = log₁₀ N (1 … 7)" },
    { key: "centred", label: "centred, [−1, 1]", m: 4, h: 3, desc: "t = (log₁₀ N − 4) / 3" },
  ],
};
const SETS: DataSet[] = [BCA, QPCR];

// ---- formatting ----
const SUP: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const pow10 = (k: number): string => k === 0 ? "1" : "10" + String(k).split("").map(c => SUP[c]).join("");
function sci(v: number, d = 2): string {
  if (!Number.isFinite(v)) return "–";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 0.01 && a < 1e4) return v.toPrecision(d);
  const [mm, e] = v.toExponential(d - 1).split("e");
  return `${mm}·10${String(Number(e)).split("").map(c => SUP[c]).join("")}`;
}
const digits = (e: number): string => !Number.isFinite(e) ? "–" : e === 0 ? "≥ 16" : Math.max(0, -Math.log10(e)).toFixed(1);

interface DegResult { d: number; kappa: number; ne: number; qr: number; neX: number[] | null; qrX: number[] | null }

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const set0 = SETS.find(s => s.key === ds.data) ?? BCA;
  const unitsOf = (s: DataSet, k: string | undefined) => s.units.find(u => u.key === k) ?? s.units[Math.min(1, s.units.length - 1)];
  const init = {
    set: set0,
    units: unitsOf(set0, ds.units ?? (set0 === BCA ? "scaled" : "raw")),
    d: Math.max(1, Math.min(set0.dmax, Math.round(Number(ds.degree) || 2))),
  };
  const st = { set: init.set, units: init.units, d: init.d, yi: init.set.y0.map(v => Math.round(v * init.set.ys)) };
  const uid = Math.random().toString(36).slice(2, 8);
  const kCache = new Map<string, number>();

  // ---- controls ----
  const row1 = htmlEl("div", "lq-row", "", host);
  const g1 = htmlEl("div", "lq-group", "", row1);
  htmlEl("span", "", "standards:", g1);
  const setBtns = presets(g1, SETS.map(s => ({ label: s.label, apply: () => chooseSet(s) })));
  const g2 = htmlEl("div", "lq-group", "", row1);
  htmlEl("span", "", "x in:", g2);
  const unitBox = htmlEl("span", "lq-group", "", g2);
  let unitBtns: HTMLButtonElement[] = [];
  const row2 = htmlEl("div", "lq-row", "", host);
  const dS = slider(row2, "degree d", 1, st.set.dmax, 1, st.d, v => { st.d = v; draw(); });
  const g3 = htmlEl("div", "lq-group", "", row2);
  button(g3, "reset points", () => { st.yi = st.set.y0.map(v => Math.round(v * st.set.ys)); draw(); });
  button(g3, "reset all", () => { chooseSet(init.set, init.units, init.d); });
  htmlEl("span", "qb-note", "Drag a standard up or down. Click a degree in (b).", row2);

  const plots = htmlEl("div", "lq-plots", "", host);
  const W = 560, H = 380;
  const svgA = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Standards, fitted polynomial and residuals" }, plots);
  const svgB = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Relative error of the coefficients against degree" }, plots);
  const ro = htmlEl("div", "qb-readout lq-read", "", host);

  // panel (a) geometry
  const aL = 64, aR = W - 12, aT = 40, aB = 232, rT = 252, rB = 318;
  const layerA = svgEl("g", {}, svgA);
  const ptLayer = svgEl("g", {}, svgA);
  let circles: SVGCircleElement[] = [];
  // panel (b) geometry
  const bL = 70, bR = W - 46, bT = 40, bB = H - 62;
  const VMIN = -17, VMAX = 1;
  const by = scale(VMIN, VMAX, bB, bT);
  const layerB = svgEl("g", {}, svgB);

  // ---- data set / units choice ----
  function chooseSet(s: DataSet, u?: Units, d?: number): void {
    st.set = s;
    st.units = u ?? unitsOf(s, st.units.key);
    st.d = Math.max(1, Math.min(s.dmax, d ?? st.d));
    st.yi = s.y0.map(v => Math.round(v * s.ys));
    dS.input.max = String(s.dmax);
    dS.set(st.d);
    setBtns.forEach((b, i) => b.classList.toggle("active", SETS[i] === s));
    unitBox.replaceChildren();
    unitBtns = presets(unitBox, s.units.map(uu => ({ label: uu.label, apply: () => { st.units = uu; draw(); } })));
    buildPoints();
    draw();
  }

  // ---- draggable standards (persistent circles, so a drag survives redraws) ----
  const ax = () => scale(st.set.xlim[0], st.set.xlim[1], aL, aR);
  const ay = () => scale(st.set.ylim[0], st.set.ylim[1], aB, aT);
  function buildPoints(): void {
    ptLayer.replaceChildren();
    circles = st.set.x.map((_, i) => {
      const c = svgEl("circle", { r: 9, fill: COL.pt, "fill-opacity": 0.75, stroke: "#fff", "stroke-width": 1.5 }, ptLayer);
      c.style.cursor = "ns-resize";
      c.style.touchAction = "none";
      c.addEventListener("pointerdown", ev => {
        ev.preventDefault();
        c.setPointerCapture(ev.pointerId);
        const move = (e: PointerEvent) => {
          const [, py] = svgPoint(svgA, e);
          const s = st.set, lo = s.ylim[0], hi = s.ylim[1];
          const v = Math.max(lo, Math.min(hi, ay().invert(py)));
          const yi = Math.round(v * s.ys);
          if (yi !== st.yi[i]) { st.yi[i] = yi; draw(); }
        };
        const up = () => {
          c.removeEventListener("pointermove", move);
          c.removeEventListener("pointerup", up);
          c.removeEventListener("pointercancel", up);
        };
        c.addEventListener("pointermove", move);
        c.addEventListener("pointerup", up);
        c.addEventListener("pointercancel", up);
      });
      return c;
    });
  }

  // click a degree in panel (b)
  svgB.addEventListener("pointerdown", ev => {
    const [px, py] = svgPoint(svgB, ev);
    if (py < bT - 10 || py > bB + 30 || px < bL - 10 || px > bR + 10) return;
    const bx = scale(0.5, st.set.dmax + 0.5, bL, bR);
    const d = Math.round(bx.invert(px));
    if (d >= 1 && d <= st.set.dmax && d !== st.d) { st.d = d; dS.set(d); draw(); }
  });
  svgB.style.cursor = "pointer";

  // ---- computation ----
  function compute(): DegResult[] {
    const s = st.set, u = st.units;
    const t = s.x.map(v => (v - u.m) / u.h);
    const y = st.yi.map(v => v / s.ys);
    const ub = s.x.map(v => BigInt(v - u.m)), ybig = st.yi.map(v => BigInt(v));
    const out: DegResult[] = [];
    for (let d = 1; d <= s.dmax; d++) {
      const A = vander(t, d);
      const key = `${s.key}/${u.key}/${d}`;
      let kappa = kCache.get(key);
      if (kappa === undefined) { const sv = singularValues(A); kappa = sv[0] / sv[sv.length - 1]; kCache.set(key, kappa); }
      const ex = exactLS(ub, ybig, d).map((r, k) => ratToNum({ p: r.p * BigInt(u.h) ** BigInt(k), q: r.q * BigInt(s.ys) }));
      const nx = Math.hypot(...ex);
      const rel = (z: number[] | null) => z ? Math.hypot(...z.map((v, i) => v - ex[i])) / nx : Infinity;
      const neX = solveNormal(A, y).x, qrX = solveQR(A, y).x;
      out.push({ d, kappa, ne: rel(neX), qr: rel(qrX), neX, qrX });
    }
    return out;
  }

  const txt = (g: Element, x: number, y: number, s: string, attrs: Record<string, string | number> = {}, fill?: string) => {
    const e = svgEl("text", { x, y, ...attrs }, g);
    e.textContent = s;
    if (fill) e.style.fill = fill;
    return e;
  };

  function draw(): void {
    const s = st.set, u = st.units, d = st.d;
    unitBtns.forEach((b, i) => b.classList.toggle("active", s.units[i] === u));
    const res = compute();
    const cur = res[d - 1];
    const y = st.yi.map(v => v / s.ys);
    const X = ax(), Y = ay();
    const tOf = (xv: number) => (xv - u.m) / u.h;
    const poly = (c: number[], xv: number) => { const t = tOf(xv); let v = 0; for (let k = c.length - 1; k >= 0; k--) v = v * t + c[k]; return v; };

    // ================= (a) =================
    layerA.replaceChildren();
    txt(layerA, 8, 24, `(a) ${s.label} standards, degree ${d}`, { class: "ttl" });
    for (let v = s.ylim[0]; v <= s.ylim[1] + 1e-9; v += s.ystep) {
      svgEl("line", { x1: aL, x2: aR, y1: Y(v), y2: Y(v), class: "grid" }, layerA);
      txt(layerA, aL - 6, Y(v) + 6, String(+v.toFixed(1)), { "text-anchor": "end", class: "sm" });
    }
    for (const v of s.xticks) {
      svgEl("line", { x1: X(v), x2: X(v), y1: aT, y2: rB, class: "grid" }, layerA);
      txt(layerA, X(v), rB + 22, String(v), { "text-anchor": "middle", class: "sm" });
    }
    svgEl("line", { x1: aL, x2: aR, y1: aB, y2: aB, class: "axis" }, layerA);
    svgEl("line", { x1: aL, x2: aL, y1: aT, y2: aB, class: "axis" }, layerA);
    txt(layerA, (aL + aR) / 2, H - 12, s.xlab, { "text-anchor": "middle" });
    const ylab = txt(layerA, 16, (aT + aB) / 2, s.ylab, { "text-anchor": "middle", transform: `rotate(-90 16 ${(aT + aB) / 2})` });
    void ylab;
    const clipId = `lq-clip-${uid}`;
    const cp = svgEl("clipPath", { id: clipId }, layerA);
    svgEl("rect", { x: aL, y: aT, width: aR - aL, height: aB - aT }, cp);
    const g = svgEl("g", { "clip-path": `url(#${clipId})` }, layerA);
    const xs = Array.from({ length: 241 }, (_, i) => s.xlim[0] + (s.xlim[1] - s.xlim[0]) * i / 240);
    const curve = (c: number[]) => xs.map((xv, i) => `${i ? "L" : "M"}${X(xv).toFixed(1)},${Y(Math.max(s.ylim[0] - 10 * s.ystep, Math.min(s.ylim[1] + 10 * s.ystep, poly(c, xv)))).toFixed(1)}`).join("");
    if (cur.qrX) svgEl("path", { d: curve(cur.qrX), fill: "none", stroke: COL.qr, "stroke-width": 3 }, g);
    if (cur.neX) svgEl("path", { d: curve(cur.neX), fill: "none", stroke: COL.ne, "stroke-width": 2.5, "stroke-dasharray": "8 6" }, g);
    // residuals: vertical segments from the standard to the QR curve
    const fit = s.x.map(xv => cur.qrX ? poly(cur.qrX, xv) : NaN);
    const r = y.map((v, i) => v - fit[i]);
    s.x.forEach((xv, i) => {
      if (!Number.isFinite(fit[i])) return;
      svgEl("line", { x1: X(xv), x2: X(xv), y1: Y(y[i]), y2: Y(fit[i]), stroke: COL.res, "stroke-width": 3 }, g);
    });
    // residual strip (jitter replicates sideways so they do not hide each other)
    const R = scale(-s.rlim, s.rlim, rB, rT);
    svgEl("rect", { x: aL, y: rT, width: aR - aL, height: rB - rT, fill: "#f6f7f9" }, layerA);
    svgEl("line", { x1: aL, x2: aR, y1: R(0), y2: R(0), stroke: COL.ref, "stroke-width": 1.5 }, layerA);
    txt(layerA, aL - 6, R(0) + 6, "0", { "text-anchor": "end", class: "sm" });
    txt(layerA, aL + 6, rT + 18, `residuals (±${s.rlim})`, { class: "sm" });
    const rep: Record<number, number> = {};
    s.x.forEach((xv, i) => {
      if (!Number.isFinite(r[i])) return;
      const j = rep[xv] = (rep[xv] ?? -1) + 1;
      const xx = X(xv) + (j - (s === BCA ? 0.5 : 1)) * 7;
      const rv = Math.max(-s.rlim, Math.min(s.rlim, r[i]));
      svgEl("line", { x1: xx, x2: xx, y1: R(0), y2: R(rv), stroke: COL.res, "stroke-width": 4 }, layerA);
      if (rv !== r[i]) svgEl("circle", { cx: xx, cy: R(rv), r: 3.5, fill: COL.res }, layerA);
    });
    // legend
    const lx = s === BCA ? aR - 190 : aL + 250, ly = s === BCA ? aB - 66 : aT + 18;
    svgEl("line", { x1: lx, x2: lx + 30, y1: ly, y2: ly, stroke: COL.qr, "stroke-width": 3 }, layerA);
    txt(layerA, lx + 36, ly + 6, "QR fit", {}, COL.qr);
    svgEl("line", { x1: lx, x2: lx + 30, y1: ly + 24, y2: ly + 24, stroke: COL.ne, "stroke-width": 2.5, "stroke-dasharray": "8 6" }, layerA);
    txt(layerA, lx + 36, ly + 30, "normal eqs", {}, COL.ne);
    svgEl("line", { x1: lx + 15, x2: lx + 15, y1: ly + 40, y2: ly + 56, stroke: COL.res, "stroke-width": 3 }, layerA);
    txt(layerA, lx + 36, ly + 54, "residual", {}, COL.res);
    // standards on top
    circles.forEach((c, i) => { c.setAttribute("cx", String(X(s.x[i]))); c.setAttribute("cy", String(Y(y[i]))); });

    // ================= (b) =================
    layerB.replaceChildren();
    const bx = scale(0.5, s.dmax + 0.5, bL, bR);
    txt(layerB, 8, 24, `(b) error of the coefficients, ${u.label}`, { class: "ttl" });
    // shade degrees where kappa(A^T A) > 1/eps
    for (const q of res) if (q.kappa * q.kappa > 1 / EPS) {
      svgEl("rect", { x: bx(q.d - 0.5), y: bT, width: bx(1) - bx(0), height: bB - bT, fill: COL.shade, "fill-opacity": 0.1 }, layerB);
    }
    if (res.some(q => q.kappa * q.kappa > 1 / EPS)) {
      const first = res.find(q => q.kappa * q.kappa > 1 / EPS)!.d;
      txt(layerB, bx(first - 0.45), bB - 10, "κ(AᵀA) > 1/ε", { class: "sm" }, COL.ne);
    }
    // current degree
    svgEl("rect", { x: bx(d - 0.5), y: bT, width: bx(1) - bx(0), height: bB - bT, fill: "#000", "fill-opacity": 0.07 }, layerB);
    for (let k = VMIN + 1; k <= VMAX; k += 2) {
      svgEl("line", { x1: bL, x2: bR, y1: by(k), y2: by(k), class: "grid", "stroke-opacity": k % 4 === 0 ? 1 : 0.5 }, layerB);
      if (k % 4 === 0) txt(layerB, bL - 6, by(k) + 6, pow10(k), { "text-anchor": "end", class: "sm" });
    }
    for (let dd = 0; dd <= 16; dd += 4) txt(layerB, bR + 6, by(-dd) + 6, String(dd), { class: "sm" });
    txt(layerB, bR + 30, bT - 8, "digits", { "text-anchor": "end", class: "sm" });
    svgEl("line", { x1: bL, x2: bR, y1: bB, y2: bB, class: "axis" }, layerB);
    svgEl("line", { x1: bL, x2: bL, y1: bT, y2: bB, class: "axis" }, layerB);
    for (let dd = 1; dd <= s.dmax; dd++) txt(layerB, bx(dd), bB + 22, String(dd), { "text-anchor": "middle", class: "sm", "font-weight": dd === d ? 700 : 400 });
    txt(layerB, (bL + bR) / 2, H - 12, "polynomial degree d", { "text-anchor": "middle" });
    txt(layerB, 16, (bT + bB) / 2, "relative error", { "text-anchor": "middle", transform: `rotate(-90 16 ${(bT + bB) / 2})` });
    const cpb = svgEl("clipPath", { id: clipId + "b" }, layerB);
    svgEl("rect", { x: bL, y: bT, width: bR - bL, height: bB - bT }, cpb);
    const gb = svgEl("g", { "clip-path": `url(#${clipId}b)` }, layerB);
    const lv = (e: number) => e <= 0 ? VMIN : Math.max(VMIN, Math.min(VMAX + 1, Math.log10(e)));
    const line = (f: (q: DegResult) => number, stroke: string, dash: string, w: number) => {
      const p = res.map((q, i) => `${i ? "L" : "M"}${bx(q.d).toFixed(1)},${by(lv(f(q))).toFixed(1)}`).join("");
      svgEl("path", { d: p, fill: "none", stroke, "stroke-width": w, "stroke-dasharray": dash }, gb);
    };
    line(q => q.kappa * q.kappa * EPS, COL.ne, "7 5", 2);
    line(q => q.kappa * EPS, COL.qr, "7 5", 2);
    line(q => q.ne, COL.ne, "", 2.5);
    line(q => q.qr, COL.qr, "", 2.5);
    for (const q of res) {
      const yN = by(lv(q.ne)), yQ = by(lv(q.qr));
      svgEl("rect", { x: bx(q.d) - 7, y: yN - 7, width: 14, height: 14, fill: COL.ne, stroke: "#fff", "stroke-width": 1.5 }, gb);
      svgEl("circle", { cx: bx(q.d), cy: yQ, r: 7.5, fill: COL.qr, stroke: "#fff", "stroke-width": 1.5 }, gb);
    }
    // legend (bottom right is empty: errors are large there only for the normal equations)
    const lgx = bL + 12, lgy = bT + 18;
    const lg = (i: number, label: string, c: string, dash: string, mark: "sq" | "ci" | "") => {
      const yy = lgy + i * 24;
      svgEl("line", { x1: lgx, x2: lgx + 30, y1: yy, y2: yy, stroke: c, "stroke-width": 2.5, "stroke-dasharray": dash }, layerB);
      if (mark === "sq") svgEl("rect", { x: lgx + 8, y: yy - 7, width: 14, height: 14, fill: c }, layerB);
      if (mark === "ci") svgEl("circle", { cx: lgx + 15, cy: yy, r: 7, fill: c }, layerB);
      txt(layerB, lgx + 38, yy + 6, label, {}, c);
    };
    lg(0, "normal equations", COL.ne, "", "sq");
    lg(1, "QR (Householder)", COL.qr, "", "ci");
    lg(2, "κ²ε", COL.ne, "7 5", "");
    lg(3, "κε", COL.qr, "7 5", "");

    // ================= readout =================
    ro.replaceChildren();
    const item = (k: string, v: string, cls = "") => {
      const sp = htmlEl("span", cls, "", ro);
      if (k) htmlEl("b", "", k + " ", sp);
      sp.append(v);
    };
    const k2 = cur.kappa * cur.kappa;
    item(`d = ${d}, ${u.desc}.`, "");
    item("κ(A) =", sci(cur.kappa, 3));
    item("κ(AᵀA) = κ² =", sci(k2, 3));
    item("digits correct: normal eqs", digits(cur.ne));
    item("QR", digits(cur.qr));
    const rss = r.reduce((a, v) => a + v * v, 0);
    item("RSS =", sci(rss, 3));
    item("residual signs:", r.map(v => v > 0 ? "+" : v < 0 ? "−" : "0").join(""));
    if (k2 > 1 / EPS) item("", "κ(AᵀA) > 1/ε ≈ 4.5·10¹⁵: R's solve() stops here (\"computationally singular\")", "lq-warn");
  }

  chooseSet(init.set, init.units, init.d);
}

autoMount(NAME, mount);

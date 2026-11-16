// Runge-Kutta stepper for the runge-kutta unit (Runge-Kutta and adaptive step size).
// Tab "One step, four slopes": logistic OD curve N' = rN(1 - N/K) (K = 1.2, N0 = 0.05,
// r = 0.8 per h) or drug clearance C' = -0.8 C (C0 = 10 mg/L).
//  (a) one step from t = 0 with Euler, midpoint, Heun or RK4: every stage as a point with
//      its slope segment, the dashed line that leads to it, and the combined step (green);
//      "next stage" reveals the stages one at a time.
//  (b) global error at t = 10 against h (or against the number of f evaluations), log-log,
//      for all four methods with reference slopes 1, 2, 4, and a marker at the current h.
// Tab "Adaptive steps": an embedded pair (Heun-Euler 2(1), Bogacki-Shampine 3(2),
// Dormand-Prince 5(4)) with the unit's controller, rtol = atol from a slider, on the
// Lotka-Volterra predator-prey model (50 years) or the repressilator (t = 0 ... 50).
//  (c) the solution with a dot at every accepted step (and a tick on the time axis),
//  (d) the step size h against t (log), rejected attempts as red crosses,
//  (e) the phase plane (hares-lynx, or LacI-TetR) against a tight reference.
// Auto-mounts on <div data-demo="rkstep">. Optional attributes:
//   data-tab="step|adapt"  data-model="logistic|clearance"  data-method="euler|midpoint|heun|rk4"
//   data-h="1"  data-system="lv|rep"  data-pair="he|bs|dp"  data-logtol="-3"
//   data-adapt="1|0" (show the Adaptive tab)  data-steptab="1|0" (show the One-step tab)
import { svgEl, htmlEl, slider, button, presets, scale, autoMount } from "../lib/ui";
import {
  METHODS, stepStages, fixedSolve, LOGISTIC, CLEARANCE, Scalar,
  PAIRS, adapt, AdaptResult, lvF, lvV, repF, REP_P1_50, F, Vec,
} from "./math";

const NAME = "rkstep";
// Okabe-Ito, as in make-runge-kutta-figs.R
const OI = {
  orange: "#E69F00", sky: "#56B4E9", green: "#009E73", yellow: "#F0E442",
  blue: "#0072B2", red: "#D55E00", purple: "#CC79A7", grey: "#999999",
};
const STAGE_COL = [OI.red, OI.orange, OI.purple, OI.blue];     // k1..k4 as in fig03
const METH_COL: Record<string, string> = { euler: OI.red, midpoint: OI.sky, heun: OI.orange, rk4: OI.blue };
const INK = "#222", REF = "#737373";

const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .rk-row { display: flex; flex-wrap: wrap; gap: .35em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .rk-group { display: flex; gap: .35em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .rk-row .qb-slider { min-width: 300px; flex: 1 1 300px; max-width: 440px; grid-template-columns: auto 1fr 4.6em; }
[data-demo="${NAME}"] .rk-plots { display: flex; flex-wrap: wrap; gap: 10px; width: 100%; }
[data-demo="${NAME}"] .rk-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 19px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .rk-read { display: flex; flex-wrap: wrap; gap: .15em 1.3em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .rk-read span b { font-weight: 600; }
[data-demo="${NAME}"] .rk-warn { color: #b00020; font-weight: 600; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; gap: .25em; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; padding: .3em .5em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- formatting ----
const sup: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const pow10 = (k: number): string => k === 0 ? "1" : "10" + String(k).split("").map(c => sup[c]).join("");
function sci(v: number, d = 2): string {
  if (!Number.isFinite(v)) return String(v);
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 0.01 && a < 1000) return v.toPrecision(d);
  const [m, e] = v.toExponential(d - 1).split("e");
  return `${m}e${Number(e)}`;
}
const fmtInt = (n: number): string => Math.round(n).toString().replace(/\B(?=(\d{3})+(?!\d))/g, " ");
const tolFmt = (l: number): string => Math.abs(l - Math.round(l)) < 1e-9 ? `1e${Math.round(l)}` : sci(10 ** l, 2);

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const uid = ++mountCount;
  const withStep = ds.steptab !== "0", withAdapt = ds.adapt !== "0";
  const stepBox = htmlEl("div", "", ""), adaptBox = htmlEl("div", "", "");
  for (const b of [stepBox, adaptBox]) Object.assign(b.style, { display: "flex", flexDirection: "column", gap: ".4em", width: "100%" });
  if (withStep && withAdapt) {
    const tabs = htmlEl("div", "rk-row", "", host);
    const tb = presets(tabs, [
      { label: "One step, four slopes", apply: () => show("step") },
      { label: "Adaptive steps", apply: () => show("adapt") },
    ]);
    const show = (t: string) => {
      stepBox.style.display = t === "step" ? "flex" : "none";
      adaptBox.style.display = t === "adapt" ? "flex" : "none";
      tb.forEach((b, i) => b.classList.toggle("active", (i === 0) === (t === "step")));
    };
    host.append(stepBox, adaptBox);
    show(ds.tab === "adapt" ? "adapt" : "step");
    buildStep(stepBox, ds, uid);
    buildAdapt(adaptBox, ds, uid);
  } else if (withAdapt) { host.append(adaptBox); buildAdapt(adaptBox, ds, uid); }
  else { host.append(stepBox); buildStep(stepBox, ds, uid); }
}

function axisLabels(svg: SVGSVGElement, xl: string, yl: string, x0: number, x1: number, top: number, bot: number, H: number): void {
  svgEl("text", { x: (x0 + x1) / 2, y: H - 6, "text-anchor": "middle" }, svg).textContent = xl;
  const t = svgEl("text", { x: 18, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 18 ${(top + bot) / 2})` }, svg);
  t.textContent = yl;
}
function colText(svg: Element, x: number, y: number, txt: string, c: string, anchor = "start", bold = true): SVGTextElement {
  const t = svgEl("text", { x, y, "text-anchor": anchor, "font-weight": bold ? 600 : 400 }, svg);
  t.style.fill = c; t.textContent = txt;
  return t;
}
/** Round-number ticks for a linear axis. */
function niceTicks(lo: number, hi: number, n = 5): number[] {
  const raw = (hi - lo) / n, p = 10 ** Math.floor(Math.log10(raw));
  const st = [1, 2, 2.5, 5, 10].map(m => m * p).find(s => s >= raw) ?? 10 * p;
  const out: number[] = [];
  for (let v = Math.ceil(lo / st) * st; v <= hi + 1e-12; v += st) out.push(Math.abs(v) < 1e-12 ? 0 : v);
  return out;
}
const tickFmt = (v: number): string => Number(v.toPrecision(6)).toString();

// =================================================================== tab 1: one step
const MKEYS = ["euler", "midpoint", "heun", "rk4"];
const HS = Array.from({ length: 10 }, (_, i) => 2 ** (1 - i));   // 2 ... 1/256

function buildStep(box: HTMLElement, ds: DOMStringMap, uid: number): void {
  const start = {
    model: ds.model === "clearance" ? CLEARANCE : LOGISTIC,
    meth: MKEYS.includes(ds.method ?? "") ? ds.method as string : "rk4",
    h: Number(ds.h) >= 0.05 && Number(ds.h) <= 2 ? Number(ds.h) : 1,
  };
  const st = { model: start.model as Scalar, meth: start.meth, h: start.h, shown: 99, xEval: false };

  // error curves per model, cached
  const curves = new Map<string, Record<string, number[]>>();
  const getCurves = (m: Scalar) => {
    let c = curves.get(m.key);
    if (!c) {
      c = {};
      for (const k of MKEYS) c[k] = HS.map(h => Math.abs(fixedSolve(METHODS[k], m.f, 0, m.y0, h, Math.round(m.T / h)) - m.exact(m.T)));
      curves.set(m.key, c);
    }
    return c;
  };

  const r1 = htmlEl("div", "rk-row", "", box);
  const g1 = htmlEl("div", "rk-group", "", r1);
  htmlEl("span", "", "model:", g1);
  const modBtns = presets(g1, [LOGISTIC, CLEARANCE].map(m => ({ label: m.label, apply: () => { st.model = m; draw(); } })));
  const g2 = htmlEl("div", "rk-group", "", r1);
  htmlEl("span", "", "method:", g2);
  const methBtns = presets(g2, MKEYS.map(k => ({ label: METHODS[k].name, apply: () => { st.meth = k; st.shown = 99; draw(); } })));
  const r2 = htmlEl("div", "rk-row", "", box);
  const hS = slider(r2, "step h", 0.05, 2, 0.05, st.h, v => { st.h = v; draw(); });
  const g3 = htmlEl("div", "rk-group", "", r2);
  button(g3, "next stage", () => {
    const s = METHODS[st.meth].c.length;
    st.shown = st.shown >= s ? 1 : st.shown + 1; draw();
  });
  button(g3, "all stages", () => { st.shown = 99; draw(); });
  const g4 = htmlEl("div", "rk-group", "", r2);
  htmlEl("span", "", "(b) against:", g4);
  const xBtns = presets(g4, [
    { label: "h", apply: () => { st.xEval = false; draw(); } },
    { label: "f evaluations", apply: () => { st.xEval = true; draw(); } },
  ]);
  button(r2, "reset", () => {
    st.model = start.model; st.meth = start.meth; st.h = start.h; st.shown = 99; st.xEval = false;
    hS.set(st.h); draw();
  });

  const plots = htmlEl("div", "rk-plots", "", box);
  const H = 360, WA = 520, WB = 600, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "One Runge-Kutta step with its stages" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Global error against step size, log-log" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro = htmlEl("div", "qb-readout rk-read", "", box);
  const top = 36, bot = H - 50;

  function draw(): void {
    const m = st.model, tb = METHODS[st.meth], h = st.h, s = tb.c.length;
    const shown = Math.min(st.shown, s);
    modBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (m.key === "logistic")));
    methBtns.forEach((b, i) => b.classList.toggle("active", MKEYS[i] === st.meth));
    xBtns.forEach((b, i) => b.classList.toggle("active", (i === 1) === st.xEval));
    const res = stepStages(tb, m.f, 0, m.y0, h);
    const yEx = m.exact(h);

    // ---------- (a) one step ----------
    svgA.replaceChildren();
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgA).textContent = `(a) one ${tb.name} step, h = ${h.toFixed(2)}`;
    const L = 80, Rr = WA - 14;
    const ts = Array.from({ length: 81 }, (_, i) => -0.08 * h + i * 1.18 * h / 80);
    const ys = [...ts.map(m.exact), ...res.stages.map(q => q.y), res.y1];
    let ylo = Math.min(...ys), yhi = Math.max(...ys);
    const pad = yhi - ylo || 1; ylo -= 0.22 * pad; yhi += 0.12 * pad;
    const PT = top + 30;
    const ax = scale(ts[0], ts[ts.length - 1], L, Rr), ay = scale(ylo, yhi, bot, PT);
    for (const v of niceTicks(ylo, yhi, 5)) {
      svgEl("line", { x1: L, x2: Rr, y1: ay(v), y2: ay(v), class: "grid" }, svgA);
      svgEl("text", { x: L - 6, y: ay(v) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = tickFmt(v);
    }
    svgEl("line", { x1: L, x2: Rr, y1: bot, y2: bot, class: "axis" }, svgA);
    for (const [v, lab] of [[0, "0"], [h / 2, "h/2"], [h, "h"]] as [number, string][]) {
      svgEl("line", { x1: ax(v), x2: ax(v), y1: PT, y2: bot + 5, class: "grid" }, svgA);
      svgEl("text", { x: ax(v), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgA).textContent = lab;
    }
    axisLabels(svgA, `time t (${m.tunit})`, m.key === "logistic" ? "OD  N(t)" : "C(t) (mg/L)", L, Rr, PT, bot, H);
    const cp = svgEl("clipPath", { id: `rk-clipa-${uid}` }, svgA);
    svgEl("rect", { x: L, y: PT - 6, width: Rr - L, height: bot - PT + 6 }, cp);
    const g = svgEl("g", { "clip-path": `url(#rk-clipa-${uid})` }, svgA);
    svgEl("polyline", { points: ts.map(t => `${ax(t)},${ay(m.exact(t))}`).join(" "), fill: "none", stroke: OI.grey, "stroke-width": 5 }, g);
    const half = 0.13 * h;
    res.stages.forEach((q, i) => {
      if (i >= shown) return;
      const c = STAGE_COL[i] ?? INK;
      if (i > 0) {   // the line that leads from (t0, y0) to the stage point
        const viaCol = tb.A[i].filter(a => a !== 0).length === 1 ? STAGE_COL[tb.A[i].findIndex(a => a !== 0)] : REF;
        svgEl("line", { x1: ax(0), y1: ay(m.y0), x2: ax(q.t), y2: ay(q.y), stroke: viaCol, "stroke-width": 2, "stroke-dasharray": "5 5" }, g);
      }
      svgEl("line", { x1: ax(q.t - half), y1: ay(q.y - half * q.k), x2: ax(q.t + half), y2: ay(q.y + half * q.k), stroke: c, "stroke-width": 5, "stroke-linecap": "round" }, g);
      svgEl("circle", { cx: ax(q.t), cy: ay(q.y), r: 7, fill: c, stroke: "#fff", "stroke-width": 1.5 }, g);
      // k1 right-below, k2 below, k3 above, k4 left-above (k3's point lies above k2's in both models)
      const dx = [14, 0, 0, -14][i] ?? 0, dy = [30, 34, -16, -14][i] ?? 0;
      const lab = colText(g, ax(q.t) + dx, ay(q.y) + dy, `k${i + 1}`, c, ["start", "middle", "middle", "end"][i] ?? "start");
      lab.style.fontSize = "20px";
    });
    if (shown >= s) {
      svgEl("line", { x1: ax(0), y1: ay(m.y0), x2: ax(h), y2: ay(res.y1), stroke: OI.green, "stroke-width": 3.5, "stroke-dasharray": "10 6" }, g);
      svgEl("rect", { x: ax(h) - 8, y: ay(res.y1) - 8, width: 16, height: 16, fill: OI.green }, g);
    }
    svgEl("circle", { cx: ax(h), cy: ay(yEx), r: 8, fill: "none", stroke: INK, "stroke-width": 2.5 }, g);
    // legend
    colText(svgA, 10, top + 12, "— exact", "#888");
    colText(svgA, 100, top + 12, "■ new value", OI.green);
    colText(svgA, 240, top + 12, "○ exact at t = h", INK, "start", false);
    if (shown < s) colText(svgA, Rr, top + 12, `stage ${shown} of ${s}`, INK, "end");

    // ---------- (b) order ----------
    svgB.replaceChildren();
    const cur = getCurves(m);
    const xE = st.xEval;
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgB).textContent =
      `(b) error in ${m.yname}(${m.T}) against ${xE ? "f evaluations" : "step h"} (log–log)`;
    const BL = 74, BR = WB - 14, EY0 = -15, EY1 = 0;
    const X0 = xE ? Math.log10(4) : Math.log10(HS[HS.length - 1]) - 0.1, X1 = xE ? Math.log10(2e4) : Math.log10(2) + 0.1;
    const bx = scale(X0, X1, BL, BR), by = scale(EY0, EY1, bot, top);
    const ly = (e: number) => Math.max(EY0 - 1, Math.log10(Math.max(e, 1e-30)));
    for (let k = EY0; k <= EY1; k += 3) {
      svgEl("line", { x1: BL, x2: BR, y1: by(k), y2: by(k), class: "grid" }, svgB);
      svgEl("text", { x: BL - 6, y: by(k) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = pow10(k);
    }
    svgEl("line", { x1: BL, x2: BR, y1: bot, y2: bot, class: "axis" }, svgB);
    const xt: [number, string][] = xE ? [1, 2, 3, 4].map(k => [k, pow10(k)])
      : [[0, "1"], [-1, "0.1"], [-2, "0.01"]];
    if (!xE) for (const k of [1, 0.5, 0.25, 0.125, 1 / 16, 1 / 64, 1 / 256]) svgEl("line", { x1: bx(Math.log10(k)), x2: bx(Math.log10(k)), y1: bot, y2: bot + 5, class: "axis" }, svgB);
    for (const [k, lab] of xt) {
      svgEl("line", { x1: bx(k), x2: bx(k), y1: top, y2: bot + 5, class: "grid" }, svgB);
      svgEl("text", { x: bx(k), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgB).textContent = lab;
    }
    axisLabels(svgB, xE ? "f evaluations from 0 to 10" : `step size h (${m.tunit})`, "absolute error", BL, BR, top, bot, H);
    const cpb = svgEl("clipPath", { id: `rk-clipb-${uid}` }, svgB);
    svgEl("rect", { x: BL, y: top - 4, width: BR - BL, height: bot - top + 8 }, cpb);
    const gb = svgEl("g", { "clip-path": `url(#rk-clipb-${uid})` }, svgB);
    const xOf = (k: string, h0: number) => xE ? Math.log10(Math.round(m.T / h0) * METHODS[k].c.length) : Math.log10(h0);
    // reference slopes, anchored near the RK4/Heun/Euler curves at h = 1/8
    for (const [k, p] of [["euler", 1], ["heun", 2], ["rk4", 4]] as [string, number][]) {
      const i0 = 4, e0 = cur[k][i0] * 4;
      const pts: string[] = [];
      for (let i = 1; i <= 8; i++) {
        const hh = HS[i];
        const e = e0 * (hh / HS[i0]) ** p;
        pts.push(`${bx(xOf(k, hh))},${by(ly(e))}`);
      }
      svgEl("polyline", { points: pts.join(" "), fill: "none", stroke: REF, "stroke-width": 2, "stroke-dasharray": "7 5" }, gb);
      const hl = HS[8], el = e0 * (hl / HS[i0]) ** p;
      const tx = svgEl("text", { x: bx(xOf(k, hl)) + (xE ? 6 : -6), y: by(ly(el)) + (xE ? 18 : 18), "text-anchor": xE ? "start" : "end", class: "sm" }, gb);
      tx.textContent = `slope ${xE ? "−" : ""}${p}`;
    }
    if (xE) {
      svgEl("line", { x1: BL, x2: BR, y1: by(-6), y2: by(-6), stroke: REF, "stroke-width": 1.5, "stroke-dasharray": "2 4" }, gb);
      colText(gb, BL + 6, by(-6) - 6, "10⁻⁶", REF, "start", false);
    }
    for (const k of MKEYS) {
      if (k === st.meth) continue;
      drawCurve(k, false);
    }
    drawCurve(st.meth, true);
    function drawCurve(k: string, strong: boolean): void {
      const pts = HS.map((hh, i) => `${bx(xOf(k, hh))},${by(ly(cur[k][i]))}`).join(" ");
      svgEl("polyline", { points: pts, fill: "none", stroke: METH_COL[k], "stroke-width": strong ? 3.5 : 2, "stroke-opacity": strong ? 1 : 0.45 }, gb);
      HS.forEach((hh, i) => svgEl("circle", { cx: bx(xOf(k, hh)), cy: by(ly(cur[k][i])), r: strong ? 4.5 : 3, fill: METH_COL[k], "fill-opacity": strong ? 1 : 0.45 }, gb));
    }
    // current h
    const n = Math.max(1, Math.round(m.T / h)), hEff = m.T / n;
    const eNow = Math.abs(fixedSolve(tb, m.f, 0, m.y0, hEff, n) - m.exact(m.T));
    const eHalf = Math.abs(fixedSolve(tb, m.f, 0, m.y0, hEff / 2, 2 * n) - m.exact(m.T));
    const xm = bx(xE ? Math.log10(n * s) : Math.log10(hEff));
    svgEl("circle", { cx: xm, cy: by(ly(eNow)), r: 10, fill: METH_COL[st.meth], stroke: INK, "stroke-width": 2.5 }, gb);
    // legend
    MKEYS.forEach((k, i) => {
      const y = xE ? bot - 104 + 26 * i : bot - 104 + 26 * i, x = xE ? BL + 14 : BR - 130;
      svgEl("line", { x1: x, x2: x + 24, y1: y - 6, y2: y - 6, stroke: METH_COL[k], "stroke-width": 4, "stroke-opacity": k === st.meth ? 1 : 0.5 }, svgB);
      const t = colText(svgB, x + 30, y, METHODS[k].name, METH_COL[k]);
      t.style.fillOpacity = k === st.meth ? "1" : "0.6";
    });

    // ---------- readout ----------
    ro.replaceChildren();
    const item = (k: string, v: string, cls = "") => {
      const sp = htmlEl("span", cls, "", ro);
      if (k) htmlEl("b", "", k + " ", sp);
      sp.append(v);
    };
    const f6 = (v: number) => v.toPrecision(7);
    res.stages.forEach((q, i) => item(`k${i + 1} = f(${Number(q.t.toFixed(4))}, ${f6(q.y)}) =`, f6(q.k)));
    const wts = s === 4 ? "(k1 + 2k2 + 2k3 + k4)/6" : st.meth === "heun" ? "(k1 + k2)/2" : st.meth === "midpoint" ? "k2" : "k1";
    item(`${m.yname}₁ = ${m.yname}₀ + h·${wts} =`, f6(res.y1));
    item("exact", f6(yEx));
    item("error of this step", sci(res.y1 - yEx, 2));
    const eul = stepStages(METHODS.euler, m.f, 0, m.y0, h).y1 - yEx;
    if (st.meth !== "euler") item("(Euler:", sci(eul, 2) + ")");
    item(`to t = ${m.T}:`, `${n} steps, ${n * s} f evaluations, error ${sci(eNow, 2)}`);
    item("halving h divides it by", eHalf > 0 ? `${(eNow / eHalf).toFixed(1)} = 2^${Math.log2(eNow / eHalf).toFixed(2)}` : "–");
  }
  draw();
}

// =================================================================== tab 2: adaptive steps
interface Sys {
  key: string; label: string; f: F; y0: Vec; T: number; h0: number; tunit: string;
  show: number[]; names: string[]; cols: string[]; phase: [number, number]; phaseNames: [string, string];
  ref?: AdaptResult;
}
const SYSTEMS: Record<string, Sys> = {
  lv: {
    key: "lv", label: "predator–prey", f: lvF, y0: [10, 5], T: 50, h0: 0.01, tunit: "years",
    show: [0, 1], names: ["hares x", "lynx y"], cols: [OI.blue, OI.red], phase: [0, 1], phaseNames: ["hares x", "lynx y"],
  },
  rep: {
    key: "rep", label: "repressilator", f: repF, y0: [1, 0, 0, 0, 0, 0], T: 50, h0: 0.01, tunit: "mRNA lifetimes",
    show: [3, 4, 5], names: ["LacI p₁", "TetR p₂", "cI p₃"], cols: [OI.blue, OI.red, OI.green], phase: [3, 4], phaseNames: ["LacI p₁", "TetR p₂"],
  },
};
const refOf = (s: Sys): AdaptResult => s.ref ?? (s.ref = adapt(PAIRS.dp, s.f, 0, s.y0, s.T, 1e-10, 1e-10, s.h0));

function buildAdapt(box: HTMLElement, ds: DOMStringMap, uid: number): void {
  const start = {
    sys: ds.system === "rep" ? "rep" : "lv",
    pair: ds.pair && PAIRS[ds.pair] ? ds.pair : "dp",
    lt: Number.isFinite(Number(ds.logtol)) && ds.logtol !== undefined ? Math.min(-2, Math.max(-9, Number(ds.logtol))) : -3,
  };
  const st = { ...start };

  const r1 = htmlEl("div", "rk-row", "", box);
  const g1 = htmlEl("div", "rk-group", "", r1);
  htmlEl("span", "", "model:", g1);
  const sysKeys = ["lv", "rep"];
  const sysBtns = presets(g1, sysKeys.map(k => ({ label: SYSTEMS[k].label, apply: () => { st.sys = k; draw(); } })));
  const g2 = htmlEl("div", "rk-group", "", r1);
  htmlEl("span", "", "pair:", g2);
  const pairKeys = ["he", "bs", "dp"];
  const pairBtns = presets(g2, pairKeys.map(k => ({ label: PAIRS[k].name, apply: () => { st.pair = k; draw(); } })));
  const r2 = htmlEl("div", "rk-row", "", box);
  const tS = slider(r2, "rtol = atol", -9, -2, 0.5, st.lt, v => { st.lt = v; draw(); }, tolFmt);
  const g3 = htmlEl("div", "rk-group", "", r2);
  const tolVals = [-2, -3, -6, -9];
  const tolBtns = presets(g3, tolVals.map(l => ({ label: `1e${l}`, apply: () => { st.lt = l; tS.set(l); draw(); } })));
  button(r2, "reset", () => { Object.assign(st, start); tS.set(st.lt); draw(); });

  const plots = htmlEl("div", "rk-plots", "", box);
  const H = 380, WA = 720, WB = 400, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Adaptive solution and step sizes" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Phase plane" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro = htmlEl("div", "qb-readout rk-read", "", box);

  const cache = new Map<string, AdaptResult>();
  function solve(sys: Sys, pk: string, lt: number): AdaptResult {
    const key = `${sys.key}|${pk}|${lt}`;
    let r = cache.get(key);
    if (!r) {
      const tol = 10 ** lt;
      r = adapt(PAIRS[pk], sys.f, 0, sys.y0, sys.T, tol, tol, sys.h0, 150000);
      if (cache.size > 40) cache.clear();
      cache.set(key, r);
    }
    return r;
  }

  function draw(): void {
    const sys = SYSTEMS[st.sys], pair = PAIRS[st.pair];
    sysBtns.forEach((b, i) => b.classList.toggle("active", sysKeys[i] === st.sys));
    pairBtns.forEach((b, i) => b.classList.toggle("active", pairKeys[i] === st.pair));
    tolBtns.forEach((b, i) => b.classList.toggle("active", Math.abs(tolVals[i] - st.lt) < 1e-9));
    const ref = refOf(sys), r = solve(sys, st.pair, st.lt);
    const nAcc = r.h.length;
    const thin = Math.max(1, Math.ceil(nAcc / 2500));   // draw at most ~2500 dots

    // ---------- (c) solution, (d) step sizes ----------
    svgA.replaceChildren();
    const L = 70, Rr = WA - 12;
    const cTop = 34, cBot = 196, dTop = 250, dBot = H - 44;
    const tx = scale(0, sys.T, L, Rr);
    let ymax = 0;
    for (const y of ref.y) for (const j of sys.show) ymax = Math.max(ymax, y[j]);
    ymax *= 1.08;
    const cy = scale(0, ymax, cBot, cTop + 10);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgA).textContent = `(c) solution: a dot at each accepted step`;
    for (const v of niceTicks(0, ymax, 3)) {
      svgEl("line", { x1: L, x2: Rr, y1: cy(v), y2: cy(v), class: "grid" }, svgA);
      svgEl("text", { x: L - 6, y: cy(v) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = tickFmt(v);
    }
    svgEl("line", { x1: L, x2: Rr, y1: cBot, y2: cBot, class: "axis" }, svgA);
    const cpa = svgEl("clipPath", { id: `rk-clipc-${uid}` }, svgA);
    svgEl("rect", { x: L, y: cTop - 2, width: Rr - L, height: cBot - cTop + 2 }, cpa);
    const gc = svgEl("g", { "clip-path": `url(#rk-clipc-${uid})` }, svgA);
    const clampC = (v: number) => cy(Math.max(-0.05 * ymax, Math.min(1.3 * ymax, v)));
    sys.show.forEach((j, ii) => {
      const pr = ref.t.map((t, i) => `${tx(t).toFixed(1)},${cy(ref.y[i][j]).toFixed(1)}`).join(" ");
      svgEl("polyline", { points: pr, fill: "none", stroke: "#cfcfcf", "stroke-width": 6 }, gc);
      const ps: string[] = [];
      for (let i = 0; i < r.t.length; i += thin) ps.push(`${tx(r.t[i]).toFixed(1)},${clampC(r.y[i][j]).toFixed(1)}`);
      svgEl("polyline", { points: ps.join(" "), fill: "none", stroke: sys.cols[ii], "stroke-width": 1.6 }, gc);
      if (ii === 0 || sys.key === "lv") for (let i = 0; i < r.t.length; i += thin)
        svgEl("circle", { cx: tx(r.t[i]).toFixed(1), cy: clampC(r.y[i][j]).toFixed(1), r: 2.8, fill: sys.cols[ii] }, gc);
    });
    // accepted-step ticks along the time axis
    const tk: string[] = [];
    for (let i = 0; i < r.t.length; i += thin) { const x = tx(r.t[i]).toFixed(1); tk.push(`M${x},${cBot}V${cBot + 9}`); }
    svgEl("path", { d: tk.join(""), stroke: INK, "stroke-width": 1 }, svgA);
    // legend
    let lx = L + 10;
    sys.show.forEach((_, ii) => {
      const t = colText(svgA, lx, cTop + 14, "● " + sys.names[ii], sys.cols[ii]);
      lx += 22 + 10.5 * sys.names[ii].length + 14;
      void t;
    });
    colText(svgA, lx, cTop + 14, "▬ reference", "#aaa");

    // (d) step size
    const hs = r.h, rh = r.rejH;
    let hmin = Infinity, hmax = 0;    // loops, not Math.min(...): the arrays can hold 300 000 steps
    for (const v of hs) { hmin = Math.min(hmin, v); hmax = Math.max(hmax, v); }
    const amin = rh.reduce((a, v) => Math.min(a, v), hmin), amax = rh.reduce((a, v) => Math.max(a, v), hmax);
    let lo = Math.floor(Math.log10(amin)), hi = Math.ceil(Math.log10(amax));
    if (hi - lo < 2) { lo -= 1; hi += 1; }
    const hy = scale(lo, hi, dBot, dTop);
    svgEl("text", { x: 8, y: dTop - 12, class: "ttl" }, svgA).textContent = "(d) step size h (log axis); × = rejected";
    const every = hi - lo > 3 ? 2 : 1;
    for (let k = lo; k <= hi; k++) {
      svgEl("line", { x1: L, x2: Rr, y1: hy(k), y2: hy(k), class: "grid" }, svgA);
      if ((hi - k) % every === 0) svgEl("text", { x: L - 6, y: hy(k) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = pow10(k);
    }
    svgEl("line", { x1: L, x2: Rr, y1: dBot, y2: dBot, class: "axis" }, svgA);
    for (const v of niceTicks(0, sys.T, 5)) {
      svgEl("line", { x1: tx(v), x2: tx(v), y1: dBot, y2: dBot + 5, class: "axis" }, svgA);
      svgEl("text", { x: tx(v), y: dBot + 22, "text-anchor": "middle", class: "sm" }, svgA).textContent = tickFmt(v);
    }
    svgEl("text", { x: (L + Rr) / 2, y: H - 4, "text-anchor": "middle" }, svgA).textContent = `time t (${sys.tunit})`;
    const hp: string[] = [];
    for (let i = 0; i < hs.length; i += thin) hp.push(`${tx(r.t[i]).toFixed(1)},${hy(Math.log10(hs[i])).toFixed(1)}`);
    svgEl("polyline", { points: hp.join(" "), fill: "none", stroke: OI.orange, "stroke-width": 2 }, svgA);
    if (hs.length / thin < 400) for (let i = 0; i < hs.length; i += thin)
      svgEl("circle", { cx: tx(r.t[i]).toFixed(1), cy: hy(Math.log10(hs[i])).toFixed(1), r: 3, fill: OI.orange }, svgA);
    const rx: string[] = [];
    for (let i = 0; i < r.rejT.length; i++) {
      const x = tx(r.rejT[i]), y = hy(Math.log10(r.rejH[i]));
      rx.push(`M${x - 6},${y - 6}L${x + 6},${y + 6}M${x - 6},${y + 6}L${x + 6},${y - 6}`);
    }
    if (rx.length) svgEl("path", { d: rx.join(""), stroke: "#b00020", "stroke-width": 2.5 }, svgA);
    for (let i = 0; i < r.rejT.length; i++) {
      const x = tx(r.rejT[i]);
      svgEl("line", { x1: x, x2: x, y1: cBot, y2: cBot + 9, stroke: "#b00020", "stroke-width": 2.5 }, svgA);
    }

    // ---------- (e) phase plane ----------
    svgB.replaceChildren();
    const [px, py] = sys.phase;
    let xm = 0, ym = 0;
    for (const y of ref.y) { xm = Math.max(xm, y[px]); ym = Math.max(ym, y[py]); }
    xm *= 1.12; ym *= 1.12;
    const PL = 66, PR = WB - 14, PT = 34, PB = H - 50;
    const sx = scale(0, xm, PL, PR), sy = scale(0, ym, PB, PT);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgB).textContent = "(e) phase plane";
    for (const v of niceTicks(0, ym, 4)) {
      svgEl("line", { x1: PL, x2: PR, y1: sy(v), y2: sy(v), class: "grid" }, svgB);
      svgEl("text", { x: PL - 6, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = tickFmt(v);
    }
    for (const v of niceTicks(0, xm, 4)) {
      svgEl("line", { x1: sx(v), x2: sx(v), y1: PT, y2: PB + 5, class: "grid" }, svgB);
      svgEl("text", { x: sx(v), y: PB + 22, "text-anchor": "middle", class: "sm" }, svgB).textContent = tickFmt(v);
    }
    svgEl("line", { x1: PL, x2: PR, y1: PB, y2: PB, class: "axis" }, svgB);
    axisLabels(svgB, sys.phaseNames[0], sys.phaseNames[1], PL, PR, PT, PB, H);
    const cpe = svgEl("clipPath", { id: `rk-clipe-${uid}` }, svgB);
    svgEl("rect", { x: PL, y: PT, width: PR - PL, height: PB - PT }, cpe);
    const ge = svgEl("g", { "clip-path": `url(#rk-clipe-${uid})` }, svgB);
    // reference: for LV one closed orbit; for the repressilator the cycle after the transient
    const refPts: string[] = [];
    for (let i = 0; i < ref.t.length; i++) {
      if (sys.key === "lv" && ref.t[i] > 6) break;
      if (sys.key === "rep" && ref.t[i] < 25) continue;
      refPts.push(`${sx(ref.y[i][px]).toFixed(1)},${sy(ref.y[i][py]).toFixed(1)}`);
    }
    svgEl("polyline", { points: refPts.join(" "), fill: "none", stroke: "#cfcfcf", "stroke-width": 7 }, ge);
    const pp: string[] = [];
    for (let i = 0; i < r.t.length; i += thin) pp.push(`${sx(r.y[i][px]).toFixed(1)},${sy(r.y[i][py]).toFixed(1)}`);
    svgEl("polyline", { points: pp.join(" "), fill: "none", stroke: OI.blue, "stroke-width": 1.4 }, ge);
    if (r.t.length / thin < 600) for (let i = 0; i < r.t.length; i += thin)
      svgEl("circle", { cx: sx(r.y[i][px]).toFixed(1), cy: sy(r.y[i][py]).toFixed(1), r: 2.6, fill: OI.blue }, ge);
    svgEl("circle", { cx: sx(sys.y0[px]), cy: sy(sys.y0[py]), r: 6, fill: INK }, ge);

    // ---------- readout ----------
    ro.replaceChildren();
    const item = (k: string, v: string, cls = "") => {
      const sp = htmlEl("span", cls, "", ro);
      if (k) htmlEl("b", "", k + " ", sp);
      sp.append(v);
    };
    item(pair.name + ", rtol = atol =", tolFmt(st.lt));
    item("accepted steps", fmtInt(nAcc));
    item("rejected", fmtInt(r.rejT.length));
    item("f evaluations", fmtInt(r.nfev));
    item("h from", `${sci(hmin, 2)} to ${sci(hmax, 2)}`);
    const yEnd = r.y[r.y.length - 1];
    if (sys.key === "lv") {
      const V0 = lvV(sys.y0[0], sys.y0[1]);
      let dv = 0;
      for (const y of r.y) dv = Math.max(dv, Math.abs(lvV(y[0], y[1]) - V0));
      item("conserved V: max |V(t) − V(0)| =", sci(dv, 2));
    } else {
      item("p₁(50) =", `${yEnd[3].toFixed(4)} (reference ${REP_P1_50}, error ${sci(yEnd[3] - REP_P1_50, 2)})`);
    }
    if (r.stopped) item("", `stopped after ${fmtInt(nAcc + r.rejT.length)} attempts at t = ${r.t[r.t.length - 1].toFixed(2)}: too many steps for this pair`, "rk-warn");
  }
  draw();
}

autoMount(NAME, mount);

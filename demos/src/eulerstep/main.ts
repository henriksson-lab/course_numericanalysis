// Euler stepper for the euler unit (Initial value problems and Euler's method).
// Tab "One equation": logistic growth N' = rN(1 - N/K) (K = 1.2, N0 = 0.05, r = 0.8 per h)
// or drug clearance C' = -kC (k = 0.8 per h, C0 = 10 mg/L).
//  (a) direction field, exact solution, Euler steps drawn as tangent segments, the next
//      step (dashed) and the exact solution through the current Euler point (local curve);
//      Step / Run / Reset walk through the steps one at a time.
//  (b) global error at t = Te (5 h logistic, 6 h clearance) against h, log-log, slope 1,
//      with the unstable range h > 2/k shaded and a marker at the current h.
// Tab "Systems": the SIR epidemic (N = 10 000, beta = 0.3, gamma = 0.1 per day) with the
// error in the peak of I against h (exact peak 3007.96), and the 24-hour oscillator
// x'' = -w^2 x in the (x, x'/w) plane, with the amplitude error after 10 days against h.
// Auto-mounts on <div data-demo="eulerstep">. Optional attributes:
//   data-tab="one|sys"  data-model="logistic|clearance|sir|osc"  data-h="1"
//   data-systems="1|0" (show the Systems tab)  data-field="1|0" (direction field)
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, autoMount } from "../lib/ui";
import {
  Scalar, LOGISTIC, CLEARANCE, euler1, eulerSys, rk4Sys, stepsTo, globalErr, errCurve,
  SIR, sirF, IMAX, sirPeak, OSC, oscF, amp,
} from "./math";

const NAME = "eulerstep";
// Okabe-Ito, as in make-euler-figs.R: exact blue, Euler vermillion, S/I/R blue/red/green
const COL = {
  exact: "#0072B2", euler: "#D55E00", local: "#999999", next: "#E69F00", grn: "#009E73",
  sky: "#56B4E9", field: "#b5b5b5", ref: "#737373", ink: "#222", S: "#0072B2", I: "#D55E00", R: "#009E73",
};
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .es-row { display: flex; flex-wrap: wrap; gap: .35em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .es-group { display: flex; gap: .35em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .es-row .qb-slider { min-width: 300px; flex: 1 1 300px; max-width: 460px; grid-template-columns: auto 1fr 4.6em; }
[data-demo="${NAME}"] .es-plots { display: flex; flex-wrap: wrap; gap: 10px; width: 100%; }
[data-demo="${NAME}"] .es-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .es-read { display: flex; flex-wrap: wrap; gap: .1em 1.3em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .es-read b { font-weight: 600; }
[data-demo="${NAME}"] .es-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .es-tabs .qb-btn { font-size: .95em; }
[data-demo="${NAME}"] .es-box { display: flex; flex-direction: column; gap: .4em; width: 100%; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; gap: .25em; }
.reveal [data-demo="${NAME}"] .es-box { gap: .25em; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 16px; padding: .3em .5em; }
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
const pow10 = (k: number): string => k === 0 ? "1" : k === 1 ? "10" : "10" + String(k).split("").map(c => sup[c]).join("");
function sci(v: number, d = 3): string {
  if (!Number.isFinite(v)) return "–";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e-3 && a < 1e5) return v.toPrecision(d);
  const [m, e] = v.toExponential(d - 1).split("e");
  return `${m}e${Number(e)}`;
}
const fix = (v: number, d: number): string => Number.isFinite(v) ? (Math.abs(v) < 1e5 ? v.toFixed(d) : sci(v, 3)) : "–";
const fmtH = (h: number): string => String(h);
const minus = (s: string): string => s.replace(/^-/, "−").replace(/e-/, "e−");

let mountCount = 0;

// ---------------------------------------------------------------- shared helpers
interface Stepper { n: number; N: number; stop(): void; full(N: number): void; }

/** Step / Run / Reset / All buttons that walk n from 0 to N. */
function stepper(parent: HTMLElement, draw: () => void): Stepper {
  const g = htmlEl("div", "es-group", "", parent);
  let timer = 0;
  const s: Stepper = {
    n: 0, N: 0,
    stop() { if (timer) { window.clearInterval(timer); timer = 0; } runB.textContent = "Run"; },
    full(N: number) { s.stop(); s.N = N; s.n = N; },
  };
  button(g, "Step", () => { s.stop(); s.n = s.n >= s.N ? 1 : s.n + 1; s.n = Math.min(s.n, s.N); draw(); });
  const runB = button(g, "Run", () => {
    if (timer) { s.stop(); return; }
    if (s.n >= s.N) s.n = 0;
    const slow = s.N <= 50, chunk = slow ? 1 : Math.ceil(s.N / 75);
    runB.textContent = "Pause";
    timer = window.setInterval(() => {
      s.n = Math.min(s.N, s.n + chunk);
      draw();
      if (s.n >= s.N) s.stop();
    }, slow ? 220 : 40);
    draw();
  });
  button(g, "Reset", () => { s.stop(); s.n = 0; draw(); });
  button(g, "All steps", () => { s.stop(); s.n = s.N; draw(); });
  return s;
}

/** Points of a polyline, thinned to at most ~2500 vertices; y clamped to stay drawable. */
function poly(xs: number[], ys: number[], fx: (v: number) => number, fy: (v: number) => number, ylo: number, yhi: number): string {
  const n = xs.length, step = Math.max(1, Math.floor(n / 2500));
  const span = yhi - ylo, out: string[] = [];
  for (let i = 0; i < n; i += step) {
    if (!Number.isFinite(ys[i])) continue;
    out.push(`${fx(xs[i]).toFixed(1)},${fy(Math.max(ylo - 5 * span, Math.min(yhi + 5 * span, ys[i]))).toFixed(1)}`);
  }
  if ((n - 1) % step && Number.isFinite(ys[n - 1])) out.push(`${fx(xs[n - 1]).toFixed(1)},${fy(Math.max(ylo - 5 * span, Math.min(yhi + 5 * span, ys[n - 1]))).toFixed(1)}`);
  return out.join(" ");
}

interface Frame { L: number; R: number; T: number; B: number; }

function axesLin(svg: SVGSVGElement, fr: Frame, sx: (v: number) => number, sy: (v: number) => number,
  xt: number[], yt: number[], xlab: string, ylab: string, yfmt: (v: number) => string, H: number): void {
  for (const v of yt) {
    svgEl("line", { x1: fr.L, x2: fr.R, y1: sy(v), y2: sy(v), class: v === 0 ? "axis" : "grid" }, svg);
    svgEl("text", { x: fr.L - 6, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = yfmt(v);
  }
  svgEl("line", { x1: fr.L, x2: fr.R, y1: fr.B, y2: fr.B, class: "axis" }, svg);
  svgEl("line", { x1: fr.L, x2: fr.L, y1: fr.T, y2: fr.B, class: "axis" }, svg);
  for (const v of xt) {
    svgEl("line", { x1: sx(v), x2: sx(v), y1: fr.B, y2: fr.B + 5, class: "axis" }, svg);
    svgEl("text", { x: sx(v), y: fr.B + 22, "text-anchor": "middle", class: "sm" }, svg).textContent = String(v);
  }
  svgEl("text", { x: (fr.L + fr.R) / 2, y: H - 6, "text-anchor": "middle" }, svg).textContent = xlab;
  const cy = (fr.T + fr.B) / 2;
  svgEl("text", { x: 16, y: cy, "text-anchor": "middle", transform: `rotate(-90 16 ${cy})` }, svg).textContent = ylab;
}

function axesLog(svg: SVGSVGElement, fr: Frame, sx: (v: number) => number, sy: (v: number) => number,
  x0: number, x1: number, y0: number, y1: number, xlab: string, ylab: string, H: number): void {
  for (let k = Math.ceil(y0); k <= y1; k++) {
    svgEl("line", { x1: fr.L, x2: fr.R, y1: sy(k), y2: sy(k), class: "grid" }, svg);
    svgEl("text", { x: fr.L - 6, y: sy(k) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = pow10(k);
  }
  svgEl("line", { x1: fr.L, x2: fr.R, y1: fr.B, y2: fr.B, class: "axis" }, svg);
  svgEl("line", { x1: fr.L, x2: fr.L, y1: fr.T, y2: fr.B, class: "axis" }, svg);
  for (let k = Math.ceil(x0); k <= x1; k++) {
    svgEl("line", { x1: sx(k), x2: sx(k), y1: fr.T, y2: fr.B, class: "grid" }, svg);
    svgEl("text", { x: sx(k), y: fr.B + 22, "text-anchor": "middle", class: "sm" }, svg).textContent = k === 0 ? "1" : pow10(k);
  }
  svgEl("text", { x: (fr.L + fr.R) / 2, y: H - 6, "text-anchor": "middle" }, svg).textContent = xlab;
  const cy = (fr.T + fr.B) / 2;
  svgEl("text", { x: 16, y: cy, "text-anchor": "middle", transform: `rotate(-90 16 ${cy})` }, svg).textContent = ylab;
}

function clipGroup(svg: SVGSVGElement, fr: Frame, id: string): SVGGElement {
  const cp = svgEl("clipPath", { id }, svg);
  svgEl("rect", { x: fr.L, y: fr.T, width: fr.R - fr.L, height: fr.B - fr.T }, cp);
  return svgEl("g", { "clip-path": `url(#${id})` }, svg);
}

/** One row of legend items [text, colour, dash, width?] starting at (x, y). */
function legend(svg: SVGSVGElement, x: number, y: number, items: [string, string, string, number?][]): void {
  for (const [txt, c, dash, wd] of items) {
    svgEl("line", { x1: x, x2: x + 24, y1: y - 6, y2: y - 6, stroke: c, "stroke-width": wd ?? 3, "stroke-dasharray": dash, "stroke-opacity": wd && wd > 4 ? 0.35 : 1 }, svg);
    const t = svgEl("text", { x: x + 30, y }, svg);
    t.textContent = txt;
    x += 30 + t.getComputedTextLength() + 18;
    if (!t.getComputedTextLength()) x += txt.length * 9;
  }
}

function readItem(ro: HTMLElement, k: string, v: string, cls = ""): void {
  const s = htmlEl("span", cls, "", ro);
  if (k) htmlEl("b", "", k + " ", s);
  s.append(v);
}

// ---------------------------------------------------------------- mount
function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const uid = ++mountCount;
  const withSys = ds.systems !== "0";
  const startModel = ds.model ?? "logistic";
  const startTab = ds.tab === "sys" || startModel === "sir" || startModel === "osc" ? "sys" : "one";
  const startH = Number(ds.h);

  const oneBox = htmlEl("div", "es-box", "");
  const sysBox = htmlEl("div", "es-box", "");
  const tabRow = htmlEl("div", "es-row es-tabs", "", host);
  let tb: HTMLButtonElement[] = [];
  const show = (t: string) => {
    oneBox.style.display = t === "one" ? "flex" : "none";
    sysBox.style.display = t === "sys" ? "flex" : "none";
    tb.forEach((b, i) => b.classList.toggle("active", (i === 0) === (t === "one")));
  };
  if (withSys) {
    tb = presets(tabRow, [
      { label: "One equation", apply: () => { sysApi.stop(); show("one"); } },
      { label: "Systems: SIR, oscillator", apply: () => { oneApi.stop(); show("sys"); } },
    ]);
  }
  host.append(oneBox, sysBox);
  const oneApi = buildOne(oneBox, withSys ? tabRow : null, uid,
    startModel === "clearance" ? CLEARANCE : LOGISTIC, startTab === "one" ? startH : NaN, ds.field !== "0");
  const sysApi = withSys ? buildSys(sysBox, tabRow, uid, startModel === "osc" ? "osc" : "sir", startTab === "sys" ? startH : NaN)
    : { stop: () => { /* no systems tab */ } };
  show(withSys ? startTab : "one");
}

/** A group for the model buttons: on the tab row (shown only while `box` is visible), or
 *  on a row of its own when there are no tabs. */
function tabGroup(tabRow: HTMLElement | null, box: HTMLElement): HTMLElement {
  const r1 = tabRow ?? htmlEl("div", "es-row", "", box);
  const mg = htmlEl("div", "es-group", "", r1);
  if (!tabRow) htmlEl("span", "", "model:", mg);
  else {
    const obs = new MutationObserver(() => { mg.style.display = box.style.display === "none" ? "none" : "flex"; });
    obs.observe(box, { attributes: true, attributeFilter: ["style"] });
  }
  return mg;
}

// ---------------------------------------------------------------- tab 1: one equation
function buildOne(box: HTMLElement, tabRow: HTMLElement | null, uid: number, m0: Scalar, h0: number, field0: boolean) {
  const pickH = (m: Scalar, h: number) => {
    const i = m.hList.findIndex(v => Math.abs(v - h) < 1e-12);
    return i >= 0 ? i : m.hList.indexOf(1);
  };
  const init = { m: m0, hi: pickH(m0, h0), span: m0.spans[0], field: field0 };
  const st = { ...init };

  // row 1: model (on the tab row if there is one)
  const mg = tabGroup(tabRow, box);
  const models = [LOGISTIC, CLEARANCE];
  const mBtns = presets(mg, models.map(mm => ({ label: mm.label, apply: () => setModel(mm) })));
  // row 2: h and stepping
  const r2 = htmlEl("div", "es-row", "", box);
  const hS = slider(r2, "step h", 0, st.m.hList.length - 1, 1, st.hi, v => { st.hi = v; rerun(); }, v => fmtH(st.m.hList[v] ?? NaN) + " h");
  const sp = stepper(r2, () => draw());
  // row 3: presets, span, field
  const r3 = htmlEl("div", "es-row", "", box);
  const pg = htmlEl("div", "es-group", "", r3);
  const pre = presets(pg, [
    { label: "reset", apply: () => { st.m = init.m; st.span = init.span; st.field = init.field; fieldCb.checked = st.field; resetSlider(); st.hi = init.hi; hS.set(st.hi); rerun(); } },
    { label: "logistic h = 1", apply: () => apply(LOGISTIC, 1, 10) },
    { label: "h = 0.1", apply: () => apply(st.m, 0.1, st.m.spans[0]) },
    { label: "logistic h = 3, 60 h", apply: () => apply(LOGISTIC, 3, 60) },
    { label: "clearance h = 2", apply: () => apply(CLEARANCE, 2, 12) },
    { label: "clearance h = 3", apply: () => apply(CLEARANCE, 3, 12) },
  ]);
  pre.forEach(b => b.addEventListener("click", () => window.setTimeout(() => b.classList.remove("active"), 600)));
  const sg = htmlEl("div", "es-group", "", r3);
  htmlEl("span", "", "time span:", sg);
  const spanBtns = presets(sg, [0, 1].map(i => ({ label: "", apply: () => { st.span = st.m.spans[i]; rerun(); } })));
  const fieldCb = checkbox(r3, "direction field", st.field, v => { st.field = v; draw(); });

  function resetSlider(): void { hS.input.max = String(st.m.hList.length - 1); }
  function setModel(mm: Scalar): void {
    const h = st.m.hList[st.hi];
    st.m = mm; st.span = mm.spans[0];
    resetSlider(); st.hi = pickH(mm, h); hS.set(st.hi); rerun();
  }
  function apply(mm: Scalar, h: number, span: number): void {
    st.m = mm; st.span = span; resetSlider(); st.hi = pickH(mm, h); hS.set(st.hi); rerun();
  }

  // plots
  const plots = htmlEl("div", "es-plots", "", box);
  const H = 370, WA = 640, WB = 460, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Direction field, exact solution and Euler steps" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Global error against step size, log-log" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro = htmlEl("div", "qb-readout es-read", "", box);

  const curves = new Map<string, { h: number; e: number }[]>();
  const curveOf = (mm: Scalar) => {
    let c = curves.get(mm.key);
    if (!c) { c = errCurve(mm, 0.004); curves.set(mm.key, c); }
    return c;
  };
  let ys: number[] = [];

  function rerun(): void {
    const h = st.m.hList[st.hi];
    const N = Math.floor(st.span / h + 1e-9);
    ys = euler1(st.m.f, st.m.y0, h, N);
    sp.full(N);
    draw();
  }

  function draw(): void {
    const m = st.m, h = m.hList[st.hi], n = sp.n, N = sp.N;
    mBtns.forEach((b, i) => b.classList.toggle("active", models[i] === m));
    spanBtns.forEach((b, i) => { b.textContent = `${m.spans[i]} h`; b.classList.toggle("active", m.spans[i] === st.span); });
    const [ylo, yhi] = m.ylim;

    // ---------- (a) ----------
    svgA.replaceChildren();
    const fr: Frame = { L: 66, R: WA - 14, T: 60, B: H - 50 };
    const sx = scale(0, st.span, fr.L, fr.R), sy = scale(ylo, yhi, fr.B, fr.T);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgA).textContent =
      `(a) ${m.key === "logistic" ? "N′ = 0.8 N (1 − N/1.2)" : "C′ = −0.8 C"},  h = ${fmtH(h)} h`;
    const tstep = st.span <= 12 ? 2 : st.span <= 24 ? 4 : 10;
    const xt: number[] = []; for (let t = 0; t <= st.span + 1e-9; t += tstep) xt.push(t);
    const yt: number[] = []; const ystep = m.key === "logistic" ? 0.4 : 5;
    for (let v = Math.ceil(ylo / ystep) * ystep; v <= yhi + 1e-9; v += ystep) yt.push(+v.toFixed(6));
    axesLin(svgA, fr, sx, sy, xt, yt, "time t (h)", m.yName, v => m.key === "logistic" ? v.toFixed(1) : String(v), H);
    const g = clipGroup(svgA, fr, `es-a-${uid}`);
    // direction field: short segments with slope f(y), in screen units
    if (st.field) {
      const nx = 18, ny = 11, len = 9;
      const kx = (fr.R - fr.L) / st.span, ky = (fr.B - fr.T) / (yhi - ylo);
      for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
        const t = (i + 0.5) * st.span / nx, y = ylo + (j + 0.5) * (yhi - ylo) / ny;
        const dx = kx, dy = -m.f(y) * ky, s = len / Math.hypot(dx, dy);
        svgEl("line", { x1: sx(t) - dx * s, y1: sy(y) - dy * s, x2: sx(t) + dx * s, y2: sy(y) + dy * s, stroke: COL.field, "stroke-width": 1.6 }, g);
      }
    }
    if (m.key === "logistic") svgEl("line", { x1: fr.L, x2: fr.R, y1: sy(1.2), y2: sy(1.2), stroke: COL.grn, "stroke-width": 1.5, "stroke-dasharray": "6 5" }, g);
    // error time
    svgEl("line", { x1: sx(m.Te), x2: sx(m.Te), y1: fr.T, y2: fr.B, stroke: COL.ref, "stroke-width": 1.2, "stroke-dasharray": "3 4" }, g);
    if (m.Te < st.span) svgEl("text", { x: sx(m.Te) + 4, y: fr.T + 18, class: "sm" }, g).textContent = `t = ${m.Te}`;
    // exact
    const tt = Array.from({ length: 301 }, (_, i) => i * st.span / 300);
    svgEl("polyline", { points: poly(tt, tt.map(t => m.exactFrom(t, 0, m.y0)), sx, sy, ylo, yhi), fill: "none", stroke: COL.exact, "stroke-width": 3 }, g);
    // Euler steps so far
    const ts = ys.map((_, i) => i * h);
    const tn = n * h, yn = ys[n];
    if (n > 0) svgEl("polyline", { points: poly(ts.slice(0, n + 1), ys.slice(0, n + 1), sx, sy, ylo, yhi), fill: "none", stroke: COL.euler, "stroke-width": 2.6 }, g);
    if (N <= 80) for (let i = 0; i <= n; i++) {
      if (ys[i] >= ylo - 1 && ys[i] <= yhi + 1) svgEl("circle", { cx: sx(ts[i]), cy: sy(ys[i]), r: 4.5, fill: COL.euler }, g);
    }
    // exact solution through the current point, and the next step
    if (n > 0 && Number.isFinite(yn)) {
      const tl = Array.from({ length: 121 }, (_, i) => tn + i * (st.span - tn) / 120);
      svgEl("polyline", { points: poly(tl, tl.map(t => m.exactFrom(t, tn, yn)), sx, sy, ylo, yhi), fill: "none", stroke: COL.local, "stroke-width": 2, "stroke-dasharray": "7 5" }, g);
    }
    if (n < N) {
      const y1 = yn + h * m.f(yn);
      const cl = (v: number) => Math.max(ylo - 5 * (yhi - ylo), Math.min(yhi + 5 * (yhi - ylo), v));
      svgEl("line", { x1: sx(tn), y1: sy(cl(yn)), x2: sx(tn + h), y2: sy(cl(y1)), stroke: COL.next, "stroke-width": 3, "stroke-dasharray": "8 5" }, g);
      if (y1 >= ylo && y1 <= yhi) svgEl("circle", { cx: sx(tn + h), cy: sy(y1), r: 5, fill: "#fff", stroke: COL.next, "stroke-width": 2.5 }, g);
    }
    if (yn >= ylo && yn <= yhi) svgEl("circle", { cx: sx(tn), cy: sy(yn), r: 7, fill: COL.euler, stroke: COL.ink, "stroke-width": 2 }, g);
    else {
      const up = yn > yhi;
      const t = svgEl("text", { x: fr.R - 8, y: up ? fr.T + 22 : fr.B - 10, "text-anchor": "end", "font-weight": 600 }, svgA);
      t.style.fill = COL.euler; t.textContent = `${up ? "↑" : "↓"} off scale: y = ${minus(sci(yn, 3))}`;
    }
    // legend
    legend(svgA, 10, 48, [
      ["exact", COL.exact, ""], ["Euler", COL.euler, ""],
      ["exact through current point", COL.local, "7 5"], ["next step", COL.next, "8 5"],
    ]);

    // ---------- (b) ----------
    svgB.replaceChildren();
    const fb: Frame = { L: 70, R: WB - 18, T: 60, B: H - 50 };
    const X0 = -2.4, X1 = Math.log10(m.Te) + 0.12, Y0 = -5, Y1 = 2;
    const bx = scale(X0, X1, fb.L, fb.R), by = scale(Y0, Y1, fb.B, fb.T);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgB).textContent = `(b) error at t = ${m.Te} h against h`;
    axesLog(svgB, fb, bx, by, X0, X1, Y0, Y1, "step h (log)", "|global error|", H);
    const gb = clipGroup(svgB, fb, `es-b-${uid}`);
    const hu = Math.log10(2 / m.rate);
    svgEl("rect", { x: bx(hu), y: fb.T, width: fb.R - bx(hu), height: fb.B - fb.T, fill: COL.euler, "fill-opacity": 0.1 }, gb);
    const ut = svgEl("text", { x: bx(hu) + 20, y: fb.B - 8, transform: `rotate(-90 ${bx(hu) + 20} ${fb.B - 8})`, class: "sm" }, gb);
    ut.textContent = m.key === "logistic" ? "h > 2/r: unstable near K" : "h > 2/k: unstable";
    const cy = (e: number) => by(Math.max(Y0 - 1, Math.min(Y1 + 1, Math.log10(Math.max(Math.abs(e), 1e-30)))));
    const cv = curveOf(m);
    const last = cv[cv.length - 1];
    const C = Math.abs(last.e) / last.h;
    svgEl("line", { x1: bx(X0), y1: by(Math.log10(C) + X0), x2: bx(X1), y2: by(Math.log10(C) + X1), stroke: COL.ref, "stroke-width": 2, "stroke-dasharray": "9 6" }, gb);
    const lx = -1.6;
    const lt = svgEl("text", { x: bx(lx) + 10, y: by(Math.log10(C) + lx) + 26, class: "sm" }, gb);
    lt.textContent = `slope 1:  ${C.toPrecision(2)}·h`;
    svgEl("polyline", { points: cv.map(p => `${bx(Math.log10(p.h)).toFixed(1)},${cy(p.e).toFixed(1)}`).join(" "), fill: "none", stroke: COL.euler, "stroke-width": 2.2 }, gb);
    cv.forEach(p => svgEl("circle", { cx: bx(Math.log10(p.h)), cy: cy(p.e), r: 3, fill: COL.euler }, gb));
    const e = globalErr(m, h), eh2 = globalErr(m, h / 2);
    if (Number.isFinite(e)) {
      svgEl("circle", { cx: bx(Math.log10(h)), cy: cy(e), r: 9, fill: COL.next, stroke: COL.ink, "stroke-width": 2.5 }, svgB);
    } else {
      svgEl("line", { x1: bx(Math.log10(h)), x2: bx(Math.log10(h)), y1: fb.T, y2: fb.B, stroke: COL.ink, "stroke-width": 1, "stroke-dasharray": "2 4" }, svgB);
      const nt = svgEl("text", { x: bx(Math.log10(h)) - 6, y: fb.T + 22, "text-anchor": "end", class: "sm" }, svgB);
      nt.textContent = `t = ${m.Te} is not a step point`;
    }

    // ---------- readout ----------
    ro.replaceChildren();
    const fac = 1 - m.rate * h;
    readItem(ro, "step", `${n} of ${N}`);
    readItem(ro, "t =", `${+tn.toFixed(4)} h`);
    readItem(ro, m.key === "logistic" ? "N =" : "C =", minus(fix(yn, 5)));
    readItem(ro, "slope f =", minus(fix(m.f(yn), 5)));
    readItem(ro, "exact here", fix(m.exactFrom(tn, 0, m.y0), 5));
    readItem(ro, `factor 1 − ${m.rateName}h =`, minus(fac.toFixed(2)) + (Math.abs(fac) > 1 ? "  (|·| > 1: unstable)" : ""), Math.abs(fac) > 1 ? "es-warn" : "");
    if (Number.isFinite(e)) {
      readItem(ro, `at t = ${m.Te}: Euler`, minus(fix(globalErr(m, h) + m.exactFrom(m.Te, 0, m.y0), 6)));
      readItem(ro, "exact", fix(m.exactFrom(m.Te, 0, m.y0), 6));
      readItem(ro, "error", minus(sci(e, 4)));
      if (Number.isFinite(eh2) && eh2 !== 0) readItem(ro, "error(h)/error(h/2) =", (e / eh2).toFixed(2));
    } else readItem(ro, `t = ${m.Te} h`, `is not a multiple of h = ${fmtH(h)}: no error at t = ${m.Te}`);
  }

  // start
  rerun();
  return { stop: () => sp.stop() };
}

// ---------------------------------------------------------------- tab 2: systems
function buildSys(box: HTMLElement, tabRow: HTMLElement, uid: number, model0: "sir" | "osc", h0: number) {
  const LIST = {
    sir: [0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 4, 5, 8],
    osc: [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2],
  };
  const SPANS = { sir: [160], osc: [48, 240] };
  const pickH = (k: "sir" | "osc", h: number) => {
    const i = LIST[k].findIndex(v => Math.abs(v - h) < 1e-12);
    return i >= 0 ? i : LIST[k].indexOf(1);
  };
  const st = { k: model0, hi: pickH(model0, h0), span: SPANS[model0][0] };

  const mg = tabGroup(tabRow, box);
  const keys: ("sir" | "osc")[] = ["sir", "osc"];
  const mBtns = presets(mg, [
    { label: "SIR epidemic", apply: () => setModel("sir") },
    { label: "24-hour oscillator", apply: () => setModel("osc") },
  ]);
  const r2 = htmlEl("div", "es-row", "", box);
  const unit = () => (st.k === "sir" ? " d" : " h");
  const hS = slider(r2, "step h", 0, LIST[st.k].length - 1, 1, st.hi, v => { st.hi = v; rerun(); }, v => String(LIST[st.k][v] ?? "") + unit());
  const sp = stepper(r2, () => draw());
  const r3 = htmlEl("div", "es-row", "", box);
  const QUICK = { sir: [4, 1, 0.1, 0.01], osc: [1, 0.1, 0.01] };
  const qg = htmlEl("div", "es-group", "", r3);
  htmlEl("span", "", "h =", qg);
  const qBtns = presets(qg, [0, 1, 2, 3].map(i => ({ label: "", apply: () => { const h = QUICK[st.k][i]; if (h !== undefined) { st.hi = pickH(st.k, h); hS.set(st.hi); rerun(); } } })));
  const sg = htmlEl("div", "es-group", "", r3);
  const spanLab = htmlEl("span", "", "time span:", sg);
  const spanBtns = presets(sg, [0, 1].map(i => ({ label: "", apply: () => { st.span = SPANS[st.k][i] ?? st.span; rerun(); } })));

  function setModel(k: "sir" | "osc"): void {
    st.k = k; st.span = SPANS[k][0];
    hS.input.max = String(LIST[k].length - 1);
    st.hi = LIST[k].indexOf(1); hS.set(st.hi); rerun();
  }

  const plots = htmlEl("div", "es-plots", "", box);
  const H = 370, WA = 640, WB = 460, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Euler solution of a system" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Error against step size, log-log" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro = htmlEl("div", "qb-readout es-read", "", box);

  // reference curves and error-vs-h curves, computed once
  let sirRef: number[][] | null = null, oscCurve: { h: number; e: number }[] | null = null, sirCurve: { h: number; e: number }[] | null = null;
  const getSirRef = () => sirRef ?? (sirRef = rk4Sys(sirF, SIR.y0, 0.05, Math.round(SIR.T / 0.05)));
  const hGrid = (a: number, b: number) => { const o: number[] = []; for (let l = a; l <= b + 1e-9; l += 0.1) o.push(10 ** l); return o; };
  const getSirCurve = () => sirCurve ?? (sirCurve = hGrid(-2, 0.9).map(h => {
    const Y = eulerSys(sirF, SIR.y0, h, Math.floor(SIR.T / h + 1e-9));
    return { h, e: sirPeak(Y, h).peak - IMAX };
  }));
  const getOscCurve = () => oscCurve ?? (oscCurve = hGrid(-2.3, 0.3).map(h => {
    const n = Math.round(OSC.T / h);
    const g = Math.sqrt(1 + OSC.w * OSC.w * h * h);
    return { h, e: g ** n - 1 };                 // exact for Euler: amplitude x g per step
  }));

  let Y: number[][] = [];
  function rerun(): void {
    const h = LIST[st.k][st.hi];
    const N = Math.floor(st.span / h + 1e-9);
    Y = eulerSys(st.k === "sir" ? sirF : oscF, st.k === "sir" ? SIR.y0 : OSC.y0, h, N);
    sp.full(N);
    draw();
  }

  function draw(): void {
    const h = LIST[st.k][st.hi], n = sp.n, N = sp.N;
    mBtns.forEach((b, i) => b.classList.toggle("active", keys[i] === st.k));
    const q = QUICK[st.k];
    qBtns.forEach((b, i) => {
      b.style.display = i < q.length ? "" : "none";
      if (i < q.length) { b.textContent = `${q[i]}${unit()}`; b.classList.toggle("active", Math.abs(q[i] - h) < 1e-12); }
    });
    const spans = SPANS[st.k];
    spanLab.style.display = spans.length > 1 ? "" : "none";
    spanBtns.forEach((b, i) => {
      b.style.display = i < spans.length && spans.length > 1 ? "" : "none";
      if (i < spans.length) { b.textContent = `${spans[i] / 24} days`; b.classList.toggle("active", spans[i] === st.span); }
    });
    if (st.k === "sir") drawSir(h, n, N); else drawOsc(h, n, N);
  }

  function drawSir(h: number, n: number, N: number): void {
    svgA.replaceChildren();
    const fr: Frame = { L: 76, R: WA - 22, T: 60, B: H - 50 };
    const sx = scale(0, SIR.T, fr.L, fr.R), sy = scale(0, 10000, fr.B, fr.T);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgA).textContent = `(a) SIR, N = 10 000, β = 0.3, γ = 0.1,  h = ${h} d`;
    axesLin(svgA, fr, sx, sy, [0, 40, 80, 120, 160], [0, 2500, 5000, 7500, 10000], "time t (days)", "people", v => String(v), H);
    const g = clipGroup(svgA, fr, `es-s-${uid}`);
    const ref = getSirRef();
    const tr = ref.map((_, i) => i * 0.05);
    const cols = [COL.S, COL.I, COL.R];
    cols.forEach((c, j) => svgEl("polyline", { points: poly(tr, ref.map(y => y[j]), sx, sy, 0, 10000), fill: "none", stroke: c, "stroke-width": 5, "stroke-opacity": 0.3 }, g));
    svgEl("line", { x1: fr.L, x2: fr.R, y1: sy(IMAX), y2: sy(IMAX), stroke: COL.I, "stroke-width": 1.5, "stroke-dasharray": "6 5" }, g);
    const te = Y.slice(0, n + 1).map((_, i) => i * h);
    cols.forEach((c, j) => {
      if (n > 0) svgEl("polyline", { points: poly(te, Y.slice(0, n + 1).map(y => y[j]), sx, sy, 0, 10000), fill: "none", stroke: c, "stroke-width": 2.4 }, g);
      if (N <= 200) for (let i = 0; i <= n; i++) svgEl("circle", { cx: sx(te[i]), cy: sy(Y[i][j]), r: 3.5, fill: c }, g);
    });
    const pk = sirPeak(Y.slice(0, n + 1), h);
    if (n > 0) svgEl("circle", { cx: sx(pk.day), cy: sy(pk.peak), r: 7, fill: "none", stroke: COL.ink, "stroke-width": 2.5 }, g);
    legend(svgA, 10, 48, [["S", COL.S, ""], ["I", COL.I, ""], ["R", COL.R, ""], ["exact peak of I", COL.I, "6 5"], ["reference", COL.ref, "", 7]]);

    // (b) peak error vs h
    svgB.replaceChildren();
    const fb: Frame = { L: 70, R: WB - 18, T: 60, B: H - 50 };
    const X0 = -2.15, X1 = 1.0, Y0 = -1, Y1 = 4;
    const bx = scale(X0, X1, fb.L, fb.R), by = scale(Y0, Y1, fb.B, fb.T);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgB).textContent = "(b) error in the peak of I against h";
    axesLog(svgB, fb, bx, by, X0, X1, Y0, Y1, "step h (days, log)", "|peak − 3008|", H);
    const gb = clipGroup(svgB, fb, `es-sb-${uid}`);
    const cv = getSirCurve();
    const cy = (e: number) => by(Math.max(Y0 - 1, Math.min(Y1 + 1, Math.log10(Math.max(Math.abs(e), 1e-30)))));
    const C = Math.abs(cv[0].e) / cv[0].h;
    svgEl("line", { x1: bx(X0), y1: by(Math.log10(C) + X0), x2: bx(X1), y2: by(Math.log10(C) + X1), stroke: COL.ref, "stroke-width": 2, "stroke-dasharray": "9 6" }, gb);
    svgEl("text", { x: bx(-1.2), y: by(Math.log10(C) - 1.2) + 30, class: "sm" }, gb).textContent = `slope 1:  ${C.toPrecision(3)}·h`;
    svgEl("polyline", { points: cv.map(p => `${bx(Math.log10(p.h)).toFixed(1)},${cy(p.e).toFixed(1)}`).join(" "), fill: "none", stroke: COL.I, "stroke-width": 2.2 }, gb);
    cv.forEach(p => svgEl("circle", { cx: bx(Math.log10(p.h)), cy: cy(p.e), r: 3, fill: COL.I }, gb));
    const full = sirPeak(Y, h);
    svgEl("circle", { cx: bx(Math.log10(h)), cy: cy(full.peak - IMAX), r: 9, fill: COL.next, stroke: COL.ink, "stroke-width": 2.5 }, svgB);

    ro.replaceChildren();
    let dev = 0;
    for (const y of Y) dev = Math.max(dev, Math.abs(y[0] + y[1] + y[2] - SIR.N));
    const y = Y[n];
    readItem(ro, "step", `${n} of ${N}`);
    readItem(ro, "t =", `${+(n * h).toFixed(4)} d`);
    readItem(ro, "S, I, R =", `${y[0].toFixed(1)}, ${y[1].toFixed(1)}, ${y[2].toFixed(1)}`);
    readItem(ro, "Euler peak", `${full.peak.toFixed(3)} on day ${+full.day.toFixed(3)}`);
    readItem(ro, "exact", IMAX.toFixed(3));
    readItem(ro, "error", (full.peak - IMAX >= 0 ? "+" : "−") + Math.abs(full.peak - IMAX).toFixed(3));
    readItem(ro, "max |S + I + R − 10 000| =", sci(dev, 2));
  }

  function drawOsc(h: number, n: number, N: number): void {
    const w = OSC.w;
    svgA.replaceChildren();
    // phase plane, left square
    const fp: Frame = { L: 66, R: 66 + 260, T: 60, B: 60 + 260 };
    const L = 2.6;
    const px = scale(-L, L, fp.L, fp.R), py = scale(-L, L, fp.B, fp.T);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgA).textContent = `(a) x″ = −ω²x,  ω = 2π/24 per h,  h = ${h} h`;
    for (const v of [-2, -1, 0, 1, 2]) {
      svgEl("line", { x1: fp.L, x2: fp.R, y1: py(v), y2: py(v), class: v === 0 ? "axis" : "grid" }, svgA);
      svgEl("line", { x1: px(v), x2: px(v), y1: fp.T, y2: fp.B, class: v === 0 ? "axis" : "grid" }, svgA);
      svgEl("text", { x: fp.L - 6, y: py(v) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = String(v);
      svgEl("text", { x: px(v), y: fp.B + 20, "text-anchor": "middle", class: "sm" }, svgA).textContent = String(v);
    }
    svgEl("rect", { x: fp.L, y: fp.T, width: fp.R - fp.L, height: fp.B - fp.T, fill: "none", class: "axis" }, svgA);
    svgEl("text", { x: (fp.L + fp.R) / 2, y: H - 6, "text-anchor": "middle" }, svgA).textContent = "x";
    const cyp = (fp.T + fp.B) / 2;
    svgEl("text", { x: 16, y: cyp, "text-anchor": "middle", transform: `rotate(-90 16 ${cyp})` }, svgA).textContent = "x′/ω";
    const g = clipGroup(svgA, fp, `es-p-${uid}`);
    const th = Array.from({ length: 181 }, (_, i) => i * 2 * Math.PI / 180);
    svgEl("polyline", { points: th.map(a => `${px(Math.cos(a)).toFixed(1)},${py(Math.sin(a)).toFixed(1)}`).join(" "), fill: "none", stroke: COL.exact, "stroke-width": 3 }, g);
    const Ys = Y.slice(0, n + 1);
    if (n > 0) svgEl("polyline", { points: poly(Ys.map(y => y[0]), Ys.map(y => y[1] / w), px, py, -L, L), fill: "none", stroke: COL.euler, "stroke-width": 2.2 }, g);
    const yn = Y[n];
    const inside = Math.abs(yn[0]) <= L && Math.abs(yn[1] / w) <= L;
    if (inside) svgEl("circle", { cx: px(yn[0]), cy: py(yn[1] / w), r: 6, fill: COL.euler, stroke: COL.ink, "stroke-width": 2 }, g);
    svgEl("circle", { cx: px(1), cy: py(0), r: 5, fill: COL.ink }, g);
    legend(svgA, 10, 48, [["exact", COL.exact, ""], ["Euler", COL.euler, ""]]);
    // x(t) on the right
    const ft: Frame = { L: fp.R + 62, R: WA - 12, T: fp.T, B: fp.B };
    const tx = scale(0, st.span, ft.L, ft.R), ty = scale(-L, L, ft.B, ft.T);
    for (const v of [-2, -1, 0, 1, 2]) {
      svgEl("line", { x1: ft.L, x2: ft.R, y1: ty(v), y2: ty(v), class: v === 0 ? "axis" : "grid" }, svgA);
      svgEl("text", { x: ft.L - 6, y: ty(v) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = String(v);
    }
    const dstep = st.span <= 48 ? 1 : 2;
    for (let d = 0; d <= st.span / 24 + 1e-9; d += dstep) {
      svgEl("line", { x1: tx(d * 24), x2: tx(d * 24), y1: ft.B, y2: ft.B + 5, class: "axis" }, svgA);
      svgEl("text", { x: tx(d * 24), y: ft.B + 20, "text-anchor": "middle", class: "sm" }, svgA).textContent = String(d);
    }
    svgEl("rect", { x: ft.L, y: ft.T, width: ft.R - ft.L, height: ft.B - ft.T, fill: "none", class: "axis" }, svgA);
    svgEl("text", { x: (ft.L + ft.R) / 2, y: H - 6, "text-anchor": "middle" }, svgA).textContent = "time (days)";
    const cyt = (ft.T + ft.B) / 2;
    svgEl("text", { x: ft.L - 34, y: cyt, "text-anchor": "middle", transform: `rotate(-90 ${ft.L - 34} ${cyt})` }, svgA).textContent = "x(t)";
    const gt = clipGroup(svgA, ft, `es-t-${uid}`);
    const tt = Array.from({ length: 601 }, (_, i) => i * st.span / 600);
    svgEl("polyline", { points: poly(tt, tt.map(t => Math.cos(w * t)), tx, ty, -L, L), fill: "none", stroke: COL.exact, "stroke-width": 2.5 }, gt);
    if (n > 0) svgEl("polyline", { points: poly(Ys.map((_, i) => i * h), Ys.map(y => y[0]), tx, ty, -L, L), fill: "none", stroke: COL.euler, "stroke-width": 2 }, gt);

    // (b) amplitude error after 10 days
    svgB.replaceChildren();
    const fb: Frame = { L: 70, R: WB - 18, T: 60, B: H - 50 };
    const X0 = -2.4, X1 = 0.4, Y0 = -3, Y1 = 4;
    const bx = scale(X0, X1, fb.L, fb.R), by = scale(Y0, Y1, fb.B, fb.T);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgB).textContent = "(b) amplitude − 1 after 10 days";
    axesLog(svgB, fb, bx, by, X0, X1, Y0, Y1, "step h (hours, log)", "amplitude error", H);
    const gb = clipGroup(svgB, fb, `es-ob-${uid}`);
    const cv = getOscCurve();
    const cy = (e: number) => by(Math.max(Y0 - 1, Math.min(Y1 + 1, Math.log10(Math.max(Math.abs(e), 1e-30)))));
    const C = cv[0].e / cv[0].h;
    svgEl("line", { x1: bx(X0), y1: by(Math.log10(C) + X0), x2: bx(X1), y2: by(Math.log10(C) + X1), stroke: COL.ref, "stroke-width": 2, "stroke-dasharray": "9 6" }, gb);
    svgEl("text", { x: bx(-1.4) + 10, y: by(Math.log10(C) - 1.4) + 28, class: "sm" }, gb).textContent = `slope 1:  ${C.toPrecision(3)}·h`;
    svgEl("polyline", { points: cv.map(p => `${bx(Math.log10(p.h)).toFixed(1)},${cy(p.e).toFixed(1)}`).join(" "), fill: "none", stroke: COL.euler, "stroke-width": 2.2 }, gb);
    cv.forEach(p => svgEl("circle", { cx: bx(Math.log10(p.h)), cy: cy(p.e), r: 3, fill: COL.euler }, gb));
    const g1 = Math.sqrt(1 + w * w * h * h);
    const a10 = g1 ** Math.round(OSC.T / h);
    svgEl("circle", { cx: bx(Math.log10(h)), cy: cy(a10 - 1), r: 9, fill: COL.next, stroke: COL.ink, "stroke-width": 2.5 }, svgB);

    ro.replaceChildren();
    readItem(ro, "step", `${n} of ${N}`);
    readItem(ro, "t =", `${+(n * h).toFixed(3)} h`);
    readItem(ro, "amplitude now", fix(amp(yn), 4));
    if (!inside) readItem(ro, "", "(off the plot)", "es-warn");
    readItem(ro, "growth per step √(1 + ω²h²) =", g1.toFixed(6));
    readItem(ro, "per day", (g1 ** (24 / h)).toFixed(4));
    readItem(ro, "amplitude after 10 days", a10 < 1e5 ? a10.toFixed(4) : sci(a10, 4));
  }

  rerun();
  return { stop: () => sp.stop() };
}

autoMount(NAME, mount);

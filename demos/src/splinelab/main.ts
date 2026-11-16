// Spline lab for the splines unit. Four data sets of the unit: the immediate-early gene
// time course (8 uneven points), 11 exact OD values of the logistic growth curve, the 41
// noisy plate-reader reads, and one gene in 40 cells along pseudotime.
//  (a) the data (drag points up and down, or select one and nudge it with − / +) and the
//      interpolants: piecewise linear, natural cubic and not-a-knot cubic splines (both by
//      a tridiagonal Thomas solve), a monotone cubic (PCHIP), the polynomial of degree
//      n − 1 through all points, and the smoothing spline (Reinsch, pentadiagonal solve)
//      with λ from a slider or chosen by GCV;
//  (b) either the slope S′(t) (the growth rate for OD data, against the true rate) or the
//      change |ΔS(t)| on a log axis since the preset data, which shows that moving one
//      point changes a spline only nearby but the polynomial everywhere.
// Click either panel (not on a point) to move the probe time t*; the readout table gives
// each curve's value, slope and change at t*, its minimum and its largest slope.
// Auto-mounts on <div data-demo="splinelab">. Optional attributes:
//   data-set="gene|od|odnoisy|pseudo"  data-lower="slope|change"
//   data-show="linear natural notaknot monotone poly smooth"  data-probe="3"
import { svgEl, htmlEl, slider, presets, button, checkbox, scale, svgPoint, autoMount } from "../lib/ui";
import { Curve, linear, fromM, naturalM, notAKnotM, pchip, polynomial, smoothing, SmoothFit } from "./math";

const NAME = "splinelab";
// Okabe-Ito, as in make-splines-figs.R (linear green, natural spline blue, polynomial purple)
const COL = {
  linear: "#009E73", natural: "#0072B2", notaknot: "#D55E00", monotone: "#E69F00",
  poly: "#CC79A7", smooth: "#000000", pt: "#222222", sel: "#D55E00", ref: "#737373", neg: "#F5DCDC",
};
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .sl-row { display: flex; flex-wrap: wrap; gap: .35em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .sl-group { display: flex; gap: .35em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .sl-row .qb-slider { min-width: 280px; flex: 1 1 300px; max-width: 420px; grid-template-columns: auto 1fr 4.6em; }
[data-demo="${NAME}"] .qb-check { font-weight: 600; }
[data-demo="${NAME}"] .sl-plots { display: flex; flex-wrap: wrap; gap: 10px; width: 100%; }
[data-demo="${NAME}"] .sl-plots svg { min-width: 0; height: auto; touch-action: none; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .qb-readout { width: 100%; box-sizing: border-box; white-space: normal; }
[data-demo="${NAME}"] table.sl-tab { border-collapse: collapse; width: 100%; }
[data-demo="${NAME}"] table.sl-tab td, [data-demo="${NAME}"] table.sl-tab th { padding: .05em .5em; text-align: right; border: none; }
[data-demo="${NAME}"] table.sl-tab th { font-weight: 600; color: #555; font-family: system-ui, sans-serif; }
[data-demo="${NAME}"] table.sl-tab td:first-child, [data-demo="${NAME}"] table.sl-tab th:first-child { text-align: left; font-weight: 600; font-family: system-ui, sans-serif; }
[data-demo="${NAME}"] .sl-info { display: flex; flex-wrap: wrap; gap: .1em 1.4em; margin-bottom: .2em; }
[data-demo="${NAME}"] .sl-warn { color: #b00020; }
@media (max-width: 760px) { [data-demo="${NAME}"] .sl-plots svg { width: 100% !important; } }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 18px; gap: .25em; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 16px; padding: .2em .5em; }
.reveal [data-demo="${NAME}"] .qb-btn { font-size: .85em; padding: .15em .55em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- the unit's data ----
const K = 1.2, N0 = 0.05, R = 0.8;
const od = (t: number): number => K / (1 + (K / N0 - 1) * Math.exp(-R * t));
const rate = (t: number): number => R * od(t) * (1 - od(t) / K);
const mPt = (t: number): number => 1 + 2 * Math.exp(-(((t - 0.6) / 0.15) ** 2));
const dPt = (t: number): number => -4 * (t - 0.6) / 0.15 ** 2 * Math.exp(-(((t - 0.6) / 0.15) ** 2));

type Method = "linear" | "natural" | "notaknot" | "monotone" | "poly" | "smooth";
const METHODS: Method[] = ["linear", "natural", "notaknot", "monotone", "poly", "smooth"];
const LABEL: Record<Method, string> = {
  linear: "piecewise linear", natural: "natural cubic", notaknot: "not-a-knot",
  monotone: "monotone (PCHIP)", poly: "polynomial", smooth: "smoothing spline",
};

interface DataSet {
  key: string; label: string; x: number[]; y: number[];
  xlim: [number, number]; xticks: number[]; ylim: [number, number]; yticks: number[];
  dlim: [number, number]; dticks: number[];
  xlab: string; ylab: string; dlab: string;
  step: number; sel: number; probe: number; show: Method[];
  truth?: (t: number) => number; dtruth?: (t: number) => number; truthLabel?: string;
}
const tOD = Array.from({ length: 41 }, (_, i) => i * 0.25);
const SETS: DataSet[] = [
  {
    key: "gene", label: "gene time course",
    x: [0, 0.5, 1, 2, 4, 8, 12, 24], y: [1, 40, 32, 6, 1, 0.5, 0.5, 0.5],
    xlim: [-0.5, 24.5], xticks: [0, 4, 8, 12, 16, 20, 24], ylim: [-12, 50], yticks: [-10, 0, 10, 20, 30, 40, 50],
    dlim: [-80, 100], dticks: [-80, -40, 0, 40, 80],
    xlab: "time after stimulation (h)", ylab: "expression (TPM)", dlab: "slope (TPM/h)",
    step: 5, sel: 1, probe: 3, show: ["linear", "natural", "poly"],
  },
  {
    key: "od", label: "OD, 11 exact values",
    x: Array.from({ length: 11 }, (_, i) => i), y: Array.from({ length: 11 }, (_, i) => od(i)),
    xlim: [-0.2, 10.2], xticks: [0, 2, 4, 6, 8, 10], ylim: [-0.1, 1.4], yticks: [0, 0.4, 0.8, 1.2],
    dlim: [-0.1, 0.35], dticks: [-0.1, 0, 0.1, 0.2, 0.3],
    xlab: "time (h)", ylab: "OD600", dlab: "growth rate (OD/h)",
    step: 0.05, sel: 5, probe: 0.5, show: ["natural", "notaknot"],
    truth: od, dtruth: rate, truthLabel: "true logistic",
  },
  {
    key: "odnoisy", label: "OD, 41 noisy reads",
    x: tOD, y: [0.040, 0.058, 0.076, 0.077, 0.108, 0.127, 0.152, 0.191, 0.200, 0.263, 0.284, 0.327, 0.382, 0.446,
      0.502, 0.556, 0.610, 0.672, 0.749, 0.794, 0.839, 0.883, 0.934, 0.958, 1.004, 1.032, 1.076, 1.097, 1.105, 1.111,
      1.144, 1.155, 1.163, 1.171, 1.167, 1.182, 1.193, 1.184, 1.177, 1.197, 1.199],
    xlim: [-0.2, 10.2], xticks: [0, 2, 4, 6, 8, 10], ylim: [-0.1, 1.4], yticks: [0, 0.4, 0.8, 1.2],
    dlim: [-0.15, 0.4], dticks: [-0.1, 0, 0.1, 0.2, 0.3, 0.4],
    xlab: "time (h)", ylab: "OD600", dlab: "growth rate (OD/h)",
    step: 0.05, sel: 20, probe: 4, show: ["natural", "smooth"],
    truth: od, dtruth: rate, truthLabel: "true logistic",
  },
  {
    key: "pseudo", label: "pseudotime, 40 cells",
    x: [0.013, 0.062, 0.108, 0.126, 0.177, 0.186, 0.202, 0.206, 0.212, 0.266, 0.267, 0.340, 0.372, 0.380, 0.382,
      0.384, 0.386, 0.411, 0.482, 0.494, 0.498, 0.573, 0.600, 0.629, 0.652, 0.661, 0.668, 0.687, 0.718, 0.724,
      0.770, 0.777, 0.794, 0.827, 0.870, 0.898, 0.908, 0.935, 0.945, 0.992],
    y: [1.32, 1.27, 1.03, 0.30, 1.22, 0.98, 0.95, 0.49, 0.84, 1.16, 1.49, 1.06, 1.33, 1.21, 0.76, 1.11, 1.12, 1.39,
      2.46, 2.48, 2.20, 2.85, 3.24, 3.12, 2.53, 2.45, 2.76, 2.70, 2.04, 2.32, 1.69, 1.28, 1.49, 0.81, 1.58, 1.73,
      0.90, 0.65, 1.21, 0.95],
    xlim: [0, 1], xticks: [0, 0.2, 0.4, 0.6, 0.8, 1], ylim: [-0.5, 4.5], yticks: [0, 1, 2, 3, 4],
    dlim: [-25, 25], dticks: [-20, -10, 0, 10, 20],
    xlab: "pseudotime", ylab: "expression", dlab: "slope (per unit pseudotime)",
    step: 0.25, sel: 22, probe: 0.6, show: ["natural", "smooth"],
    truth: mPt, dtruth: dPt, truthLabel: "true mean",
  },
];

// ---- formatting ----
function fmt(v: number, d = 4): string {
  if (!Number.isFinite(v)) return "–";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e5 || a < 1e-3) { const [m, e] = v.toExponential(d - 2).split("e"); return `${m}e${Number(e)}`; }
  return String(+v.toPrecision(d));
}
const fmtT = (t: number): string => String(+t.toFixed(3));
const sup: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const pow10 = (k: number): string => k === 0 ? "1" : "10" + String(k).split("").map(c => sup[c]).join("");

// λ slider: log10 of R's lambda (t rescaled to [0, 1]); Python's lam = lambda · (range t)³
const LL0 = -9, LL1 = 1, LSTEP = 0.05;

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const uid = ++mountCount;
  const startSet = SETS.find(s => s.key === ds.set) ?? SETS[0];
  const showAttr = ds.show ? ds.show.split(/[\s,]+/).filter((m): m is Method => (METHODS as string[]).includes(m)) : null;
  const st = {
    set: startSet,
    y: startSet.y.slice(),
    show: new Set<Method>(showAttr ?? startSet.show),
    lower: (ds.lower === "change" ? "change" : "slope") as "slope" | "change",
    ll: 0,
    sel: startSet.sel,
    probe: Number.isFinite(Number(ds.probe)) && ds.probe ? Number(ds.probe) : startSet.probe,
  };

  // ---- controls ----
  const row1 = htmlEl("div", "sl-row", "", host);
  const g1 = htmlEl("div", "sl-group", "", row1);
  htmlEl("span", "", "data:", g1);
  const setBtns = presets(g1, SETS.map(s => ({ label: s.label, apply: () => loadSet(s, true) })));
  button(row1, "reset points", () => { st.y = st.set.y.slice(); draw(); });

  const row2 = htmlEl("div", "sl-row", "", host);
  const checks = {} as Record<Method, HTMLInputElement>;
  for (const m of METHODS) {
    checks[m] = checkbox(row2, LABEL[m], st.show.has(m), v => { if (v) st.show.add(m); else st.show.delete(m); draw(); });
    (checks[m].parentElement as HTMLElement).style.color = COL[m];
  }

  const row3 = htmlEl("div", "sl-row", "", host);
  const lamS = slider(row3, "smoothing log₁₀ λ", LL0, LL1, LSTEP, st.ll, v => { st.ll = v; draw(); }, v => v.toFixed(2));
  lamS.input.addEventListener("input", () => { st.ll = lamS.get(); });
  button(row3, "GCV", () => { st.ll = gcvBest(); lamS.set(st.ll); draw(); });
  const g3 = htmlEl("div", "sl-group", "", row3);
  htmlEl("span", "", "panel (b):", g3);
  const lowBtns = presets(g3, [
    { label: "slope", apply: () => { st.lower = "slope"; draw(); } },
    { label: "change", apply: () => { st.lower = "change"; draw(); } },
  ]);
  const g4 = htmlEl("div", "sl-group", "", row3);
  htmlEl("span", "", "point:", g4);
  button(g4, "◀", () => { st.sel = Math.max(0, st.sel - 1); draw(); });
  button(g4, "▶", () => { st.sel = Math.min(st.set.x.length - 1, st.sel + 1); draw(); });
  const minusB = button(g4, "−", () => nudge(-1));
  const plusB = button(g4, "+", () => nudge(1));

  // ---- plots ----
  const plots = htmlEl("div", "sl-plots", "", host);
  const H = 330, WA = 640, WB = 500, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Data and interpolating curves" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Slope or change of the curves" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro = htmlEl("div", "qb-readout", "", host);

  const top = 36, bot = H - 50, aL = 66, aR = WA - 12, bL = 74, bR = WB - 12;
  const layerA = svgEl("g", {}, svgA);
  const ptLayer = svgEl("g", {}, svgA);
  for (const [svg, l, r] of [[svgA, aL, aR], [svgB, bL, bR]] as [SVGSVGElement, number, number][]) {
    const cp = svgEl("clipPath", { id: `sl-clip-${uid}-${svg === svgA ? "a" : "b"}` }, svg);
    svgEl("rect", { x: l, y: top, width: r - l, height: bot - top }, cp);
  }
  const ax = () => scale(st.set.xlim[0], st.set.xlim[1], aL, aR);
  const ay = () => scale(st.set.ylim[0], st.set.ylim[1], bot, top);
  const bx = () => scale(st.set.xlim[0], st.set.xlim[1], bL, bR);

  let circles: SVGCircleElement[] = [];
  function buildPoints(): void {
    ptLayer.replaceChildren();
    circles = st.set.x.map((_, i) => {
      const c = svgEl("circle", { r: st.set.x.length > 20 ? 7 : 9, fill: COL.pt, stroke: "#fff", "stroke-width": 1.5 }, ptLayer);
      c.style.cursor = "ns-resize";
      c.style.touchAction = "none";
      c.addEventListener("pointerdown", ev => {
        ev.preventDefault();
        ev.stopPropagation();
        c.setPointerCapture(ev.pointerId);
        st.sel = i; draw();
        const move = (e: PointerEvent) => {
          const [, py] = svgPoint(svgA, e);
          const [lo, hi] = st.set.ylim;
          const v = Math.max(lo, Math.min(hi, ay().invert(py)));
          st.y[i] = Math.round(v / (st.set.step / 50)) * (st.set.step / 50);
          draw();
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
  // probe: click the background of either panel
  const setProbe = (svg: SVGSVGElement, sx: ReturnType<typeof scale>, ev: PointerEvent) => {
    const [px] = svgPoint(svg, ev);
    const t = sx.invert(px);
    const [x0, x1] = [st.set.x[0], st.set.x[st.set.x.length - 1]];
    st.probe = Math.max(x0, Math.min(x1, Math.round(t * 100) / 100));
    draw();
  };
  svgA.addEventListener("pointerdown", ev => setProbe(svgA, ax(), ev));
  svgB.addEventListener("pointerdown", ev => setProbe(svgB, bx(), ev));
  svgA.style.cursor = svgB.style.cursor = "crosshair";

  function nudge(sign: number): void {
    const [lo, hi] = st.set.ylim;
    st.y[st.sel] = Math.max(lo, Math.min(hi, +(st.y[st.sel] + sign * st.set.step).toFixed(6)));
    draw();
  }

  const range3 = () => (st.set.x[st.set.x.length - 1] - st.set.x[0]) ** 3;
  const lamOf = (ll: number) => 10 ** ll * range3();     // Python units (t as given)
  function gcvBest(): number {
    let best = Infinity, arg = 0;
    for (let k = 0; k <= Math.round((LL1 - LL0) / LSTEP); k++) {
      const ll = LL0 + k * LSTEP;
      const g = smoothing(st.set.x, st.y, lamOf(ll)).gcv;
      if (g < best) { best = g; arg = ll; }
    }
    // refine between the neighbouring grid points by golden-section search
    const g = (ll: number) => smoothing(st.set.x, st.y, lamOf(ll)).gcv;
    let a = Math.max(LL0, arg - LSTEP), b = Math.min(LL1, arg + LSTEP);
    const r = (Math.sqrt(5) - 1) / 2;
    let c = b - r * (b - a), d = a + r * (b - a), gc = g(c), gd = g(d);
    for (let k = 0; k < 30; k++) {
      if (gc < gd) { b = d; d = c; gd = gc; c = b - r * (b - a); gc = g(c); }
      else { a = c; c = d; gc = gd; d = a + r * (b - a); gd = g(d); }
    }
    return (a + b) / 2;
  }

  function loadSet(s: DataSet, resetShow: boolean): void {
    st.set = s;
    st.y = s.y.slice();
    st.sel = s.sel;
    st.probe = s.probe;
    if (resetShow) {
      st.show = new Set(s.show);
      for (const m of METHODS) checks[m].checked = st.show.has(m);
    }
    st.ll = gcvBest();
    lamS.set(st.ll);
    buildPoints();
    draw();
  }

  function curves(y: number[]): Partial<Record<Method, Curve>> & { smoothFit?: SmoothFit } {
    const x = st.set.x, out: Partial<Record<Method, Curve>> & { smoothFit?: SmoothFit } = {};
    for (const m of st.show) {
      if (m === "linear") out.linear = linear(x, y);
      else if (m === "natural") out.natural = fromM(x, y, naturalM(x, y));
      else if (m === "notaknot") out.notaknot = fromM(x, y, notAKnotM(x, y));
      else if (m === "monotone") out.monotone = pchip(x, y);
      else if (m === "poly") out.poly = polynomial(x, y);
      else { const f = smoothing(x, y, lamOf(st.ll)); out.smooth = f; out.smoothFit = f; }
    }
    return out;
  }

  // axes helper
  function axes(svg: Element, sx: (v: number) => number, sy: (v: number) => number, L: number, Rr: number,
    yt: number[], ylab: (v: number) => string, title: string, xlab: string, ytitle: string): void {
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svg).textContent = title;
    for (const v of yt) {
      svgEl("line", { x1: L, x2: Rr, y1: sy(v), y2: sy(v), class: "grid" }, svg);
      svgEl("text", { x: L - 6, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = ylab(v);
    }
    svgEl("line", { x1: L, x2: Rr, y1: bot, y2: bot, class: "axis" }, svg);
    for (const t of st.set.xticks) {
      svgEl("line", { x1: sx(t), x2: sx(t), y1: bot, y2: bot + 5, class: "axis" }, svg);
      svgEl("text", { x: sx(t), y: bot + 22, "text-anchor": "middle", class: "sm" }, svg).textContent = String(t);
    }
    svgEl("text", { x: (L + Rr) / 2, y: H - 6, "text-anchor": "middle" }, svg).textContent = xlab;
    const yc = (top + bot) / 2;
    svgEl("text", { x: 16, y: yc, "text-anchor": "middle", transform: `rotate(-90 16 ${yc})` }, svg).textContent = ytitle;
  }

  function path(ts: number[], vs: number[], sx: (v: number) => number, sy: (v: number) => number): string {
    let d = "", pen = false;
    ts.forEach((t, i) => {
      const v = vs[i];
      if (!Number.isFinite(v)) { pen = false; return; }
      const py = Math.max(-2000, Math.min(H + 2000, sy(v)));
      d += (pen ? "L" : "M") + sx(t).toFixed(1) + "," + py.toFixed(1);
      pen = true;
    });
    return d;
  }

  function draw(): void {
    const s = st.set, x = s.x, n = x.length;
    setBtns.forEach(b => b.classList.toggle("active", b.textContent === s.label));
    lowBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.lower === "slope")));
    const sx = ax(), sy = ay(), sbx = bx();
    const x0 = x[0], x1 = x[n - 1];
    const NG = 700;
    const grid = Array.from({ length: NG + 1 }, (_, i) => x0 + (x1 - x0) * i / NG);
    // add the knots, so piecewise curves are drawn exactly
    const ts = Array.from(new Set([...grid, ...x])).sort((a, b) => a - b);
    const cur = curves(st.y);
    const base = curves(s.y);
    const moved = st.y.map((v, i) => Math.abs(v - s.y[i]) > 1e-12);
    const anyMoved = moved.some(Boolean);
    const order = METHODS.filter(m => st.show.has(m));
    // draw the smoothing spline last (on top), the polynomial first (underneath)
    const drawOrder = ["poly", "linear", "monotone", "notaknot", "natural", "smooth"].filter(m => st.show.has(m as Method)) as Method[];

    // ---------- (a) ----------
    layerA.replaceChildren();
    axes(layerA, sx, sy, aL, aR, s.yticks, v => String(v), "(a) data and curves", s.xlab, s.ylab);
    if (s.ylim[0] < 0) {
      svgEl("rect", { x: aL, y: sy(0), width: aR - aL, height: bot - sy(0), fill: COL.neg, "fill-opacity": 0.6 }, layerA);
      svgEl("line", { x1: aL, x2: aR, y1: sy(0), y2: sy(0), stroke: "#999", "stroke-width": 1.2 }, layerA);
    }
    const ca = svgEl("g", { "clip-path": `url(#sl-clip-${uid}-a)` }, layerA);
    if (s.truth) {
      const tv = ts.map(s.truth);
      svgEl("path", { d: path(ts, tv, sx, sy), fill: "none", stroke: COL.ref, "stroke-width": 2, "stroke-dasharray": "7 5" }, ca);
    }
    // the probe
    svgEl("line", { x1: sx(st.probe), x2: sx(st.probe), y1: top, y2: bot, stroke: "#555", "stroke-width": 1.2, "stroke-dasharray": "2 4" }, ca);
    for (const m of drawOrder) {
      const c = cur[m]!;
      if (anyMoved && st.lower === "change") {
        const b = base[m]!;
        svgEl("path", { d: path(ts, ts.map(t => b.f(t)), sx, sy), fill: "none", stroke: COL[m], "stroke-width": 1.5, "stroke-opacity": 0.4, "stroke-dasharray": "4 4" }, ca);
      }
      svgEl("path", { d: path(ts, ts.map(t => c.f(t)), sx, sy), fill: "none", stroke: COL[m], "stroke-width": m === "smooth" ? 3.2 : 2.6 }, ca);
    }
    // ghost markers at the preset positions of moved points
    moved.forEach((mv, i) => {
      if (mv) svgEl("circle", { cx: sx(x[i]), cy: sy(s.y[i]), r: 6, fill: "none", stroke: COL.ref, "stroke-width": 1.5, "stroke-dasharray": "3 2" }, layerA);
    });
    svgEl("text", { x: sx(st.probe) + 4, y: top + 16, class: "sm" }, layerA).textContent = "t*";
    // legend for the truth
    if (s.truth) {
      const lx = s.key === "pseudo" ? aL + 12 : aR - 190, ly = s.key === "pseudo" ? top + 14 : bot - 14;
      svgEl("line", { x1: lx, x2: lx + 26, y1: ly - 6, y2: ly - 6, stroke: COL.ref, "stroke-width": 2, "stroke-dasharray": "7 5" }, layerA);
      svgEl("text", { x: lx + 32, y: ly, class: "sm" }, layerA).textContent = s.truthLabel ?? "truth";
    }
    circles.forEach((c, i) => {
      c.setAttribute("cx", String(sx(x[i])));
      c.setAttribute("cy", String(sy(Math.max(s.ylim[0], Math.min(s.ylim[1], st.y[i])))));
      const isSel = i === st.sel;
      c.setAttribute("fill", isSel ? COL.sel : moved[i] ? "#555" : COL.pt);
      c.setAttribute("stroke", isSel ? "#000" : "#fff");
      c.setAttribute("stroke-width", isSel ? "2.5" : "1.5");
    });

    // ---------- (b) ----------
    svgB.replaceChildren(svgB.querySelector("clipPath")!);
    const cb = () => svgEl("g", { "clip-path": `url(#sl-clip-${uid}-b)` }, svgB);
    if (st.lower === "slope") {
      const sdy = scale(s.dlim[0], s.dlim[1], bot, top);
      axes(svgB, sbx, sdy, bL, bR, s.dticks, v => String(v), "(b) slope S′(t)", s.xlab, s.dlab);
      svgEl("line", { x1: bL, x2: bR, y1: sdy(0), y2: sdy(0), stroke: "#999", "stroke-width": 1.2 }, svgB);
      const g = cb();
      if (s.dtruth) svgEl("path", { d: path(ts, ts.map(s.dtruth), sbx, sdy), fill: "none", stroke: COL.ref, "stroke-width": 2, "stroke-dasharray": "7 5" }, g);
      svgEl("line", { x1: sbx(st.probe), x2: sbx(st.probe), y1: top, y2: bot, stroke: "#555", "stroke-width": 1.2, "stroke-dasharray": "2 4" }, g);
      for (const m of drawOrder) {
        const c = cur[m]!;
        // piecewise linear: a staircase, drawn with breaks at the knots
        const tt = m === "linear" ? ts.flatMap(t => x.includes(t) && t !== x0 && t !== x1 ? [t - 1e-9 * (x1 - x0), t] : [t]) : ts;
        svgEl("path", { d: path(tt, tt.map(t => c.d(t)), sbx, sdy), fill: "none", stroke: COL[m], "stroke-width": m === "smooth" ? 3.2 : 2.4 }, g);
      }
    } else {
      // log10 |change| since the preset data
      const span = Math.max(...s.y) - Math.min(...s.y);
      const e1 = Math.ceil(Math.log10(span)) + 3, e0 = e1 - 10;
      const sly = scale(e0, e1, bot, top);
      const yt: number[] = [];
      for (let k = e0; k <= e1; k += 2) yt.push(k);
      axes(svgB, sbx, sly, bL, bR, yt, v => pow10(v), "(b) change |ΔS(t)| (log)", s.xlab, "|change| since preset");
      const g = cb();
      moved.forEach((mv, i) => {
        if (mv) svgEl("line", { x1: sbx(x[i]), x2: sbx(x[i]), y1: top, y2: bot, stroke: COL.sel, "stroke-width": 1.5, "stroke-dasharray": "5 4" }, g);
      });
      svgEl("line", { x1: sbx(st.probe), x2: sbx(st.probe), y1: top, y2: bot, stroke: "#555", "stroke-width": 1.2, "stroke-dasharray": "2 4" }, g);
      if (!anyMoved) {
        const t1 = svgEl("text", { x: (bL + bR) / 2, y: (top + bot) / 2 - 12, "text-anchor": "middle", class: "sm" }, svgB);
        t1.textContent = "move a point (drag, or select";
        const t2 = svgEl("text", { x: (bL + bR) / 2, y: (top + bot) / 2 + 14, "text-anchor": "middle", class: "sm" }, svgB);
        t2.textContent = "it and press − / +) to see how far";
        const t3 = svgEl("text", { x: (bL + bR) / 2, y: (top + bot) / 2 + 40, "text-anchor": "middle", class: "sm" }, svgB);
        t3.textContent = "the change spreads";
      } else {
        for (const m of drawOrder) {
          const c = cur[m]!, b = base[m]!;
          const vv = ts.map(t => { const d = Math.abs(c.f(t) - b.f(t)); return d > 0 ? Math.max(e0 - 1, Math.log10(d)) : e0 - 1; });
          svgEl("path", { d: path(ts, vv, sbx, sly), fill: "none", stroke: COL[m], "stroke-width": m === "smooth" ? 3.2 : 2.4 }, g);
        }
      }
    }

    // ---------- readout ----------
    ro.replaceChildren();
    const info = htmlEl("div", "sl-info", "", ro);
    const dy = st.y[st.sel] - s.y[st.sel];
    htmlEl("span", "", `point ${st.sel + 1} of ${n}: t = ${fmtT(x[st.sel])}, y = ${fmt(st.y[st.sel])}` +
      (Math.abs(dy) > 1e-12 ? ` (moved ${dy > 0 ? "+" : ""}${fmt(dy, 3)})` : ""), info);
    htmlEl("span", "", `probe t* = ${fmtT(st.probe)} (click a panel to move it)`, info);
    const sf = cur.smoothFit;
    minusB.title = plusB.title = `move the selected point by ${s.step}`;
    if (sf) {
      const lamR = 10 ** st.ll;
      htmlEl("span", "", `smoothing: λ = ${fmt(lamR, 3)} (R) = ${fmt(lamOf(st.ll), 3)} (Python), df = ${sf.df.toFixed(2)}, GCV = ${fmt(sf.gcv, 4)}`, info);
    }
    const fine = Array.from({ length: 4001 }, (_, i) => x0 + (x1 - x0) * i / 4000);
    if (order.length) {
      const tab = htmlEl("table", "sl-tab", "", ro);
      const hr = htmlEl("tr", "", "", tab);
      const showTrue = !!s.truth;
      for (const h of ["curve", "S(t*)", "S′(t*)", "change at t*", "min S", "max S′ (at t)"]) htmlEl("th", "", h, hr);
      for (const m of order) {
        const c = cur[m]!, b = base[m]!;
        let mn = Infinity, mx = -Infinity, tmx = 0;
        for (const t of fine) {
          const v = c.f(t); if (v < mn) mn = v;
          const d = c.d(t); if (d > mx) { mx = d; tmx = t; }
        }
        const tr = htmlEl("tr", "", "", tab);
        const td0 = htmlEl("td", "", LABEL[m], tr); td0.style.color = COL[m];
        htmlEl("td", "", fmt(c.f(st.probe)), tr);
        htmlEl("td", "", fmt(c.d(st.probe)), tr);
        htmlEl("td", "", anyMoved ? fmt(c.f(st.probe) - b.f(st.probe), 3) : "–", tr);
        htmlEl("td", mn < 0 && s.key === "gene" ? "sl-warn" : "", fmt(mn), tr);
        htmlEl("td", "", `${fmt(mx)} (${fmtT(tmx)})`, tr);
      }
      if (showTrue && s.truth && s.dtruth) {
        let mx = -Infinity, tmx = 0;
        for (const t of fine) { const d = s.dtruth(t); if (d > mx) { mx = d; tmx = t; } }
        const tr = htmlEl("tr", "", "", tab);
        const td0 = htmlEl("td", "", s.truthLabel ?? "truth", tr); td0.style.color = COL.ref;
        htmlEl("td", "", fmt(s.truth(st.probe)), tr);
        htmlEl("td", "", fmt(s.dtruth(st.probe)), tr);
        htmlEl("td", "", "", tr);
        htmlEl("td", "", "", tr);
        htmlEl("td", "", `${fmt(mx)} (${fmtT(tmx)})`, tr);
      }
    } else htmlEl("div", "", "Tick at least one curve.", ro);
  }

  loadSet(startSet, !showAttr);
  if (showAttr) for (const m of METHODS) checks[m].checked = st.show.has(m);
}

autoMount(NAME, mount);

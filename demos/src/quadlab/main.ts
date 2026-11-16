// Quadrature lab for the quadrature unit (Numerical integration).
//  Tab "Rules": an integrand (the lecture's oral PK model on [0, 6] h, a chromatography
//    peak with sigma = 3 s, a constant-rate infusion with a kink at 1 h, sqrt(x) on [0, 1])
//    drawn with the panels of the chosen rule (left Riemann, midpoint, trapezoid, Simpson,
//    Gauss-Legendre), an n slider, the error, and a log-log plot of |error| against n for
//    n = 2 ... 1024 with reference slopes 2 and 4. A fifth integrand is the 10 blood
//    samples of the lecture: linear trapezoid against linear-up/log-down, per-interval
//    errors, and AUC to infinity.
//  Tab "Adaptive": adaptive Simpson (the lecture's adapt, tol 1e-6) or a QUADPACK-style
//    21-point Gauss-Kronrod rule that bisects the worst piece (R's integrate defaults:
//    rel.tol = abs.tol = 1.22e-4, 100 subdivisions), run step by step on the PK model on
//    [0, 24], a chromatogram peak (position and width sliders) on [0, 600] s, or dnorm on
//    [0, L]. Initial pieces and "split at 10" show the fix.
// Auto-mounts on <div data-demo="quadlab">. Optional attributes:
//   data-tab="rules|adaptive"  data-fn="pk|peak|kink|sqrt|data"  data-rule="left|mid|trap|simp|gauss"
//   data-n="4"  data-afn="pk|chrom|dnorm"  data-method="simpson|gk"  data-l="20000"
//   data-pos="340"  data-sigma="3"  data-tabs="1|0" (0 hides the tab buttons)
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, autoMount } from "../lib/ui";

const NAME = "quadlab";
// Okabe-Ito, as in make-quadrature-figs.R
const COL = {
  orange: "#E69F00", sky: "#56B4E9", green: "#009E73", blue: "#0072B2",
  red: "#D55E00", purple: "#CC79A7", grey: "#999999", ink: "#222", ref: "#737373",
};
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .ql-row { display: flex; flex-wrap: wrap; gap: .35em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .ql-group { display: flex; gap: .4em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .ql-row .qb-slider { min-width: 260px; flex: 1 1 260px; max-width: 440px; grid-template-columns: auto 1fr 4.2em; }
[data-demo="${NAME}"] .ql-plots { display: flex; flex-wrap: wrap; gap: 8px; width: 100%; }
[data-demo="${NAME}"] .ql-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .ql-read { display: flex; flex-wrap: wrap; gap: .15em 1.3em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .ql-read span b { font-weight: 600; }
[data-demo="${NAME}"] .ql-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .ql-ok { color: #00704f; font-weight: 600; }
[data-demo="${NAME}"] .ql-hide { display: none !important; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; padding: .3em .5em; }
.reveal [data-demo="${NAME}"] .qb-btn { font-size: .78em; padding: .2em .55em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---------------------------------------------------------------- the lecture's numbers
const DOSE = 500, VOL = 50, KA = 1.5, KE = 0.2;
const pk = (t: number): number => DOSE / VOL * KA / (KA - KE) * (Math.exp(-KE * t) - Math.exp(-KA * t));
const pkInt = (a: number, b: number): number => {
  const F = (t: number) => DOSE / VOL * KA / (KA - KE) * (-Math.exp(-KE * t) / KE + Math.exp(-KA * t) / KA);
  return F(b) - F(a);
};
const T_H = [0, 0.5, 1, 2, 3, 4, 6, 8, 12, 24];
const CONC = [0, 4.99, 6.87, 7.16, 6.20, 5.16, 3.47, 2.33, 1.05, 0.09];
// constant-rate infusion: 500 mg over 1 h, V = 50 L, ke = 0.2/h; kink at t = 1 h
const T_INF = 1, R_INF = 500;
const C_END = R_INF / (VOL * KE) * (1 - Math.exp(-KE * T_INF));
const infusion = (t: number): number => t <= T_INF
  ? R_INF / (VOL * KE) * (1 - Math.exp(-KE * t))
  : C_END * Math.exp(-KE * (t - T_INF));
const infInt = (b: number): number => R_INF / (VOL * KE) * (T_INF - (1 - Math.exp(-KE * T_INF)) / KE)
  + C_END * (1 - Math.exp(-KE * (b - T_INF))) / KE;
const SQ2PI = Math.sqrt(2 * Math.PI);
const dnorm = (x: number): number => Math.exp(-x * x / 2) / SQ2PI;

// ---------------------------------------------------------------- Gauss-Legendre nodes
const glCache = new Map<number, { x: number[]; w: number[] }>();
export function gaussLegendre(n: number): { x: number[]; w: number[] } {
  const c = glCache.get(n);
  if (c) return c;
  const x = new Array<number>(n), w = new Array<number>(n);
  const m = Math.floor((n + 1) / 2);
  for (let i = 0; i < m; i++) {
    let z = Math.cos(Math.PI * (i + 0.75) / (n + 0.5)), dp = 0;
    for (let it = 0; it < 100; it++) {
      let p0 = 1, p1 = z;
      for (let k = 2; k <= n; k++) { const p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
      dp = n * (z * p1 - p0) / (z * z - 1);
      const dz = p1 / dp;
      z -= dz;
      if (Math.abs(dz) < 1e-16) break;
    }
    let p0 = 1, p1 = z;
    for (let k = 2; k <= n; k++) { const p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
    dp = n * (z * p1 - p0) / (z * z - 1);
    x[i] = -z; x[n - 1 - i] = z;
    w[i] = w[n - 1 - i] = 2 / ((1 - z * z) * dp * dp);
  }
  const r = { x, w };
  glCache.set(n, r);
  return r;
}

// ---------------------------------------------------------------- tab (a): integrands and rules
type FnKey = "pk" | "peak" | "kink" | "sqrt";
type RuleKey = "left" | "mid" | "trap" | "simp" | "gauss";
interface Integrand { label: string; f: (x: number) => number; a: number; b: number; exact: number; ymax: number; xlab: string; ylab: string; xt: number[] }
export const FNS: Record<FnKey, Integrand> = {
  pk: { label: "PK curve, 0–6 h", f: pk, a: 0, b: 6, exact: pkInt(0, 6), ymax: 8, xlab: "time after dose (h)", ylab: "mg/L", xt: [0, 1, 2, 3, 4, 5, 6] },
  peak: { label: "peak, σ = 3 s", f: t => 1000 * Math.exp(-((t - 300) ** 2) / 18), a: 270, b: 330, exact: 3000 * SQ2PI, ymax: 1100, xlab: "retention time (s)", ylab: "signal", xt: [270, 285, 300, 315, 330] },
  kink: { label: "infusion (kink at 1 h)", f: infusion, a: 0, b: 6, exact: infInt(6), ymax: 10, xlab: "time (h); 1 h infusion", ylab: "mg/L", xt: [0, 1, 2, 3, 4, 5, 6] },
  sqrt: { label: "√x on [0, 1]", f: Math.sqrt, a: 0, b: 1, exact: 2 / 3, ymax: 1.1, xlab: "x", ylab: "√x", xt: [0, 0.25, 0.5, 0.75, 1] },
};
const RULES: Record<RuleKey, { label: string; col: string; order: string }> = {
  left: { label: "left Riemann", col: COL.grey, order: "1" },
  mid: { label: "midpoint", col: COL.green, order: "2" },
  trap: { label: "trapezoid", col: COL.blue, order: "2" },
  simp: { label: "Simpson", col: COL.red, order: "4" },
  gauss: { label: "Gauss–Legendre", col: COL.purple, order: "2n−1 (degree)" },
};
const RULE_KEYS: RuleKey[] = ["left", "mid", "trap", "simp", "gauss"];

export function applyRule(r: RuleKey, F: Integrand, n: number): number {
  const { f, a, b } = F, h = (b - a) / n;
  let s = 0;
  switch (r) {
    case "left": for (let i = 0; i < n; i++) s += f(a + i * h); return h * s;
    case "mid": for (let i = 0; i < n; i++) s += f(a + (i + 0.5) * h); return h * s;
    case "trap": for (let i = 1; i < n; i++) s += f(a + i * h); return h * (s + (f(a) + f(b)) / 2);
    case "simp": for (let i = 1; i < n; i++) s += (i % 2 ? 4 : 2) * f(a + i * h); return h / 3 * (s + f(a) + f(b));
    case "gauss": {
      const g = gaussLegendre(n), c = (a + b) / 2, hl = (b - a) / 2;
      for (let i = 0; i < n; i++) s += g.w[i] * f(c + hl * g.x[i]);
      return hl * s;
    }
  }
}
const NS = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];
const errCache = new Map<string, number[]>();
function errCurve(fk: FnKey, r: RuleKey): number[] {
  const key = fk + r;
  let e = errCache.get(key);
  if (!e) { e = NS.map(n => Math.abs(applyRule(r, FNS[fk], n) - FNS[fk].exact)); errCache.set(key, e); }
  return e;
}

// AUC from samples
const aucLin = (i: number): number => (T_H[i + 1] - T_H[i]) * (CONC[i] + CONC[i + 1]) / 2;
const isLogDown = (i: number): boolean => CONC[i + 1] < CONC[i] && CONC[i + 1] > 0;
const aucLul = (i: number): number => isLogDown(i)
  ? (T_H[i + 1] - T_H[i]) * (CONC[i] - CONC[i + 1]) / Math.log(CONC[i] / CONC[i + 1]) : aucLin(i);
const NINT = T_H.length - 1;
const sumOver = (g: (i: number) => number): number => { let s = 0; for (let i = 0; i < NINT; i++) s += g(i); return s; };
// terminal slope: straight line through log C of the last three samples
export const LZ = (() => {
  const xs = T_H.slice(-3), ys = CONC.slice(-3).map(Math.log);
  const mx = xs.reduce((p, q) => p + q) / 3, my = ys.reduce((p, q) => p + q) / 3;
  let sxy = 0, sxx = 0;
  for (let i = 0; i < 3; i++) { sxy += (xs[i] - mx) * (ys[i] - my); sxx += (xs[i] - mx) ** 2; }
  return -sxy / sxx;
})();

// ---------------------------------------------------------------- formatting
const sup: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const pow10 = (k: number): string => k === 0 ? "1" : "10" + String(k).split("").map(c => sup[c]).join("");
function sci(v: number, d = 2): string {
  if (v === 0) return "0";
  if (!Number.isFinite(v)) return String(v);
  const av = Math.abs(v);
  if (av >= 1e-3 && av < 1e5) return v.toPrecision(d);
  const [m, e] = v.toExponential(d - 1).split("e");
  return `${m}e${Number(e)}`;
}
const fix = (v: number, d: number): string => {
  if (v === 0) return "0";
  const av = Math.abs(v);
  if (av < 1e-4 || av >= 1e7) return sci(v, d);
  return v.toPrecision(d);
};

// ---------------------------------------------------------------- shared plotting helpers
function axesFrame(svg: SVGSVGElement, x0: number, x1: number, y0: number, y1: number,
  xt: number[], yt: number[], xlab: string, ylab: string, H: number,
  sx: (v: number) => number, sy: (v: number) => number, xfmt: (v: number) => string = String, yfmt: (v: number) => string = String): void {
  for (const v of yt) {
    svgEl("line", { x1: x0, x2: x1, y1: sy(v), y2: sy(v), class: "grid" }, svg);
    svgEl("text", { x: x0 - 6, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = yfmt(v);
  }
  svgEl("line", { x1: x0, x2: x1, y1: y0, y2: y0, class: "axis" }, svg);
  for (const v of xt) {
    svgEl("line", { x1: sx(v), x2: sx(v), y1: y0, y2: y0 + 5, class: "axis" }, svg);
    svgEl("text", { x: sx(v), y: y0 + 22, "text-anchor": "middle", class: "sm" }, svg).textContent = xfmt(v);
  }
  svgEl("text", { x: (x0 + x1) / 2, y: H - 6, "text-anchor": "middle" }, svg).textContent = xlab;
  const yl = svgEl("text", { x: 16, y: (y0 + y1) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(y0 + y1) / 2})` }, svg);
  yl.textContent = ylab;
}
const niceTicks = (lo: number, hi: number, k = 5): number[] => {
  const span = hi - lo, raw = span / k, p = 10 ** Math.floor(Math.log10(raw));
  const st = [1, 2, 2.5, 5, 10].map(m => m * p).find(s => span / s <= k + 0.5) ?? 10 * p;
  const out: number[] = [];
  for (let v = Math.ceil(lo / st - 1e-9) * st; v <= hi + 1e-9 * span; v += st) out.push(Math.abs(v) < 1e-12 * st ? 0 : +v.toPrecision(12));
  return out;
};
const tickFmt = (v: number): string => Math.abs(v) >= 1e4 ? String(Math.round(v)) : String(+v.toPrecision(6));

let mountCount = 0;

// ================================================================== mount
function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const uid = ++mountCount;
  const rulesBox = htmlEl("div", "", ""), adaptBox = htmlEl("div", "", "");
  for (const b of [rulesBox, adaptBox]) Object.assign(b.style, { display: "flex", flexDirection: "column", gap: ".4em", width: "100%" });
  const tabRow = htmlEl("div", "ql-row", "", host);
  const tb = presets(tabRow, [
    { label: "(a) Rules and order", apply: () => show("rules") },
    { label: "(b) Adaptive quadrature", apply: () => show("adaptive") },
  ]);
  if (ds.tabs === "0") tabRow.classList.add("ql-hide");
  host.append(rulesBox, adaptBox);
  const show = (t: string) => {
    rulesBox.style.display = t === "rules" ? "flex" : "none";
    adaptBox.style.display = t === "adaptive" ? "flex" : "none";
    tb.forEach((b, i) => b.classList.toggle("active", (i === 0) === (t === "rules")));
  };
  buildRules(rulesBox, ds, uid);
  buildAdaptive(adaptBox, ds, uid);
  show(ds.tab === "adaptive" ? "adaptive" : "rules");
}

// ================================================================== tab (a)
function buildRules(box: HTMLElement, ds: DOMStringMap, uid: number): void {
  type Mode = FnKey | "data";
  const startFn: Mode = (["pk", "peak", "kink", "sqrt", "data"] as const).find(k => k === ds.fn) ?? "pk";
  const startRule: RuleKey = RULE_KEYS.find(k => k === ds.rule) ?? "trap";
  const nIn = Math.round(Number(ds.n));
  const startN = Number.isFinite(nIn) && nIn >= 2 ? Math.min(64, 2 * Math.round(nIn / 2)) : 4;
  const st = { fn: startFn as Mode, rule: startRule, n: startN, lul: false, logy: false };

  const r1 = htmlEl("div", "ql-row", "", box);
  const g1 = htmlEl("div", "ql-group", "", r1);
  htmlEl("span", "", "integrand:", g1);
  const fnKeys: Mode[] = ["pk", "peak", "kink", "sqrt", "data"];
  const fnBtns = presets(g1, fnKeys.map(k => ({
    label: k === "data" ? "PK: 10 samples" : FNS[k].label,
    apply: () => { st.fn = k; draw(); },
  })));

  const r2 = htmlEl("div", "ql-row", "", box);
  const g2 = htmlEl("div", "ql-group", "", r2);
  htmlEl("span", "", "rule:", g2);
  const ruleBtns = presets(g2, RULE_KEYS.map(k => ({ label: RULES[k].label, apply: () => { st.rule = k; draw(); } })));
  ruleBtns.forEach((b, i) => { b.style.borderBottom = `4px solid ${RULES[RULE_KEYS[i]].col}`; });
  const nS = slider(r2, "n", 2, 64, 2, st.n, v => { st.n = v; draw(); });
  // data-mode controls
  const g3 = htmlEl("div", "ql-group", "", r2);
  htmlEl("span", "", "rule:", g3);
  const dBtns = presets(g3, [
    { label: "linear trapezoid", apply: () => { st.lul = false; draw(); } },
    { label: "linear-up / log-down", apply: () => { st.lul = true; draw(); } },
  ]);
  const logCb = checkbox(g3, "log scale", st.logy, v => { st.logy = v; draw(); });
  const rst = button(r2, "reset", () => {
    st.fn = startFn; st.rule = startRule; st.n = startN; st.lul = false; st.logy = false;
    nS.set(st.n); logCb.checked = false; draw();
  });
  void rst;

  const plots = htmlEl("div", "ql-plots", "", box);
  const H = 350, WA = 620, WB = 530, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Integrand with the panels of the rule" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 4px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Error against n, log-log" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 4px)`;
  const ro = htmlEl("div", "qb-readout ql-read", "", box);
  const item = (k: string, v: string, cls = "") => {
    const s = htmlEl("span", cls, "", ro);
    if (k) htmlEl("b", "", k + " ", s);
    s.append(v);
  };

  const top = 36, bot = H - 50;

  function draw(): void {
    const data = st.fn === "data";
    fnBtns.forEach((b, i) => b.classList.toggle("active", fnKeys[i] === st.fn));
    ruleBtns.forEach((b, i) => b.classList.toggle("active", RULE_KEYS[i] === st.rule));
    dBtns.forEach((b, i) => b.classList.toggle("active", (i === 1) === st.lul));
    g2.classList.toggle("ql-hide", data);
    nS.input.parentElement!.classList.toggle("ql-hide", data);
    g3.classList.toggle("ql-hide", !data);
    svgA.replaceChildren(); svgB.replaceChildren(); ro.replaceChildren();
    if (data) drawData(); else drawRule(st.fn as FnKey);
  }

  // ---------------- formula mode
  function drawRule(fk: FnKey): void {
    const F = FNS[fk], r = st.rule, n = st.n, col = RULES[r].col;
    const L = 64, R = WA - 12;
    const sx = scale(F.a, F.b, L, R), sy = scale(0, F.ymax, bot, top);
    const yt = niceTicks(0, F.ymax * 0.95, 4);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgA).textContent = `(a) ${RULES[r].label}, n = ${n}`;
    axesFrame(svgA, L, R, bot, top, F.xt, yt, F.xlab, F.ylab, H, sx, sy, tickFmt, tickFmt);
    const fill = COL.sky, edge = COL.blue;
    const g = svgEl("g", { "fill-opacity": 0.4 }, svgA);
    const h = (F.b - F.a) / n, f = F.f;
    const poly = (pts: [number, number][]) => svgEl("polygon", {
      points: pts.map(([x, y]) => `${sx(x)},${sy(y)}`).join(" "), fill, stroke: edge, "stroke-width": n > 32 ? 0.8 : 1.5, "stroke-opacity": 1,
    }, g);
    const dots: [number, number][] = [];
    if (r === "left" || r === "mid") {
      for (let i = 0; i < n; i++) {
        const xa = F.a + i * h, xe = r === "left" ? xa : xa + h / 2, y = f(xe);
        poly([[xa, 0], [xa, y], [xa + h, y], [xa + h, 0]]);
        dots.push([xe, y]);
      }
    } else if (r === "trap") {
      for (let i = 0; i < n; i++) {
        const xa = F.a + i * h;
        poly([[xa, 0], [xa, f(xa)], [xa + h, f(xa + h)], [xa + h, 0]]);
      }
      for (let i = 0; i <= n; i++) dots.push([F.a + i * h, f(F.a + i * h)]);
    } else if (r === "simp") {
      for (let i = 0; i < n; i += 2) {
        const x0 = F.a + i * h, x1 = x0 + h, x2 = x0 + 2 * h, f0 = f(x0), f1 = f(x1), f2 = f(x2);
        const pts: [number, number][] = [[x0, 0]];
        for (let k = 0; k <= 24; k++) {
          const s = k / 12 - 1;   // -1..1
          pts.push([x1 + s * h, f1 + s * (f2 - f0) / 2 + s * s * (f0 - 2 * f1 + f2) / 2]);
        }
        pts.push([x2, 0]);
        poly(pts);
      }
      for (let i = 0; i <= n; i++) dots.push([F.a + i * h, f(F.a + i * h)]);
    } else {
      const gl = gaussLegendre(n), c = (F.a + F.b) / 2, hl = (F.b - F.a) / 2;
      for (let i = 0; i < n; i++) {
        const x = c + hl * gl.x[i], w = hl * gl.w[i], y = f(x);
        poly([[x - w / 2, 0], [x - w / 2, y], [x + w / 2, y], [x + w / 2, 0]]);
        dots.push([x, y]);
      }
    }
    // true curve
    const NP = 600, cp: string[] = [];
    for (let k = 0; k <= NP; k++) { const x = F.a + (F.b - F.a) * k / NP; cp.push(`${sx(x)},${sy(f(x))}`); }
    if (fk === "kink") cp.splice(Math.round(NP * (T_INF - F.a) / (F.b - F.a)), 0, `${sx(T_INF)},${sy(f(T_INF))}`);
    svgEl("polyline", { points: cp.join(" "), fill: "none", stroke: COL.red, "stroke-width": 3 }, svgA);
    const rd = n > 32 ? 2.5 : 4.5;
    for (const [x, y] of dots) svgEl("circle", { cx: sx(x), cy: sy(y), r: rd, fill: COL.ink }, svgA);
    if (r === "gauss") svgEl("text", { x: R, y: top + 20, "text-anchor": "end", class: "sm" }, svgA).textContent = "bar width = weight";
    if (fk === "kink") {
      const t = svgEl("text", { x: sx(T_INF) + 8, y: sy(f(T_INF)) - 8, class: "sm" }, svgA);
      t.textContent = "kink";
    }

    // ---- (b) log-log
    const L2 = 74, R2 = WB - 12;
    const bx = scale(1, 10, L2, R2);                       // log2 n
    const EY0 = -16, EY1 = 4, by = scale(EY0, EY1, bot, top);
    const cl = (e: number) => Math.max(EY0 - 0.3, Math.min(EY1, Math.log10(Math.max(e, 1e-30))));
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgB).textContent = "(b) |error| against n (log–log)";
    const yt2: number[] = []; for (let k = EY0; k <= EY1; k += 4) yt2.push(k);
    axesFrame(svgB, L2, R2, bot, top, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], yt2, "n (intervals; Gauss: points)", "|error|", H,
      bx, by, v => String(2 ** v), v => pow10(v));
    const clip = svgEl("clipPath", { id: `ql-clip-${uid}` }, svgB);
    svgEl("rect", { x: L2, y: top - 6, width: R2 - L2, height: bot - top + 8 }, clip);
    const pg = svgEl("g", { "clip-path": `url(#ql-clip-${uid})` }, svgB);
    // reference slopes, anchored at the trapezoid / Simpson error at n = 8
    const anchor = (rk: RuleKey, p: number) => {
      const e8 = errCurve(fk, rk)[2] || 1e-6;
      const ly = (l2: number) => Math.log10(e8) - p * (l2 - 3) * Math.log10(2);
      const yA = ly(1), yB = ly(10);
      svgEl("line", { x1: bx(1), y1: by(yA + 0.5), x2: bx(10), y2: by(yB + 0.5), stroke: COL.ref, "stroke-width": 2, "stroke-dasharray": "8 6" }, pg);
      return (l2: number) => ly(l2) + 0.5;
    };
    const l2 = anchor("trap", 2), l4 = anchor("simp", 4);
    const lab = (txt: string, x: number, y: number) => {
      const t = svgEl("text", { x, y, class: "sm", "text-anchor": "end" }, pg); t.style.fill = COL.ref; t.textContent = txt;
    };
    const labAt = (f: (l: number) => number, txt: string) => {
      let x = 9.6;
      while (x > 1.5 && f(x) < EY0 + 1.2) x -= 0.1;
      lab(txt, bx(x), by(f(x)) - 10);
    };
    labAt(l2, "slope −2"); labAt(l4, "slope −4");
    // all rules, current one strong
    for (const rk of RULE_KEYS) {
      const e = errCurve(fk, rk), strong = rk === r;
      const pts = NS.map((nn, i) => `${bx(Math.log2(nn))},${by(cl(e[i]))}`).join(" ");
      svgEl("polyline", { points: pts, fill: "none", stroke: RULES[rk].col, "stroke-width": strong ? 3 : 1.6, "stroke-opacity": strong ? 1 : 0.35 }, pg);
      NS.forEach((nn, i) => svgEl("circle", { cx: bx(Math.log2(nn)), cy: by(cl(e[i])), r: strong ? 4 : 2.5, fill: RULES[rk].col, "fill-opacity": strong ? 1 : 0.35 }, pg));
    }
    // current n
    const val = applyRule(r, F, n), err = Math.abs(val - F.exact);
    svgEl("line", { x1: bx(Math.log2(n)), x2: bx(Math.log2(n)), y1: top, y2: bot, stroke: COL.ink, "stroke-width": 1, "stroke-dasharray": "2 3" }, svgB);
    svgEl("circle", { cx: bx(Math.log2(n)), cy: by(cl(err)), r: 8, fill: col, stroke: COL.ink, "stroke-width": 2.5 }, svgB);
    if (err === 0 || err < 1e-16 * Math.abs(F.exact) * 4)
      svgEl("text", { x: R2 - 4, y: bot - 8, "text-anchor": "end", class: "sm" }, svgB).textContent = "bottom: rounding level";

    // readout
    const e2 = Math.abs(applyRule(r, F, 2 * n) - F.exact);
    const order = err > 0 && e2 > 0 ? Math.log2(err / e2) : NaN;
    const evals = r === "left" || r === "mid" || r === "gauss" ? n : n + 1;
    item("estimate", val.toPrecision(10));
    item("exact", F.exact.toPrecision(10));
    item("|error|", sci(err, 2));
    item("rel. error", err / F.exact >= 1e-3 ? (100 * err / F.exact).toPrecision(3) + " %" : sci(err / F.exact, 2));
    item("function values:", String(evals));
    item("order n → 2n: log₂(err(n)/err(2n)) =", Number.isFinite(order) ? order.toFixed(2) : "–");
    if (fk === "kink" && r !== "gauss" && Math.abs(T_INF / h - Math.round(T_INF / h)) < 1e-9)
      item("", "a grid point sits on the kink", "ql-ok");
  }

  // ---------------- data mode: the 10 blood samples
  function drawData(): void {
    const L = 64, R = WA - 12;
    const sx = scale(0, 24, L, R);
    const LY0 = -2, LY1 = 1;
    const sy = st.logy ? (v: number) => scale(LY0, LY1, bot, top)(Math.log10(Math.max(v, 10 ** (LY0 - 0.2))))
      : scale(0, 8.5, bot, top);
    const yt = st.logy ? [0.01, 0.1, 1, 10] : [0, 2, 4, 6, 8];
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgA).textContent =
      `(a) 10 samples, ${st.lul ? "linear-up / log-down" : "linear trapezoid"}`;
    axesFrame(svgA, L, R, bot, top, [0, 2, 4, 6, 8, 12, 24], yt, "time after dose (h)", "concentration (mg/L)", H,
      sx, sy, String, v => String(v));
    const g = svgEl("g", {}, svgA);
    const ybase = st.logy ? bot : sy(0);
    for (let i = 0; i < NINT; i++) {
      const ta = T_H[i], tb = T_H[i + 1], ca = CONC[i], cb = CONC[i + 1];
      const pts: string[] = [`${sx(ta)},${ybase}`];
      const logd = st.lul && isLogDown(i);
      for (let k = 0; k <= 30; k++) {
        const s = k / 30, t = ta + s * (tb - ta);
        const c = logd ? ca * (cb / ca) ** s : ca + s * (cb - ca);
        pts.push(`${sx(t)},${c > 0 || !st.logy ? sy(c) : ybase}`);
      }
      pts.push(`${sx(tb)},${ybase}`);
      svgEl("polygon", { points: pts.join(" "), fill: logd ? COL.green : COL.sky, "fill-opacity": 0.35, stroke: logd ? COL.green : COL.blue, "stroke-width": 1.5 }, g);
    }
    const cp: string[] = [];
    for (let k = 0; k <= 480; k++) { const t = 24 * k / 480; if (!st.logy || t > 0.02) cp.push(`${sx(t)},${sy(pk(t))}`); }
    svgEl("polyline", { points: cp.join(" "), fill: "none", stroke: COL.red, "stroke-width": 3 }, svgA);
    T_H.forEach((t, i) => { if (!st.logy || CONC[i] > 0) svgEl("circle", { cx: sx(t), cy: sy(CONC[i]), r: 6, fill: COL.ink }, svgA); });
    const lg = (txt: string, y: number, c: string) => { const t = svgEl("text", { x: R - 4, y, "text-anchor": "end", class: "sm" }, svgA); t.style.fill = c; t.style.fontWeight = "600"; t.textContent = txt; };
    lg("— model (true curve)", top + 20, COL.red);
    lg("■ straight line", top + 42, COL.blue);
    if (st.lul) lg("■ exponential (log-down)", top + 64, COL.green);

    // ---- (b) per-interval error
    const L2 = 74, R2 = WB - 12;
    const errL = Array.from({ length: NINT }, (_, i) => aucLin(i) - pkInt(T_H[i], T_H[i + 1]));
    const errU = Array.from({ length: NINT }, (_, i) => aucLul(i) - pkInt(T_H[i], T_H[i + 1]));
    const by = scale(-0.8, 2.4, bot, top);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgB).textContent = "(b) error per interval (mg·h/L)";
    const labs = Array.from({ length: NINT }, (_, i) => i);
    const bw = (R2 - L2) / NINT;
    const bxc = (i: number) => L2 + (i + 0.5) * bw;
    axesFrame(svgB, L2, R2, by(0), top, [], [-0.5, 0, 0.5, 1, 1.5, 2], "", "rule − true", H, bxc, by, String, String);
    svgEl("text", { x: (L2 + R2) / 2, y: H - 6, "text-anchor": "middle" }, svgB).textContent = "interval (h)";
    labs.forEach(i => {
      const t = svgEl("text", { x: bxc(i), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgB);
      t.textContent = `${T_H[i]}–${T_H[i + 1]}`;
      t.style.fontSize = "14px";
      const bar = (v: number, dx: number, c: string, strong: boolean) => {
        const y0 = by(0), y1 = by(Math.max(-0.8, Math.min(2.4, v)));
        svgEl("rect", { x: bxc(i) + dx, y: Math.min(y0, y1), width: bw * 0.36, height: Math.max(1, Math.abs(y1 - y0)), fill: c, "fill-opacity": strong ? 0.95 : 0.35 }, svgB);
      };
      bar(errL[i], -bw * 0.38, COL.blue, !st.lul);
      bar(errU[i], bw * 0.02, COL.green, st.lul);
    });
    const big = errL.reduce((m, v, i) => Math.abs(v) > Math.abs(errL[m]) ? i : m, 0);
    const tt = svgEl("text", { x: bxc(big) - 6, y: by(errL[big]) - 8, "text-anchor": "end", class: "sm" }, svgB);
    tt.textContent = `linear: +${errL[big].toFixed(2)}`;
    const lgB = (txt: string, y: number, c: string) => {
      svgEl("rect", { x: L2 + 12, y: y - 14, width: 16, height: 16, fill: c }, svgB);
      const t = svgEl("text", { x: L2 + 34, y }, svgB); t.textContent = txt;
    };
    lgB("linear trapezoid", top + 22, COL.blue);
    lgB("linear-up / log-down", top + 46, COL.green);

    // readout
    const aL = sumOver(aucLin), aU = sumOver(aucLul), tru = pkInt(0, 24);
    const ext = CONC[NINT] / LZ, a = st.lul ? aU : aL;
    item("AUC 0–24: linear", aL.toPrecision(6));
    item("lin-up/log-down", aU.toPrecision(7));
    item("model", tru.toPrecision(7));
    item(st.lul ? "error (log-down)" : "error (linear)", `${a > tru ? "+" : ""}${(100 * (a - tru) / tru).toFixed(1)} %`);
    item("λz (last 3 samples) =", LZ.toFixed(4) + " /h");
    item("tail C_last/λz =", ext.toFixed(3));
    item("AUC 0–∞ =", (a + ext).toPrecision(6));
    item("model Dose/(V·ke) =", (DOSE / VOL / KE).toFixed(0));
  }

  draw();
}

// ================================================================== tab (b): adaptive quadrature
type AFn = "pk" | "chrom" | "dnorm";
type Method = "simpson" | "gk";
const SIMP_TOL = 1e-6;
const GK_TOL = 2 ** -13;          // R's integrate: rel.tol = abs.tol = .Machine$double.eps^0.25
const GK_LIMIT = 100;
// 21-point Gauss-Kronrod (QUADPACK qk21): nodes on [0, 1], Kronrod weights, Gauss weights for odd i
const XGK = [0.995657163025808080735527280689003, 0.973906528517171720077964012084452, 0.930157491355708226001207180059508,
  0.865063366688984510732096688423493, 0.780817726586416897063717578345042, 0.679409568299024406234327365114874,
  0.562757134668604683339000099272694, 0.433395394129247190799265943165784, 0.294392862701460198131126603103866,
  0.148874338981631210884826001129720, 0];
const WGK = [0.011694638867371874278064396062192, 0.032558162307964727478818972459390, 0.054755896574351996031381300244580,
  0.075039674810919952767043140916190, 0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
  0.123491976262065851077208937799534, 0.134709217311473325928054001771707, 0.142775938577060080797094273138717,
  0.147739104901338491374841515972068, 0.149445554002916905664936468389821];
const WG = [0.066671344308688137593568809893332, 0.149451349150580593145776339657697, 0.219086362515982043995534934228163,
  0.269266719309996355091226921569469, 0.295524224714752870173892994651338];
const EPM = 2.220446049250313e-16, UFLOW = 2.2250738585072014e-308;

function qk21(f: (x: number) => number, a: number, b: number, rec: (x: number) => void): { res: number; err: number; asc: number } {
  const c = (a + b) / 2, hl = (b - a) / 2;
  const fc = f(c); rec(c);
  let resk = WGK[10] * fc, resg = 0, resabs = Math.abs(resk);
  const fv1: number[] = [], fv2: number[] = [];
  for (let j = 0; j < 10; j++) {
    const x1 = c - hl * XGK[j], x2 = c + hl * XGK[j];
    const f1 = f(x1), f2 = f(x2); rec(x1); rec(x2);
    fv1.push(f1); fv2.push(f2);
    resk += WGK[j] * (f1 + f2);
    resabs += WGK[j] * (Math.abs(f1) + Math.abs(f2));
    if (j % 2 === 1) resg += WG[(j - 1) / 2] * (f1 + f2);
  }
  const reskh = resk * 0.5;
  let resasc = WGK[10] * Math.abs(fc - reskh);
  for (let j = 0; j < 10; j++) resasc += WGK[j] * (Math.abs(fv1[j] - reskh) + Math.abs(fv2[j] - reskh));
  const res = resk * hl;
  resabs *= Math.abs(hl); resasc *= Math.abs(hl);
  let err = Math.abs((resk - resg) * hl);
  if (resasc !== 0 && err !== 0) err = resasc * Math.min(1, (200 * err / resasc) ** 1.5);
  if (resabs > UFLOW / (50 * EPM)) err = Math.max(EPM * 50 * resabs, err);
  return { res, err, asc: resasc };
}

interface Piece { a: number; b: number; depth: number; val: number; err: number; status: "pending" | "ok" | "limit";
  fa?: number; fm?: number; fb?: number; tol?: number }
interface RunState { pieces: Piece[]; nodes: number[]; done: boolean; note: string; first?: { val: number; err: number }; steps: number }

export function initRun(m: Method, f: (x: number) => number, breaks: number[]): RunState {
  const nodes: number[] = [];
  const rec = (x: number) => { nodes.push(x); };
  const A = breaks[0], B = breaks[breaks.length - 1];
  const pieces: Piece[] = [];
  if (m === "simpson") {
    for (let i = 0; i + 1 < breaks.length; i++) {
      const a = breaks[i], b = breaks[i + 1], mm = (a + b) / 2;
      const fa = i === 0 ? f(a) : pieces[i - 1].fb!;
      if (i === 0) rec(a);
      const fm = f(mm), fb = f(b); rec(mm); rec(b);
      pieces.push({ a, b, depth: 0, fa, fm, fb, val: (b - a) / 6 * (fa + 4 * fm + fb), err: NaN, status: "pending", tol: SIMP_TOL * (b - a) / (B - A) });
    }
    return { pieces, nodes, done: false, note: "", steps: 0 };
  }
  let asc = 0;
  for (let i = 0; i + 1 < breaks.length; i++) {
    const q = qk21(f, breaks[i], breaks[i + 1], rec);
    pieces.push({ a: breaks[i], b: breaks[i + 1], depth: 0, val: q.res, err: q.err, status: "ok" });
    asc += q.asc;
  }
  const st: RunState = { pieces, nodes, done: false, note: "", steps: 0 };
  st.first = { val: total(st), err: errTotal(st) };
  // QUADPACK (dqagse) after the first pass: stop if the error estimate is 0, or if it is
  // below the tolerance and not "saturated" (equal to resasc, i.e. the rule could not judge)
  const e = errTotal(st), bound = Math.max(GK_TOL, GK_TOL * Math.abs(total(st)));
  if (e === 0) { st.done = true; st.note = "error estimate 0 after the first pass: accept"; }
  else if (e <= bound && e !== asc) { st.done = true; st.note = "error estimate ≤ tolerance: accept"; }
  else if (e <= bound) st.note = "estimate saturated, not trusted: bisect";
  return st;
}
export const total = (s: RunState): number => s.pieces.reduce((p, q) => p + q.val, 0);
export const errTotal = (s: RunState): number => s.pieces.reduce((p, q) => p + (Number.isFinite(q.err) ? q.err : 0), 0);
function checkGK(s: RunState): void {
  const bound = Math.max(GK_TOL, GK_TOL * Math.abs(total(s)));
  if (errTotal(s) <= bound) { s.done = true; s.note = "error estimate ≤ tolerance: accept"; }
  else if (s.pieces.length >= GK_LIMIT) { s.done = true; s.note = "maximum number of subdivisions (100) reached"; }
}

/** One step: Simpson processes the leftmost pending piece; Gauss-Kronrod bisects the worst piece. */
export function stepRun(m: Method, f: (x: number) => number, s: RunState): void {
  if (s.done) return;
  s.steps++;
  const rec = (x: number) => { s.nodes.push(x); };
  if (m === "simpson") {
    const i = s.pieces.findIndex(p => p.status === "pending");
    if (i < 0) { s.done = true; s.note = "every piece accepted"; return; }
    const p = s.pieces[i], mm = (p.a + p.b) / 2;
    const fl = f((p.a + mm) / 2), fr = f((mm + p.b) / 2); rec((p.a + mm) / 2); rec((mm + p.b) / 2);
    const S1 = p.val, S2 = (p.b - p.a) / 12 * (p.fa! + 4 * fl + 2 * p.fm! + 4 * fr + p.fb!);
    if (Math.abs(S2 - S1) < 15 * p.tol! || p.depth >= 45 || s.nodes.length > 40000) {
      p.val = S2 + (S2 - S1) / 15; p.err = Math.abs(S2 - S1) / 15;
      p.status = Math.abs(S2 - S1) < 15 * p.tol! ? "ok" : "limit";
    } else {
      const L: Piece = { a: p.a, b: mm, depth: p.depth + 1, fa: p.fa, fm: fl, fb: p.fm, val: (mm - p.a) / 6 * (p.fa! + 4 * fl + p.fm!), err: NaN, status: "pending", tol: p.tol! / 2 };
      const R: Piece = { a: mm, b: p.b, depth: p.depth + 1, fa: p.fm, fm: fr, fb: p.fb, val: (p.b - mm) / 6 * (p.fm! + 4 * fr + p.fb!), err: NaN, status: "pending", tol: p.tol! / 2 };
      s.pieces.splice(i, 1, L, R);
    }
    if (!s.pieces.some(q => q.status === "pending")) { s.done = true; s.note = "every piece accepted"; }
    return;
  }
  s.note = "";
  let w = 0;
  s.pieces.forEach((p, i) => { if (p.err > s.pieces[w].err) w = i; });
  const p = s.pieces[w], mm = (p.a + p.b) / 2;
  const l = qk21(f, p.a, mm, rec), r = qk21(f, mm, p.b, rec);
  s.pieces.splice(w, 1,
    { a: p.a, b: mm, depth: p.depth + 1, val: l.res, err: l.err, status: "ok" },
    { a: mm, b: p.b, depth: p.depth + 1, val: r.res, err: r.err, status: "ok" });
  checkGK(s);
}

/** Accurate reference value: composite 20-point Gauss-Legendre on many pieces where f lives. */
export function reference(f: (x: number) => number, a: number, b: number, lo: number, hi: number): number {
  const parts: [number, number, number][] = [];
  const l = Math.max(a, lo), h = Math.min(b, hi);
  if (l < h) {
    if (l > a) parts.push([a, l, 50]);
    parts.push([l, h, 400]);
    if (h < b) parts.push([h, b, 50]);
  } else parts.push([a, b, 400]);
  const g = gaussLegendre(20);
  let s = 0;
  for (const [u, v, k] of parts) {
    const w = (v - u) / k;
    for (let j = 0; j < k; j++) {
      const c = u + (j + 0.5) * w;
      let t = 0;
      for (let i = 0; i < 20; i++) t += g.w[i] * f(c + w / 2 * g.x[i]);
      s += t * w / 2;
    }
  }
  return s;
}

function buildAdaptive(box: HTMLElement, ds: DOMStringMap, uid: number): void {
  const num = (v: string | undefined, d: number, lo: number, hi: number) => {
    const x = Number(v); return v !== undefined && Number.isFinite(x) ? Math.min(hi, Math.max(lo, x)) : d;
  };
  const start = {
    fn: ((["pk", "chrom", "dnorm"] as const).find(k => k === ds.afn) ?? "dnorm") as AFn,
    m: (ds.method === "simpson" ? "simpson" : "gk") as Method,
    L: num(ds.l, 20000, 10, 20000), pos: num(ds.pos, 340, 0, 600), sig: num(ds.sigma, 3, 0.5, 40),
    k: 1, split: false,
  };
  const st = { ...start };

  // ---- controls
  const r1 = htmlEl("div", "ql-row", "", box);
  const g1 = htmlEl("div", "ql-group", "", r1);
  htmlEl("span", "", "integrand:", g1);
  const fnKeys: AFn[] = ["pk", "chrom", "dnorm"];
  const fnBtns = presets(g1, [
    { label: "PK model, 0–24 h", apply: () => { st.fn = "pk"; recompute(); } },
    { label: "chromatogram peak", apply: () => { st.fn = "chrom"; recompute(); } },
    { label: "dnorm on [0, L]", apply: () => { st.fn = "dnorm"; recompute(); } },
  ]);
  const g2 = htmlEl("div", "ql-group", "", r1);
  htmlEl("span", "", "method:", g2);
  const mBtns = presets(g2, [
    { label: "adaptive Simpson", apply: () => { st.m = "simpson"; recompute(); } },
    { label: "Gauss–Kronrod 21 (integrate/quad)", apply: () => { st.m = "gk"; recompute(); } },
  ]);

  const r2 = htmlEl("div", "ql-row", "", box);
  const lS = slider(r2, "L", 1, Math.log10(20000), 0.005, Math.log10(st.L), v => { st.L = Math.round(10 ** v); recompute(); },
    v => String(Math.round(10 ** v)));
  const lOut = lS.input.parentElement!.querySelector("output")!;
  const setL = (L: number) => { st.L = L; lS.set(Math.log10(L)); lOut.value = String(L); recompute(); };
  const gL = htmlEl("div", "ql-group", "", r2);
  const Ls = [20, 2000, 2300, 4300, 4400, 5000, 20000];
  const lBtns = presets(gL, Ls.map(L => ({ label: String(L), apply: () => setL(L) })));
  const posS = slider(r2, "peak at (s)", 0, 600, 1, st.pos, v => { st.pos = v; recompute(); });
  const sigS = slider(r2, "σ (s)", 0.5, 40, 0.5, st.sig, v => { st.sig = v; recompute(); });

  const r3 = htmlEl("div", "ql-row", "", box);
  const kS = slider(r3, "initial pieces", 1, 60, 1, st.k, v => { st.k = v; recompute(); });
  const splitCb = checkbox(r3, "split at x = 10", st.split, v => { st.split = v; recompute(); });
  const gB = htmlEl("div", "ql-group", "", r3);
  const bStart = button(gB, "start over", () => { stop(); run = initRun(st.m, fcur, breaks); draw(); });
  const bStep = button(gB, "step", () => { stop(); stepRun(st.m, fcur, run); draw(); });
  const bRun = button(gB, "run", () => { if (timer) stop(); else play(); draw(); });
  const bEnd = button(gB, "finish", () => { stop(); finish(run); draw(); });
  button(gB, "reset", () => {
    stop(); Object.assign(st, start);
    lS.set(Math.log10(st.L)); lOut.value = String(st.L); posS.set(st.pos); sigS.set(st.sig); kS.set(st.k); splitCb.checked = st.split;
    recompute();
  });
  void bStart; void bStep; void bEnd;

  const plots = htmlEl("div", "ql-plots", "", box);
  const H = 330, WA = 740, WB = 410, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Integrand, evaluation points and pieces" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 4px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Zoom on the feature" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 4px)`;
  const ro = htmlEl("div", "qb-readout ql-read", "", box);

  // ---- integrand
  let fcur: (x: number) => number = pk, A = 0, B = 24, exact = 0, breaks: number[] = [0, 24];
  let run: RunState;
  let timer: number | undefined;
  let feature: [number, number] = [0, 3];       // where the function lives (dense curve sampling, zoom)
  let zoom: [number, number] = [0, 3];

  function recompute(): void {
    stop();
    const { fn } = st;
    if (fn === "pk") { fcur = pk; A = 0; B = 24; exact = pkInt(0, 24); feature = [0, 24]; zoom = [0, 4]; }
    else if (fn === "chrom") {
      const p = st.pos, s = st.sig;
      fcur = t => 1000 * Math.exp(-((t - p) ** 2) / (2 * s * s));
      A = 0; B = 600; feature = [p - 8 * s, p + 8 * s];
      exact = reference(fcur, A, B, p - 40 * s, p + 40 * s);
      zoom = [Math.max(0, p - 30), Math.min(600, p + 30)];
      if (s > 6) zoom = [Math.max(0, p - 5 * s), Math.min(600, p + 5 * s)];
    } else {
      fcur = dnorm; A = 0; B = st.L; feature = [0, 8];
      exact = reference(dnorm, 0, B, 0, 40);
      zoom = [0, Math.min(B, 60)];
    }
    const bs = Array.from({ length: st.k + 1 }, (_, i) => A + (B - A) * i / st.k);
    if (fn === "dnorm" && st.split && B > 10 && !bs.includes(10)) { bs.push(10); bs.sort((x, y) => x - y); }
    breaks = bs;
    run = initRun(st.m, fcur, breaks);
    finish(run);
    draw();
  }
  function finish(s: RunState): void { let guard = 0; while (!s.done && guard++ < 100000) stepRun(st.m, fcur, s); }
  function stop(): void { if (timer !== undefined) { clearInterval(timer); timer = undefined; } }
  function play(): void {
    if (run.done) run = initRun(st.m, fcur, breaks);
    const probe = initRun(st.m, fcur, breaks); finish(probe);
    const per = Math.max(1, Math.ceil(probe.steps / 70));
    timer = window.setInterval(() => {
      for (let i = 0; i < per && !run.done; i++) stepRun(st.m, fcur, run);
      if (run.done) stop();
      draw();
    }, 110);
  }

  const top = 34, bot = H - 92;

  function panel(svg: SVGSVGElement, W: number, x0: number, x1: number, title: string, main: boolean): void {
    const L = main ? 70 : 62, R = W - (main ? 30 : 22);
    const sx = scale(x0, x1, L, R);
    const ymax = st.fn === "pk" ? 8 : st.fn === "chrom" ? 1100 : 0.44;
    const sy = scale(0, ymax, bot, top);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svg).textContent = title;
    const xt = niceTicks(x0, x1, main ? 6 : 4), yt = niceTicks(0, ymax * 0.95, 4);
    for (const v of yt) {
      svgEl("line", { x1: L, x2: R, y1: sy(v), y2: sy(v), class: "grid" }, svg);
      svgEl("text", { x: L - 6, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = tickFmt(v);
    }
    svgEl("line", { x1: L, x2: R, y1: bot, y2: bot, class: "axis" }, svg);
    for (const v of xt) {
      svgEl("line", { x1: sx(v), x2: sx(v), y1: bot, y2: bot + 5, class: "axis" }, svg);
      svgEl("text", { x: sx(v), y: bot + 21, "text-anchor": "middle", class: "sm" }, svg).textContent = tickFmt(v);
    }
    const clipId = `ql-ac-${uid}-${main ? 1 : 2}`;
    const cp = svgEl("clipPath", { id: clipId }, svg);
    svgEl("rect", { x: L - 2, y: top - 10, width: R - L + 4, height: H - top }, cp);
    const g = svgEl("g", { "clip-path": `url(#${clipId})` }, svg);
    // piece boundaries
    const vis = run.pieces.filter(p => p.b >= x0 && p.a <= x1);
    if (vis.length <= 80) for (const p of vis) {
      if (p.a >= x0) svgEl("line", { x1: sx(p.a), x2: sx(p.a), y1: top, y2: bot, stroke: "#bbb", "stroke-width": 1 }, g);
    }
    // curve: uniform samples plus a dense set over the feature
    const xs: number[] = [];
    for (let k = 0; k <= 700; k++) xs.push(x0 + (x1 - x0) * k / 700);
    const fa = Math.max(x0, feature[0]), fb = Math.min(x1, feature[1]);
    if (fb > fa) for (let k = 0; k <= 400; k++) xs.push(fa + (fb - fa) * k / 400);
    xs.sort((u, v) => u - v);
    svgEl("polyline", { points: xs.map(x => `${sx(x)},${sy(fcur(x))}`).join(" "), fill: "none", stroke: COL.red, "stroke-width": 3 }, g);
    // nodes
    const nd = run.nodes.filter(x => x >= x0 && x <= x1);
    const r = nd.length > 120 ? 0 : nd.length > 50 ? 2.5 : 4;
    for (const x of nd) {
      svgEl("line", { x1: sx(x), x2: sx(x), y1: bot + 28, y2: bot + 46, stroke: COL.blue, "stroke-width": 1.2 }, g);
      if (r > 0) svgEl("circle", { cx: sx(x), cy: sy(fcur(x)), r, fill: COL.blue }, g);
    }
    // piece bars: green accepted, orange still to check, red stopped at the limit
    for (const p of vis) {
      const c = run.done ? (p.status === "limit" || run.note.startsWith("maximum") ? COL.red : COL.green)
        : st.m === "simpson" ? (p.status === "pending" ? COL.orange : p.status === "ok" ? COL.green : COL.red)
          : (p === worst() ? COL.orange : COL.green);
      const xa = Math.max(L, sx(p.a)), xb = Math.min(R, sx(p.b));
      svgEl("rect", { x: xa + 0.5, y: bot + 50, width: Math.max(1, xb - xa - 1), height: 12, fill: c }, g);
    }
    if (main) svgEl("text", { x: (L + R) / 2, y: H - 4, "text-anchor": "middle" }, svg).textContent =
      st.fn === "pk" ? "time after dose (h)" : st.fn === "chrom" ? "retention time (s)" : "x";
    else svgEl("text", { x: (L + R) / 2, y: H - 4, "text-anchor": "middle", class: "sm" }, svg).textContent =
      "ticks: evaluation points; bars: pieces";
  }
  const worst = (): Piece | undefined => run.pieces.reduce<Piece | undefined>((w, p) => !w || p.err > w.err ? p : w, undefined);

  function draw(): void {
    fnBtns.forEach((b, i) => b.classList.toggle("active", fnKeys[i] === st.fn));
    mBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.m === "simpson")));
    lBtns.forEach((b, i) => b.classList.toggle("active", Ls[i] === st.L));
    const dn = st.fn === "dnorm", ch = st.fn === "chrom";
    lS.input.parentElement!.classList.toggle("ql-hide", !dn);
    gL.classList.toggle("ql-hide", !dn);
    splitCb.parentElement!.classList.toggle("ql-hide", !dn);
    posS.input.parentElement!.classList.toggle("ql-hide", !ch);
    sigS.input.parentElement!.classList.toggle("ql-hide", !ch);
    bRun.textContent = timer !== undefined ? "pause" : "run";

    svgA.replaceChildren(); svgB.replaceChildren();
    const nm = st.fn === "pk" ? "PK model" : st.fn === "chrom" ? `peak at ${st.pos} s, σ = ${st.sig} s` : `dnorm on [0, ${B}]`;
    panel(svgA, WA, A, B, nm, true);
    panel(svgB, WB, zoom[0], zoom[1], `zoom: [${tickFmt(zoom[0])}, ${tickFmt(zoom[1])}]`, false);

    // readout
    ro.replaceChildren();
    const item = (k: string, v: string, cls = "") => {
      const s = htmlEl("span", cls, "", ro); if (k) htmlEl("b", "", k + " ", s); s.append(v);
    };
    const val = total(run), est = errTotal(run), err = Math.abs(val - exact);
    const pend = run.pieces.filter(p => p.status === "pending").length;
    item("status:", run.done ? (run.note.startsWith("maximum") ? "stopped: " : "done: ") + run.note
      : `step ${run.steps}` + (st.m === "simpson" ? `, ${pend} piece(s) still to check` : `, ${run.note ? run.note + "; " : ""}next: bisect the piece with the largest error (orange)`));
    item("tolerance:", st.m === "simpson" ? "1e-6 (absolute)" : `abs = rel = ${sci(GK_TOL, 3)}`);
    if (st.m === "gk" && run.first) item("first pass (21 points):", `${sci(run.first.val, 4)} ± ${sci(run.first.err, 2)}`);
    item("estimate", fix(val, 7));
    item("error estimate", st.m === "simpson" && !run.done ? "–" : sci(est, 2));
    item("true", fix(exact, 7));
    item("actual error", sci(err, 2));
    item("evaluations", String(run.nodes.length));
    item("pieces", String(run.pieces.length));
    if (run.done && err > Math.max(10 * est, 1e-3 * Math.abs(exact)) && err > 1e-9)
      item("", `silent failure: the routine is confident (error estimate ${sci(est, 2)}), but ${(100 * err / exact).toPrecision(3)} % of the area is missing`, "ql-warn");
  }

  recompute();
}

autoMount(NAME, mount);

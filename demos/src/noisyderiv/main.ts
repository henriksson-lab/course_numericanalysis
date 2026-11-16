// Growth rate from noisy plate-reader data, for the finite-differences unit.
// The data are the lecture's 61 OD values (logistic K = 1.2, N0 = 0.05, r = 0.8 /h,
// every 10 min, noise sd 0.005, 3 decimals; hard-coded, so the numbers match R).
// Three estimators of N'(t_i), each with a half-width k (in samples, h = k * dt):
//   central     (y[i+k] - y[i-k]) / (2 k dt)                    k = 1: raw neighbours
//   local line  least-squares slope over y[i-k .. i+k]          (C5, Savitzky-Golay)
//   Richardson  (4 D(h) - D(2h)) / 3 from two central differences
//  (a) the data and the true curve, with the points used at the selected time and the
//      estimated slope there (click a panel to move the selected time);
//  (b) the estimated growth rate against the true N'(t) = r N (1 - N/K);
//  (c) rms error against the half-width h (log-log): measured, the noise model and the
//      truncation error (same estimator on the exact OD values).
// A data switch replaces the measurements by the exact OD values (truncation only).
// Auto-mounts on <div data-demo="noisyderiv">. Optional attributes:
//   data-method="central|local|richardson"  data-k="1"  data-data="noisy|exact"  data-t="4"
import { svgEl, htmlEl, slider, presets, scale, svgPoint, autoMount } from "../lib/ui";

const NAME = "noisyderiv";
// Okabe-Ito, as in make-finite-differences-figs.R
const COL = { cen: "#D55E00", loc: "#0072B2", ric: "#CC79A7", grn: "#009E73", ref: "#737373", ink: "#222", raw: "#b5b5b5" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .45em; }
[data-demo="${NAME}"] .nd-row { display: flex; flex-wrap: wrap; gap: .4em 1.3em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .nd-group { display: flex; gap: .4em; align-items: center; }
[data-demo="${NAME}"] .nd-row .qb-slider { min-width: 300px; flex: 1 1 300px; max-width: 460px; grid-template-columns: auto 1fr 7.5em; }
[data-demo="${NAME}"] .nd-plots { display: flex; gap: 8px; width: 100%; }
[data-demo="${NAME}"] .nd-plots svg { min-width: 0; height: auto; cursor: crosshair; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .nd-read { display: flex; flex-wrap: wrap; gap: .2em 1.4em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .nd-read span b { font-weight: 600; }
[data-demo="${NAME}"] .nd-warn { color: #b00020; font-weight: 600; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 20px; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; }
@media (max-width: 800px) { [data-demo="${NAME}"] .nd-plots { flex-wrap: wrap; } }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- the running example (same numbers as the lecture and exercise C4) ----
const K = 1.2, N0 = 0.05, R = 0.8, SIGMA = 0.005, DT = 1 / 6;
const Y = [
  0.047, 0.058, 0.060, 0.081, 0.084, 0.090, 0.108, 0.123, 0.137, 0.150,
  0.177, 0.192, 0.210, 0.226, 0.269, 0.292, 0.322, 0.359, 0.393, 0.428,
  0.466, 0.504, 0.540, 0.569, 0.622, 0.659, 0.698, 0.730, 0.772, 0.812,
  0.851, 0.876, 0.909, 0.935, 0.955, 0.985, 1.007, 1.029, 1.054, 1.069,
  1.079, 1.092, 1.109, 1.120, 1.123, 1.132, 1.145, 1.154, 1.155, 1.166,
  1.168, 1.167, 1.175, 1.171, 1.187, 1.192, 1.183, 1.181, 1.191, 1.189,
  1.203];
const NP = Y.length;
const T = Y.map((_, i) => i / 6);
const od = (t: number): number => K / (1 + (K / N0 - 1) * Math.exp(-R * t));
const gr = (t: number): number => { const n = od(t); return R * n * (1 - n / K); };
const YEX = T.map(od);
const KMAX = 12;

type Method = "central" | "local" | "richardson";
type Data = "noisy" | "exact";
const METHODS: Method[] = ["central", "local", "richardson"];
const MCOL: Record<Method, string> = { central: COL.cen, local: COL.loc, richardson: COL.ric };
const MLAB: Record<Method, string> = { central: "central", local: "local line", richardson: "Richardson" };

/** Points lost at each end. */
const reach = (m: Method, k: number): number => m === "richardson" ? 2 * k : k;

/** Estimate N'(t_i) from data y with half-width k; NaN where the stencil does not fit. */
function estimate(y: number[], m: Method, k: number, i: number): number {
  if (i - reach(m, k) < 0 || i + reach(m, k) >= NP) return NaN;
  const cen = (kk: number) => (y[i + kk] - y[i - kk]) / (2 * kk * DT);
  if (m === "central") return cen(k);
  if (m === "richardson") return (4 * cen(k) - cen(2 * k)) / 3;
  let num = 0, den = 0;
  for (let j = -k; j <= k; j++) { num += j * y[i + j]; den += j * j; }
  return num / (DT * den);
}
/** Weights of the estimator on y[i+j], j = -reach..reach (for drawing and the noise sd). */
function weights(m: Method, k: number): Map<number, number> {
  const w = new Map<number, number>();
  const add = (j: number, v: number) => w.set(j, (w.get(j) ?? 0) + v);
  if (m === "central") { add(k, 1 / (2 * k * DT)); add(-k, -1 / (2 * k * DT)); }
  else if (m === "richardson") {
    add(k, 4 / 3 / (2 * k * DT)); add(-k, -4 / 3 / (2 * k * DT));
    add(2 * k, -1 / 3 / (4 * k * DT)); add(-2 * k, 1 / 3 / (4 * k * DT));
  } else {
    let den = 0;
    for (let j = -k; j <= k; j++) den += j * j;
    for (let j = -k; j <= k; j++) add(j, j / (DT * den));
  }
  return w;
}
const noiseSd = (m: Method, k: number): number => {
  let s = 0;
  weights(m, k).forEach(v => { s += v * v; });
  return SIGMA * Math.sqrt(s);
};
interface Stats { rms: number; n: number; max: number; tmax: number; i0: number; i1: number }
function stats(y: number[], m: Method, k: number): Stats {
  const r = reach(m, k);
  let s = 0, n = 0, max = -Infinity, tmax = NaN;
  for (let i = r; i < NP - r; i++) {
    const d = estimate(y, m, k, i), e = d - gr(T[i]);
    s += e * e; n++;
    if (d > max) { max = d; tmax = T[i]; }
  }
  return { rms: Math.sqrt(s / n), n, max, tmax, i0: r, i1: NP - 1 - r };
}
// precomputed rms curves for panel (c)
const KS = Array.from({ length: KMAX }, (_, i) => i + 1);
const RMS: Record<Data, Record<Method, number[]>> = {
  noisy: { central: [], local: [], richardson: [] },
  exact: { central: [], local: [], richardson: [] },
};
for (const m of METHODS) for (const k of KS) {
  RMS.noisy[m].push(stats(Y, m, k).rms);
  RMS.exact[m].push(stats(YEX, m, k).rms);
}

const sup: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const pow10 = (k: number): string => k === 0 ? "1" : "10" + String(k).split("").map(c => sup[c]).join("");
function sci(v: number, d = 2): string {
  if (!Number.isFinite(v)) return "–";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e-3 && a < 1e3) return v.toPrecision(d);
  const [mm, e] = v.toExponential(d - 1).split("e");
  return `${mm}e${Number(e)}`;
}
const fmtT = (t: number): string => t.toFixed(2) + " h";

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const start = {
    m: (METHODS as string[]).includes(ds.method ?? "") ? ds.method as Method : "central" as Method,
    k: Math.max(1, Math.min(KMAX, Math.round(Number(ds.k) || 1))),
    data: (ds.data === "exact" ? "exact" : "noisy") as Data,
    i: Math.max(0, Math.min(NP - 1, Math.round(6 * (Number.isFinite(Number(ds.t)) && ds.t !== undefined ? Number(ds.t) : 4)))),
  };
  const st = { ...start };

  // ---- controls ----
  const row1 = htmlEl("div", "nd-row", "", host);
  const mg = htmlEl("div", "nd-group", "", row1);
  htmlEl("span", "", "estimator:", mg);
  const mBtns = presets(mg, METHODS.map(m => ({ label: MLAB[m], apply: () => { st.m = m; draw(); } })));
  const dg = htmlEl("div", "nd-group", "", row1);
  htmlEl("span", "", "data:", dg);
  const dBtns = presets(dg, [
    { label: "plate reader (σ = 0.005)", apply: () => { st.data = "noisy"; draw(); } },
    { label: "exact OD", apply: () => { st.data = "exact"; draw(); } },
  ]);
  const row2 = htmlEl("div", "nd-row", "", host);
  const kS = slider(row2, "half-width k", 1, KMAX, 1, st.k, v => { st.k = v; draw(); },
    v => `${v} (${v * 10} min)`);
  const pg = htmlEl("div", "nd-group", "", row2);
  const pre: { label: string; m: Method; k: number }[] = [
    { label: "neighbours", m: "central", k: 1 },
    { label: "central k = 4", m: "central", k: 4 },
    { label: "local line k = 4", m: "local", k: 4 },
    { label: "Richardson k = 1", m: "richardson", k: 1 },
  ];
  const pBtns = presets(pg, [
    ...pre.map(p => ({ label: p.label, apply: () => { st.m = p.m; st.k = p.k; st.data = "noisy"; kS.set(st.k); draw(); } })),
    { label: "reset", apply: () => { Object.assign(st, start); kS.set(st.k); draw(); } },
  ]);

  const plots = htmlEl("div", "nd-plots", "", host);
  const H = 350, WA = 360, WB = 400, WC = 390, WT = WA + WB + WC;
  const mk = (w: number, lab: string) => {
    const s = svgEl("svg", { viewBox: `0 0 ${w} ${H}`, role: "img", "aria-label": lab }, plots);
    s.style.width = `calc(${(100 * w / WT).toFixed(2)}% - 6px)`;
    return s;
  };
  const svgA = mk(WA, "OD data with the stencil at the selected time");
  const svgB = mk(WB, "Estimated and true growth rate");
  const svgC = mk(WC, "rms error against half-width, log-log");
  const ro = htmlEl("div", "qb-readout nd-read", "", host);

  const top = 36, bot = H - 50, L = 60;
  const ax = scale(0, 10, L, WA - 10), ay = scale(0, 1.3, bot, top);
  const bx = scale(0, 10, L + 8, WB - 10), GY0 = -0.06, GY1 = 0.32, byy = scale(GY0, GY1, bot, top);
  const CX0 = Math.log10(1 / 6) - 0.05, CX1 = Math.log10(2) + 0.05, CY0 = -7, CY1 = -1;
  const cx = scale(CX0, CX1, L + 12, WC - 12), cy = scale(CY0, CY1, bot, top);

  // click (or drag) in (a) or (b) selects the time point
  const pick = (svg: SVGSVGElement, xs: ReturnType<typeof scale>) => {
    const go = (ev: PointerEvent) => {
      const [x] = svgPoint(svg, ev);
      const i = Math.max(0, Math.min(NP - 1, Math.round(6 * xs.invert(x))));
      if (i !== st.i) { st.i = i; draw(); }
    };
    svg.addEventListener("pointerdown", ev => {
      go(ev);
      const mv = (e: PointerEvent) => go(e);
      const up = () => { svg.removeEventListener("pointermove", mv); window.removeEventListener("pointerup", up); };
      svg.addEventListener("pointermove", mv);
      window.addEventListener("pointerup", up);
    });
  };
  pick(svgA, ax); pick(svgB, bx);

  const txt = (svg: SVGSVGElement, x: number, y: number, s: string, a: Record<string, string | number> = {}, fill = "") => {
    const t = svgEl("text", { x, y, ...a }, svg);
    if (fill) t.style.fill = fill;
    t.textContent = s;
    return t;
  };
  const ylab = (svg: SVGSVGElement, s: string) =>
    txt(svg, 15, (top + bot) / 2, s, { "text-anchor": "middle", transform: `rotate(-90 15 ${(top + bot) / 2})` });

  function draw(): void {
    const { m, k, data, i } = st;
    const y = data === "noisy" ? Y : YEX;
    const col = MCOL[m];
    mBtns.forEach((b, j) => b.classList.toggle("active", METHODS[j] === m));
    dBtns.forEach((b, j) => b.classList.toggle("active", (j === 0) === (data === "noisy")));
    pBtns.forEach((b, j) => b.classList.toggle("active", j < pre.length && pre[j].m === m && pre[j].k === k && data === "noisy"));
    const S = stats(y, m, k);
    const di = estimate(y, m, k, i);

    // ---------- (a) data and stencil ----------
    svgA.replaceChildren();
    txt(svgA, 8, 22, "(a) OD data and stencil", { class: "ttl" });
    for (const v of [0, 0.4, 0.8, 1.2]) {
      svgEl("line", { x1: L, x2: WA - 10, y1: ay(v), y2: ay(v), class: v === 0 ? "axis" : "grid" }, svgA);
      txt(svgA, L - 6, ay(v) + 6, v.toFixed(1), { "text-anchor": "end", class: "sm" });
    }
    for (let t = 0; t <= 10; t += 2) {
      svgEl("line", { x1: ax(t), x2: ax(t), y1: bot, y2: bot + 5, class: "axis" }, svgA);
      txt(svgA, ax(t), bot + 22, String(t), { "text-anchor": "middle", class: "sm" });
    }
    txt(svgA, (L + WA) / 2, H - 6, "time t (h)", { "text-anchor": "middle" });
    ylab(svgA, "OD600");
    const r = reach(m, k);
    const lo = Math.max(0, i - r), hi = Math.min(NP - 1, i + r);
    svgEl("rect", { x: ax(T[lo]) - 3, y: top, width: ax(T[hi]) - ax(T[lo]) + 6, height: bot - top, fill: col, "fill-opacity": 0.1 }, svgA);
    const curve = Array.from({ length: 201 }, (_, j) => j / 20).map(t => `${ax(t)},${ay(od(t))}`).join(" ");
    svgEl("polyline", { points: curve, fill: "none", stroke: COL.ref, "stroke-width": 2 }, svgA);
    const w = weights(m, k);
    for (let j = 0; j < NP; j++) {
      const used = Number.isFinite(di) && w.has(j - i) && Math.abs(w.get(j - i)!) > 0;
      svgEl("circle", { cx: ax(T[j]), cy: ay(y[j]), r: used ? 5.5 : 3, fill: used ? col : COL.ink, "fill-opacity": used ? 1 : 0.75, stroke: used ? COL.ink : "none", "stroke-width": 1.2 }, svgA);
    }
    // estimated slope through the mean of the window, and the true tangent
    const seg = (slope: number, y0: number, half: number, stroke: string, wdt: number, dash = "") => {
      const ta = T[i] - half, tb = T[i] + half;
      svgEl("line", { x1: ax(ta), y1: ay(y0 - slope * half), x2: ax(tb), y2: ay(y0 + slope * half), stroke, "stroke-width": wdt, "stroke-dasharray": dash }, svgA);
    };
    const half = Math.max(1.0, (r + 1) * DT);
    if (Number.isFinite(di)) {
      let ym = 0, c = 0;
      for (let j = lo; j <= hi; j++) { ym += y[j]; c++; }
      seg(di, ym / c, half, col, 4);
    }
    seg(gr(T[i]), od(T[i]), half, COL.grn, 2.5, "7 5");
    svgEl("line", { x1: ax(T[i]), x2: ax(T[i]), y1: top, y2: bot, stroke: COL.ink, "stroke-width": 1, "stroke-dasharray": "2 3" }, svgA);
    txt(svgA, WA - 14, ay(0.22), "– – true tangent", { "font-weight": 600, "text-anchor": "end" }, COL.grn);
    txt(svgA, WA - 14, ay(0.22) + 24, `—  ${MLAB[m]} slope`, { "font-weight": 600, "text-anchor": "end" }, col);

    // ---------- (b) growth rate ----------
    svgB.replaceChildren();
    txt(svgB, 8, 22, "(b) growth rate N′(t)", { class: "ttl" });
    for (const v of [0, 0.1, 0.2, 0.3]) {
      svgEl("line", { x1: L + 8, x2: WB - 10, y1: byy(v), y2: byy(v), class: v === 0 ? "axis" : "grid" }, svgB);
      txt(svgB, L + 2, byy(v) + 6, v.toFixed(1), { "text-anchor": "end", class: "sm" });
    }
    for (let t = 0; t <= 10; t += 2) {
      svgEl("line", { x1: bx(t), x2: bx(t), y1: bot, y2: bot + 5, class: "axis" }, svgB);
      txt(svgB, bx(t), bot + 22, String(t), { "text-anchor": "middle", class: "sm" });
    }
    txt(svgB, (L + WB) / 2, H - 6, "time t (h)", { "text-anchor": "middle" });
    ylab(svgB, "OD/h");
    const cyB = (v: number) => byy(Math.max(GY0, Math.min(GY1, v)));
    // raw neighbours, faint
    for (let j = 1; j < NP - 1; j++)
      svgEl("circle", { cx: bx(T[j]), cy: cyB(estimate(y, "central", 1, j)), r: 3, fill: COL.raw }, svgB);
    const tru = Array.from({ length: 201 }, (_, j) => j / 20).map(t => `${bx(t)},${byy(gr(t))}`).join(" ");
    svgEl("polyline", { points: tru, fill: "none", stroke: COL.ink, "stroke-width": 4 }, svgB);
    const est: string[] = [];
    for (let j = S.i0; j <= S.i1; j++) est.push(`${bx(T[j])},${cyB(estimate(y, m, k, j))}`);
    svgEl("polyline", { points: est.join(" "), fill: "none", stroke: col, "stroke-width": 2.5 }, svgB);
    for (let j = S.i0; j <= S.i1; j++)
      svgEl("circle", { cx: bx(T[j]), cy: cyB(estimate(y, m, k, j)), r: 3.2, fill: col }, svgB);
    svgEl("line", { x1: bx(T[i]), x2: bx(T[i]), y1: top, y2: bot, stroke: COL.ink, "stroke-width": 1, "stroke-dasharray": "2 3" }, svgB);
    if (Number.isFinite(di)) svgEl("circle", { cx: bx(T[i]), cy: cyB(di), r: 7.5, fill: col, stroke: COL.ink, "stroke-width": 2 }, svgB);
    const lgB = (s: string, yy: number, c: string, dot = false) => {
      const x = bx(6.0);
      if (dot) svgEl("circle", { cx: x + 10, cy: yy - 6, r: 4, fill: c }, svgB);
      else svgEl("line", { x1: x, x2: x + 22, y1: yy - 6, y2: yy - 6, stroke: c, "stroke-width": 4 }, svgB);
      txt(svgB, x + 28, yy, s, {}, c === COL.raw ? "#777" : c);
    };
    lgB("true", top + 20, COL.ink);
    lgB(MLAB[m], top + 42, col);
    lgB("neighbours", top + 64, COL.raw, true);

    // ---------- (c) rms error vs half-width ----------
    svgC.replaceChildren();
    txt(svgC, 8, 22, "(c) rms error against half-width", { class: "ttl" });
    for (let e = CY0; e <= CY1; e++) {
      svgEl("line", { x1: L + 12, x2: WC - 12, y1: cy(e), y2: cy(e), class: "grid" }, svgC);
      txt(svgC, L + 6, cy(e) + 6, pow10(e), { "text-anchor": "end", class: "sm" });
    }
    svgEl("line", { x1: L + 12, x2: WC - 12, y1: bot, y2: bot, class: "axis" }, svgC);
    for (const kk of [1, 2, 4, 8, 12]) {
      const xx = cx(Math.log10(kk * DT));
      svgEl("line", { x1: xx, x2: xx, y1: bot, y2: bot + 5, class: "axis" }, svgC);
      txt(svgC, xx, bot + 22, String(kk * 10), { "text-anchor": "middle", class: "sm" });
    }
    txt(svgC, (L + WC) / 2, H - 6, "half-width h = k·10 min (min)", { "text-anchor": "middle" });
    ylab(svgC, "rms error (OD/h)");
    const cl = (v: number) => cy(Math.max(CY0 - 0.3, Math.min(CY1 + 0.3, Math.log10(v))));
    const xk = (kk: number) => cx(Math.log10(kk * DT));
    const poly = (vals: number[], stroke: string, wdt: number, dash = "", op = 1) =>
      svgEl("polyline", { points: KS.map((kk, j) => `${xk(kk)},${cl(vals[j])}`).join(" "), fill: "none", stroke, "stroke-width": wdt, "stroke-dasharray": dash, "stroke-opacity": op }, svgC);
    // models for the current estimator
    if (data === "noisy") poly(KS.map(kk => noiseSd(m, kk)), COL.ref, 2.2, "9 6");
    poly(RMS.exact[m], COL.ref, 2.4, "3 4");
    // measured: other estimators faint, the current one strong
    for (const mm of METHODS) if (mm !== m) poly(RMS[data][mm], MCOL[mm], 2, "", 0.3);
    poly(RMS[data][m], col, 3);
    KS.forEach((kk, j) => svgEl("circle", { cx: xk(kk), cy: cl(RMS[data][m][j]), r: 3.5, fill: col }, svgC));
    svgEl("circle", { cx: xk(k), cy: cl(RMS[data][m][k - 1]), r: 9, fill: col, stroke: COL.ink, "stroke-width": 2.5 }, svgC);
    // legend, lower right (empty in every mode)
    const lgC = (label: string, yy: number, stroke: string, dash: string, wdt: number) => {
      const x = cx(Math.log10(2.5 * DT));
      svgEl("line", { x1: x, x2: x + 34, y1: yy - 6, y2: yy - 6, stroke, "stroke-width": wdt, "stroke-dasharray": dash }, svgC);
      txt(svgC, x + 42, yy, label, {}, stroke === COL.ref ? "#444" : stroke);
    };
    let ly = cy(-5.2);
    lgC(data === "noisy" ? "measured" : "exact OD", ly, col, "", 3);
    if (data === "noisy") lgC("noise model", ly += 24, COL.ref, "9 6", 2.2);
    lgC("truncation", ly += 24, COL.ref, "3 4", 2.4);

    // ---------- readout ----------
    ro.replaceChildren();
    const item = (kk: string, v: string, cls = "") => {
      const s = htmlEl("span", cls, "", ro);
      htmlEl("b", "", kk + " ", s); s.append(v);
    };
    const h = k * 10;
    item(`${MLAB[m]}, h =`, m === "richardson" ? `${h} and ${2 * h} min` : `${h} min`);
    item("rms error", `${sci(S.rms, 3)} OD/h`);
    if (data === "noisy") item("noise sd", sci(noiseSd(m, k), 2));
    item("truncation rms", sci(RMS.exact[m][k - 1], 2));
    if (m !== "central") item("central, same k:", sci(RMS[data].central[k - 1], 3));
    item("max", `${S.max.toFixed(3)} at ${fmtT(S.tmax)}`);
    item("estimates for", `t = ${fmtT(T[S.i0])} … ${fmtT(T[S.i1])} (${r} lost each end)`);
    if (Number.isFinite(di)) item(`at t = ${fmtT(T[i])}:`, `${di.toFixed(4)} (true ${gr(T[i]).toFixed(4)})`);
    else item(`at t = ${fmtT(T[i])}:`, "no estimate, stencil runs past the data", "nd-warn");
  }
  draw();
}

autoMount(NAME, mount);

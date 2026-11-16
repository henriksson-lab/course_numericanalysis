// Step size explorer for the convergence-cost unit. The growth rate of the lecture's
// logistic OD curve, N(t) = K / (1 + (K/N0 - 1) e^{-rt}), K = 1.2, N0 = 0.05, r = 0.8,
// at t = 3 h, by a forward or a central difference.
//  (a) the OD curve with the tangent (true slope) and the secant through t, t+h
//      (or t-h, t+h) for the current h;
//  (b) log-log error against h for h = 1 ... 1e-16: the actual errors computed here in
//      JS doubles, the model lines C h^p (truncation) and eps |N| / h (rounding), and a
//      marker at the current h.
// A second tab ("Cost and memory", hidden with data-cost="0") shows operation counts
// 2n, 2n^2, 2n^3 and n^2 m as times at a fixed flop rate, and the memory of a dense and a
// sparse genes x cells matrix and of the gene x gene correlation matrix.
// Auto-mounts on <div data-demo="stepsize">. Optional attributes:
//   data-method="forward|central"  data-logh="-1"  (log10 of the starting h)
//   data-cost="1|0"  data-tab="step|cost"  data-flops="1e9"
import { svgEl, htmlEl, slider, presets, scale, autoMount } from "../lib/ui";

const NAME = "stepsize";
// Okabe-Ito, as in make-convergence-cost-figs.R
const COL = { fwd: "#0072B2", cen: "#D55E00", grn: "#009E73", pur: "#CC79A7", ref: "#737373", ink: "#222" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .45em; }
[data-demo="${NAME}"] .ss-row { display: flex; flex-wrap: wrap; gap: .4em 1.3em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .ss-group { display: flex; gap: .4em; align-items: center; }
[data-demo="${NAME}"] .ss-row .qb-slider { min-width: 360px; flex: 1 1 360px; max-width: 520px; grid-template-columns: auto 1fr 5.2em; }
[data-demo="${NAME}"] .ss-plots { display: flex; flex-wrap: wrap; gap: 10px; width: 100%; }
[data-demo="${NAME}"] .ss-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .ss-read { display: flex; flex-wrap: wrap; gap: .2em 1.4em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .ss-read span b { font-weight: 600; }
[data-demo="${NAME}"] .ss-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .ss-tabs .qb-btn { font-size: .95em; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 20px; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 18px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- the running example (same numbers as the slides) ----
const K = 1.2, N0 = 0.05, R = 0.8, T0 = 3;
const EPS = 2.220446049250313e-16;
const od = (t: number): number => K / (1 + (K / N0 - 1) * Math.exp(-R * t));
const N = od(T0);
const D1 = R * N * (1 - N / K);                                // N'  (logistic ODE)
const D2 = R * (1 - 2 * N / K) * D1;                           // N''
const D3 = -2 * R / K * D1 * D1 + R * (1 - 2 * N / K) * D2;    // N'''

type Method = "forward" | "central";
const est = (m: Method, h: number): number =>
  m === "forward" ? (od(T0 + h) - od(T0)) / h : (od(T0 + h) - od(T0 - h)) / (2 * h);
const err = (m: Method, h: number): number => Math.abs(est(m, h) - D1);
const MODEL = {
  forward: { C: Math.abs(D2) / 2, p: 1 },
  central: { C: Math.abs(D3) / 6, p: 2 },
};
const hOpt = (m: Method): number => m === "forward"
  ? Math.sqrt(EPS * N / MODEL.forward.C)
  : Math.cbrt(EPS * N / (2 * MODEL.central.C));

// error curves, precomputed on a grid of log10 h
const LH0 = -16, LH1 = 0, GRID = 0.1;
const grid = Array.from({ length: Math.round((LH1 - LH0) / GRID) + 1 }, (_, i) => LH0 + i * GRID);
const curves: Record<Method, number[]> = {
  forward: grid.map(l => err("forward", 10 ** l)),
  central: grid.map(l => err("central", 10 ** l)),
};

// ---- formatting ----
const sup: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const pow10 = (k: number): string => k === 0 ? "1" : "10" + String(k).split("").map(c => sup[c]).join("");
function sci(v: number, d = 2): string {
  if (v === 0) return "0";
  const s = v.toExponential(d - 1);
  const [m, e] = s.split("e");
  const ei = Number(e);
  if (ei >= -2 && ei <= 0) return v.toPrecision(d);
  return `${m}e${ei}`;
}
const hFmt = (lh: number): string => { const h = 10 ** lh; return Math.abs(lh - Math.round(lh)) < 1e-9 ? `1e${Math.round(lh)}`.replace("1e0", "1") : sci(h, 2); };
function fmtTime(s: number): string {
  const f = (v: number, u: string) => `${v >= 100 ? v.toFixed(0) : v >= 10 ? v.toFixed(1) : v.toFixed(2)} ${u}`;
  if (s < 1e-6) return f(s * 1e9, "ns");
  if (s < 1e-3) return f(s * 1e6, "µs");
  if (s < 1) return f(s * 1e3, "ms");
  if (s < 60) return f(s, "s");
  if (s < 3600) return f(s / 60, "min");
  if (s < 86400) return f(s / 3600, "h");
  if (s < 365 * 86400) return f(s / 86400, "days");
  return f(s / (365 * 86400), "years");
}
const fmtGB = (b: number): string => { const g = b / 1e9; return g >= 100 ? g.toFixed(0) + " GB" : g >= 1 ? g.toFixed(1) + " GB" : g >= 0.001 ? (g * 1000).toPrecision(2) + " MB" : (g * 1e6).toPrecision(2) + " kB"; };
const fmtInt = (n: number): string => Math.round(n).toString().replace(/\B(?=(\d{3})+(?!\d))/g, " ");

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const st = {
    m: (ds.method === "central" ? "central" : "forward") as Method,
    lh: Number.isFinite(Number(ds.logh)) && ds.logh !== undefined ? Math.min(LH1, Math.max(LH0, Number(ds.logh))) : -1,
  };
  const startM = st.m, startLh = st.lh;
  const uid = ++mountCount;
  const withCost = ds.cost !== "0";
  const flops = Number(ds.flops) > 0 ? Number(ds.flops) : 1e9;

  // ---- tabs ----
  const stepBox = htmlEl("div", "", "");
  const costBox = htmlEl("div", "", "");
  for (const b of [stepBox, costBox]) Object.assign(b.style, { display: "flex", flexDirection: "column", gap: ".45em", width: "100%" });
  if (withCost) {
    const tabs = htmlEl("div", "ss-row ss-tabs", "", host);
    const tb = presets(tabs, [
      { label: "Step size: accuracy", apply: () => show("step") },
      { label: "Cost and memory", apply: () => show("cost") },
    ]);
    const show = (t: string) => {
      stepBox.style.display = t === "step" ? "flex" : "none";
      costBox.style.display = t === "cost" ? "flex" : "none";
      tb.forEach((b, i) => b.classList.toggle("active", (i === 0) === (t === "step")));
    };
    host.append(stepBox, costBox);
    show(ds.tab === "cost" ? "cost" : "step");
    buildCost(costBox, flops);
  } else host.append(stepBox);

  // ---- step-size controls ----
  const ctr = htmlEl("div", "ss-row", "", stepBox);
  const mg = htmlEl("div", "ss-group", "", ctr);
  htmlEl("span", "", "formula:", mg);
  const mBtns = presets(mg, [
    { label: "forward", apply: () => { st.m = "forward"; draw(); } },
    { label: "central", apply: () => { st.m = "central"; draw(); } },
  ]);
  const hS = slider(ctr, "step h", LH0, LH1, 0.05, st.lh, v => { st.lh = v; draw(); }, hFmt);
  const ctr2 = htmlEl("div", "ss-row", "", stepBox);
  const hg = htmlEl("div", "ss-group", "", ctr2);
  htmlEl("span", "", "h =", hg);
  const hVals = [0, -1, -4, -5, -8, -12, -15];
  const hBtns = presets(hg, hVals.map(l => ({ label: pow10(l), apply: () => { st.lh = l; hS.set(l); draw(); } })));
  const optBtn = htmlEl("button", "qb-btn", "best h*", hg);
  optBtn.type = "button";
  optBtn.addEventListener("click", () => { st.lh = Math.round(Math.log10(hOpt(st.m)) / 0.05) * 0.05; hS.set(st.lh); draw(); optBtn.blur(); });
  const rst = htmlEl("button", "qb-btn", "reset", ctr2);
  rst.type = "button";
  rst.addEventListener("click", () => { st.m = startM; st.lh = startLh; hS.set(st.lh); draw(); rst.blur(); });

  const plots = htmlEl("div", "ss-plots", "", stepBox);
  const H = 380, WA = 450, WB = 650, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "OD curve with tangent and secant" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Error against step size, log-log" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro = htmlEl("div", "qb-readout ss-read", "", stepBox);

  // panel (a) scales
  const top = 40, bot = H - 52;
  const ax = scale(0, 10, 62, WA - 12), ay = scale(0, 1.25, bot, top);
  // panel (b) scales: log10 h on x, log10 err on y
  const EY0 = -17, EY1 = 0;
  const bx = scale(LH0, LH1, 72, WB - 14), by = scale(EY0, EY1, bot, top);
  const clampY = (e: number) => Math.max(EY0, Math.min(EY1 + 0.6, Math.log10(Math.max(e, 1e-30))));

  function draw(): void {
    const m = st.m, h = 10 ** st.lh;
    mBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (m === "forward")));
    hBtns.forEach((b, i) => b.classList.toggle("active", Math.abs(hVals[i] - st.lh) < 1e-9));
    const D = est(m, h), e = Math.abs(D - D1);
    const col = m === "forward" ? COL.fwd : COL.cen;

    // ---------- (a) ----------
    svgA.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent = "(a) OD curve, tangent and secant";
    for (let v = 0; v <= 1.2 + 1e-9; v += 0.4) {
      svgEl("line", { x1: 62, x2: WA - 12, y1: ay(v), y2: ay(v), class: v === 0 ? "axis" : "grid" }, svgA);
      svgEl("text", { x: 56, y: ay(v) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = v.toFixed(1);
    }
    for (let t = 0; t <= 10; t += 2) {
      svgEl("line", { x1: ax(t), x2: ax(t), y1: bot, y2: bot + 5, class: "axis" }, svgA);
      svgEl("text", { x: ax(t), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgA).textContent = String(t);
    }
    svgEl("text", { x: (62 + WA) / 2, y: H - 8, "text-anchor": "middle" }, svgA).textContent = "time t (h)";
    const yl = svgEl("text", { x: 16, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(top + bot) / 2})` }, svgA);
    yl.textContent = "OD600";
    const pts = Array.from({ length: 201 }, (_, i) => i * 10 / 200).map(t => `${ax(t)},${ay(od(t))}`).join(" ");
    svgEl("polyline", { points: pts, fill: "none", stroke: COL.ink, "stroke-width": 3 }, svgA);
    const clip = (t: number, y: number): [number, number] => [ax(t), ay(Math.max(-0.05, Math.min(1.3, y)))];
    const line = (slope: number, t0: number, y0: number, stroke: string, dash: string, w: number) => {
      const ta = 1.2, tb = 5.6;
      const [x1, y1] = clip(ta, y0 + slope * (ta - t0)), [x2, y2] = clip(tb, y0 + slope * (tb - t0));
      svgEl("line", { x1, y1, x2, y2, stroke, "stroke-width": w, "stroke-dasharray": dash }, svgA);
    };
    line(D1, T0, N, COL.grn, "", 3);
    const tl = m === "forward" ? T0 : T0 - h, tr = T0 + h;
    line(D, tl, od(tl), col, "8 6", 3);
    // h bracket
    const yb = od(tl);
    if (ax(tr) - ax(tl) > 6) {
      svgEl("line", { x1: ax(tl), x2: ax(tr), y1: ay(yb), y2: ay(yb), stroke: "#666", "stroke-width": 1.5 }, svgA);
      svgEl("line", { x1: ax(tr), x2: ax(tr), y1: ay(yb), y2: ay(od(tr)), stroke: "#666", "stroke-width": 1.5 }, svgA);
      svgEl("text", { x: (ax(tl) + ax(tr)) / 2, y: ay(yb) + 22, "text-anchor": "middle", class: "sm" }, svgA).textContent =
        m === "forward" ? "h" : "2h";
    }
    for (const t of [tl, tr]) svgEl("circle", { cx: ax(t), cy: ay(od(t)), r: 6, fill: col, stroke: "#fff", "stroke-width": 1.5 }, svgA);
    svgEl("circle", { cx: ax(T0), cy: ay(N), r: 4, fill: COL.grn }, svgA);
    const lg = (txt: string, y: number, c: string) => {
      const t = svgEl("text", { x: 70, y, "font-weight": 600 }, svgA);
      t.style.fill = c; t.textContent = txt;
    };
    lg("— tangent, slope N′(3)", top + 24, COL.grn);
    lg(`- - ${m} difference`, top + 48, col);
    if (h < 0.05) {
      svgEl("text", { x: WA - 16, y: ay(0.3), class: "sm", "text-anchor": "end" }, svgA).textContent = "h is too small to see;";
      svgEl("text", { x: WA - 16, y: ay(0.3) + 20, class: "sm", "text-anchor": "end" }, svgA).textContent = "dashed slope = estimate";
    }

    // ---------- (b) ----------
    svgB.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgB).textContent = "(b) error against step h (log–log)";
    for (let k = EY0 + 1; k <= EY1; k += 2) {
      svgEl("line", { x1: 72, x2: WB - 14, y1: by(k), y2: by(k), class: "grid" }, svgB);
      svgEl("text", { x: 66, y: by(k) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = pow10(k);
    }
    svgEl("line", { x1: 72, x2: WB - 14, y1: bot, y2: bot, class: "axis" }, svgB);
    for (let k = LH0; k <= LH1; k += 2) {
      svgEl("line", { x1: bx(k), x2: bx(k), y1: bot, y2: bot + 5, class: "axis" }, svgB);
      svgEl("text", { x: bx(k), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgB).textContent = pow10(k);
    }
    svgEl("text", { x: (72 + WB) / 2, y: H - 8, "text-anchor": "middle" }, svgB).textContent = "step h";
    const yl2 = svgEl("text", { x: 16, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(top + bot) / 2})` }, svgB);
    yl2.textContent = "absolute error";
    const g = svgEl("g", {}, svgB);
    const cp = svgEl("clipPath", { id: "ss-clip-" + uid }, g);
    svgEl("rect", { x: 72, y: top - 4, width: WB - 86, height: bot - top + 4 }, cp);
    const pg = svgEl("g", { "clip-path": `url(#ss-clip-${uid})` }, g);
    // model lines for the current method
    const md = MODEL[m];
    const lineLog = (f: (l: number) => number, stroke: string, dash: string) => {
      const p = grid.map(l => `${bx(l)},${by(clampY(f(l)))}`).join(" ");
      svgEl("polyline", { points: p, fill: "none", stroke, "stroke-width": 2.2, "stroke-dasharray": dash }, pg);
    };
    lineLog(l => md.C * 10 ** (md.p * l), COL.ref, "9 6");
    lineLog(l => EPS * N / 10 ** l, COL.ref, "3 5");
    // h* line
    const ho = Math.log10(hOpt(m));
    svgEl("line", { x1: bx(ho), x2: bx(ho), y1: top + 30, y2: bot, stroke: COL.grn, "stroke-width": 2 }, pg);
    const hoT = svgEl("text", { x: bx(ho) + (m === "forward" ? 6 : -6), y: top + 46, class: "sm", "text-anchor": m === "forward" ? "start" : "end" }, pg);
    hoT.style.fill = COL.grn; hoT.textContent = (m === "forward" ? "h* ≈ √ε" : "h* ≈ ε^(1/3)") + ` ≈ ${sci(10 ** ho, 1)}`;
    // actual curves: the other method faint, the current one strong
    const other: Method = m === "forward" ? "central" : "forward";
    for (const mm of [other, m]) {
      const cc = mm === "forward" ? COL.fwd : COL.cen, strong = mm === m;
      const p = grid.map((l, i) => `${bx(l)},${by(clampY(curves[mm][i]))}`).join(" ");
      svgEl("polyline", { points: p, fill: "none", stroke: cc, "stroke-width": strong ? 2.2 : 1.5, "stroke-opacity": strong ? 1 : 0.35 }, pg);
      if (strong) grid.forEach((l, i) => {
        if (i % 2) return;
        svgEl("circle", { cx: bx(l), cy: by(clampY(curves[mm][i])), r: 2.8, fill: cc }, pg);
      });
    }
    // labels of model lines
    const lab = (x: number, y: number, txt: string, anchor = "start") => {
      const t = svgEl("text", { x, y, "text-anchor": anchor, class: "sm" }, svgB);
      t.style.fill = "#444"; t.textContent = txt;
    };
    lab(WB - 18, by(m === "forward" ? -6 : -9.5),
      m === "forward" ? `truncation ≈ ${sci(md.C, 2)}·h` : `truncation ≈ ${sci(md.C, 2)}·h²`, "end");
    lab(WB - 18, by(m === "forward" ? -7 : -10.5), `(slope ${md.p})`, "end");
    lab(bx(-15.7), by(-7), "rounding ≈ ε|N|/h", "start");
    lab(bx(-15.7), by(-8), "(slope −1)", "start");
    // current h
    const ey = clampY(e);
    svgEl("line", { x1: bx(st.lh), x2: bx(st.lh), y1: top, y2: bot, stroke: COL.ink, "stroke-width": 1, "stroke-dasharray": "2 3" }, svgB);
    svgEl("circle", { cx: bx(st.lh), cy: by(ey), r: 9, fill: col, stroke: COL.ink, "stroke-width": 2.5 }, svgB);
    // legend
    const lgB = (txt: string, y: number, c: string, op = 1) => {
      const x = bx(-3);
      svgEl("line", { x1: x, x2: x + 22, y1: y - 6, y2: y - 6, stroke: c, "stroke-width": 3, "stroke-opacity": op }, svgB);
      const t = svgEl("text", { x: x + 28, y }, svgB);
      t.style.fill = c; t.style.fillOpacity = String(op); t.textContent = txt;
    };
    lgB("forward", by(-10.3), COL.fwd, m === "forward" ? 1 : 0.45);
    lgB("central", by(-11.6), COL.cen, m === "central" ? 1 : 0.45);

    // ---------- readout ----------
    const eHalf = err(m, h / 2);
    const slope = e > 0 && eHalf > 0 ? Math.log2(e / eHalf) : NaN;
    const same = m === "forward" ? od(T0 + h) === od(T0) : od(T0 + h) === od(T0 - h);
    ro.replaceChildren();
    const item = (k: string, v: string, cls = "") => {
      const s = htmlEl("span", cls, "", ro);
      htmlEl("b", "", k + " ", s); s.append(v);
    };
    item("h =", sci(h, 2));
    item("estimate", D.toPrecision(10));
    item("true N′(3) =", D1.toPrecision(10));
    item("abs. error", sci(e, 2));
    item("rel. error", e / D1 >= 0.001 ? (100 * e / D1).toPrecision(3) + " %" : sci(e / D1, 2));
    item("slope from h to h/2:", Number.isFinite(slope) ? slope.toFixed(2) : "–");
    if (same) item("", (T0 + h === T0 ? "t + h rounds to t: " : "") + "the two OD values are the same double, their difference is 0", "ss-warn");
  }
  draw();
}

// ---------------------------------------------------------------- cost tab
function buildCost(box: HTMLElement, flops: number): void {
  const st = { ln: Math.log10(20000), lm: 5, dens: 10 };
  const ctr = htmlEl("div", "ss-row", "", box);
  const nS = slider(ctr, "genes n", 2, Math.log10(20000), 0.01, st.ln, v => { st.ln = v; draw(); }, v => fmtInt(10 ** v));
  const ng = htmlEl("div", "ss-group", "", ctr);
  const setN = (n: number) => { st.ln = Math.max(2, Math.min(Math.log10(20000), Math.log10(n))); nS.set(st.ln); draw(); };
  const b1 = htmlEl("button", "qb-btn", "n ÷ 2", ng); b1.type = "button";
  b1.addEventListener("click", () => { setN(10 ** st.ln / 2); b1.blur(); });
  const b2 = htmlEl("button", "qb-btn", "n × 2", ng); b2.type = "button";
  b2.addEventListener("click", () => { setN(10 ** st.ln * 2); b2.blur(); });
  presets(ng, [100, 1000, 20000].map(n => ({ label: fmtInt(n), apply: () => setN(n) })));
  const ctr2 = htmlEl("div", "ss-row", "", box);
  const mS = slider(ctr2, "cells m", 3, 5, 0.01, st.lm, v => { st.lm = v; draw(); }, v => fmtInt(10 ** v));
  const dS = slider(ctr2, "non-zero", 1, 70, 1, st.dens, v => { st.dens = v; draw(); }, v => v + " %");
  void mS; void dS;

  const plots = htmlEl("div", "ss-plots", "", box);
  const H = 330, WA = 600, WB = 500, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Time for operation counts" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Memory" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro = htmlEl("div", "qb-readout", "", box);

  const L0 = -7, L1 = 8.5;   // log10 seconds
  const sx = scale(L0, L1, 20, WA - 20);
  const ticks: [number, string][] = [[-6, "1 µs"], [-3, "1 ms"], [0, "1 s"], [Math.log10(3600), "1 h"], [Math.log10(86400), "1 d"], [Math.log10(86400 * 365), "1 yr"]];

  function draw(): void {
    const n = 10 ** st.ln, m = 10 ** st.lm, d = st.dens / 100;
    const rows: [string, number, string, string][] = [
      ["dot product  2n", 2 * n, COL.fwd, "× 2"],
      ["matrix–vector  2n²", 2 * n * n, COL.pur, "× 4"],
      ["gene–gene correlation  n²·m", n * n * m, COL.grn, "× 4"],
      ["matrix–matrix  2n³", 2 * n ** 3, COL.cen, "× 8"],
    ];
    // ---- (c) time ----
    svgA.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent = `(c) time at ${sci(flops, 1)} operations/s (log axis)`;
    const y0 = 52, rowH = 62, yb = y0 + rows.length * rowH;
    for (const [l, lab] of ticks) {
      svgEl("line", { x1: sx(l), x2: sx(l), y1: y0 - 6, y2: yb, class: "grid", "stroke-dasharray": "4 4" }, svgA);
      svgEl("text", { x: sx(l), y: yb + 20, "text-anchor": "middle", class: "sm" }, svgA).textContent = lab;
    }
    rows.forEach(([name, ops, c, dbl], i) => {
      const y = y0 + i * rowH, secs = ops / flops;
      svgEl("text", { x: 20, y: y + 16 }, svgA).textContent = name;
      const w = Math.max(2, sx(Math.min(L1, Math.log10(secs))) - sx(L0));
      svgEl("rect", { x: 20, y: y + 24, width: w, height: 22, fill: c, rx: 3 }, svgA);
      const t = svgEl("text", { x: 20 + w + 8, y: y + 41, "font-weight": 600 }, svgA);
      const txt = `${fmtTime(secs)}  (${sci(ops, 2)} ops)`;
      if (20 + w + 8 + txt.length * 9.5 > WA) { t.setAttribute("x", String(20 + w - 8)); t.setAttribute("text-anchor", "end"); t.style.fill = "#fff"; }
      t.textContent = txt;
      const dt = svgEl("text", { x: WA - 20, y: y + 16, "text-anchor": "end", class: "sm" }, svgA);
      dt.textContent = `n × 2 → time ${dbl}`;
    });

    // ---- (d) memory ----
    svgB.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgB).textContent = "(d) memory, 8 bytes per double";
    const mem: [string, number, string][] = [
      [`dense genes × cells`, 8 * n * m, COL.cen],
      [`sparse, ${st.dens} % non-zero`, 12 * d * n * m + 4 * (m + 1), COL.fwd],
      [`gene × gene correlations`, 8 * n * n, COL.pur],
    ];
    const LAP = 16e9;
    const xmax = Math.max(20e9, ...mem.map(r => r[1])) * 1.05;
    const mx = scale(0, xmax, 20, WB - 20);
    const ybm = y0 + mem.length * 72;
    const step = xmax > 100e9 ? 50e9 : xmax > 40e9 ? 20e9 : 5e9;
    for (let v = 0; v <= xmax; v += step) {
      svgEl("line", { x1: mx(v), x2: mx(v), y1: y0 - 6, y2: ybm, class: v === 0 ? "axis" : "grid" }, svgB);
      svgEl("text", { x: mx(v), y: ybm + 20, "text-anchor": "middle", class: "sm" }, svgB).textContent = String(v / 1e9);
    }
    svgEl("line", { x1: mx(LAP), x2: mx(LAP), y1: y0 - 6, y2: ybm, stroke: COL.ref, "stroke-width": 2, "stroke-dasharray": "6 4" }, svgB);
    svgEl("text", { x: mx(LAP), y: y0 - 12, "text-anchor": "middle", class: "sm" }, svgB).textContent = "16 GB laptop";
    mem.forEach(([name, bytes, c], i) => {
      const y = y0 + i * 72;
      svgEl("text", { x: 20, y: y + 18 }, svgB).textContent = name;
      const w = Math.max(2, mx(bytes) - mx(0));
      svgEl("rect", { x: 20, y: y + 26, width: w, height: 24, fill: c, rx: 3 }, svgB);
      const t = svgEl("text", { x: 20 + w + 8, y: y + 44, "font-weight": 600 }, svgB);
      const txt = fmtGB(bytes);
      if (20 + w + 8 + txt.length * 10 > WB) { t.setAttribute("x", String(20 + w - 8)); t.setAttribute("text-anchor", "end"); t.style.fill = "#fff"; }
      t.textContent = txt;
    });
    svgEl("text", { x: WB / 2, y: H - 6, "text-anchor": "middle", class: "sm" }, svgB).textContent = "GB";

    ro.textContent = `n = ${fmtInt(n)} genes, m = ${fmtInt(m)} cells.  Sparse: 12 bytes per non-zero (value + row index); ` +
      `it beats dense below ${(100 * 8 / 12).toFixed(0)} % non-zero.  ${sci(flops, 1)} operations/s is one core doing simple arithmetic; ` +
      `an optimised library may be 10–100× faster, but the exponent stays.`;
  }
  draw();
}

autoMount(NAME, mount);

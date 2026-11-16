// Monte Carlo lab for the monte-carlo unit.
// Tab 1, "Monte Carlo integration": add random draws in batches (+10 ... +100 000) to one of
// three problems with a known answer, and watch
//  (a) the running estimate with its +-2 SE band closing in on the exact value (log N axis);
//  (b) the relative error |estimate - exact| / exact against N (log-log), with the
//      line sd / sqrt(N) of slope -1/2.
// Problems: area under the OD curve on [0, 10] h (1-D integral, 7.244448); the motif-hit
// overlap test P(K >= 18), 200 hits placed uniformly in 1 Mb with 50 promoters (0.01314);
// the rare event P(X > 5), X ~ N(0, 1) (2.87e-7), plain or by importance sampling from N(5, 1).
// Tab 2, "Curse of dimensionality": the lecture's toy marginal likelihood in d dimensions.
//  (c) relative error at a fixed budget N: grid at prior quantiles (n = floor(N^(1/d)) per
//      axis) against Monte Carlo with N prior draws (expected SE and one seeded run per d);
//  (d) evaluations needed for 1 % relative error, grid against Monte Carlo, with times at 1 µs.
// Auto-mounts on <div data-demo="mclab">. Optional attributes:
//   data-tab="mc|curse"  data-curse="1|0"  data-problem="auc|motif|rare"  data-is="0|1"
//   data-seed="1"  data-d="10"  data-budget="1e5"
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, rng, autoMount } from "../lib/ui";
import {
  ProbKey, Problem, problem, Running, perAxis, gridErr, mcRelSE, mcRun, gridNeed, mcNeed,
} from "./math";

const NAME = "mclab";
// Okabe-Ito, as in make-monte-carlo-figs.R
const COL = { a: "#0072B2", b: "#D55E00", c: "#009E73", pur: "#CC79A7", ref: "#737373", ink: "#222" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .45em; }
[data-demo="${NAME}"] .mc-row { display: flex; flex-wrap: wrap; gap: .4em 1.3em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .mc-group { display: flex; gap: .4em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .mc-row .qb-slider { min-width: 220px; flex: 0 1 300px; grid-template-columns: auto 1fr 2.6em; }
[data-demo="${NAME}"] .mc-plots { display: flex; flex-wrap: wrap; gap: 10px; width: 100%; }
[data-demo="${NAME}"] .mc-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .mc-read { display: flex; flex-wrap: wrap; gap: .2em 1.4em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .mc-read span b { font-weight: 600; }
[data-demo="${NAME}"] .mc-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .mc-ok { color: #00704f; font-weight: 600; }
[data-demo="${NAME}"] .qb-check.off { opacity: .45; }
[data-demo="${NAME}"] .mc-tabs .qb-btn { font-size: .95em; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- formatting ----
const SUP: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const pow10 = (k: number): string => k === 0 ? "1" : k === 1 ? "10" : "10" + String(k).split("").map(c => SUP[c]).join("");
function sci(v: number, d = 3): string {
  if (!Number.isFinite(v)) return "–";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e-3 && a < Math.max(1e5, 10 ** d)) { const t = v.toPrecision(d); if (!t.includes("e")) return t; }
  const [m, e] = v.toExponential(d - 1).split("e");
  return `${m}e${Number(e)}`;
}
const fmtInt = (n: number): string => Math.round(n).toString().replace(/\B(?=(\d{3})+(?!\d))/g, " ");
const pct = (v: number): string => !Number.isFinite(v) ? "–" : v >= 0.1 ? (100 * v).toFixed(0) + " %" : v >= 0.001 ? (100 * v).toPrecision(2) + " %" : sci(100 * v, 2) + " %";
function fmtTime(s: number): string {
  const f = (v: number, u: string) => `${v >= 100 ? v.toFixed(0) : v >= 10 ? v.toFixed(1) : v.toFixed(2)} ${u}`;
  if (!Number.isFinite(s)) return "–";
  if (s < 1e-3) return f(s * 1e6, "µs");
  if (s < 1) return f(s * 1e3, "ms");
  if (s < 60) return f(s, "s");
  if (s < 3600) return f(s / 60, "min");
  if (s < 86400) return f(s / 3600, "h");
  if (s < 365 * 86400) return f(s / 86400, "days");
  const y = s / (365 * 86400);
  return y < 1e4 ? f(y, "years") : `${sci(y, 2)} years`;
}

/** About `k` round tick values covering [lo, hi]. */
function ticks(lo: number, hi: number, k = 5): number[] {
  const raw = (hi - lo) / k, p = Math.pow(10, Math.floor(Math.log10(raw)));
  const st = [1, 2, 2.5, 5, 10].map(m => m * p).find(s => s >= raw) ?? 10 * p;
  const out: number[] = [];
  for (let v = Math.ceil(lo / st) * st; v <= hi + st * 1e-9; v += st) out.push(Math.abs(v) < st * 1e-9 ? 0 : v);
  return out;
}

const NMAXDRAWS = 1e6;
let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const uid = ++mountCount;
  const withCurse = ds.curse !== "0";
  const startProb: ProbKey = ds.problem === "motif" || ds.problem === "rare" ? ds.problem : "auc";
  const startIS = ds.is === "1";
  const startSeed = Math.max(1, Math.min(20, Math.round(Number(ds.seed) || 1)));

  const mcBox = htmlEl("div", "", "");
  const cuBox = htmlEl("div", "", "");
  for (const b of [mcBox, cuBox]) Object.assign(b.style, { display: "flex", flexDirection: "column", gap: ".45em", width: "100%" });
  let curseShown: () => void = () => undefined;
  if (withCurse) {
    const tabs = htmlEl("div", "mc-row mc-tabs", "", host);
    const tb = presets(tabs, [
      { label: "Monte Carlo integration", apply: () => show("mc") },
      { label: "Curse of dimensionality", apply: () => show("curse") },
    ]);
    const show = (t: string) => {
      mcBox.style.display = t === "mc" ? "flex" : "none";
      cuBox.style.display = t === "curse" ? "flex" : "none";
      tb.forEach((b, i) => b.classList.toggle("active", (i === 0) === (t === "mc")));
      if (t === "curse") curseShown();
    };
    host.append(mcBox, cuBox);
    curseShown = buildCurse(cuBox, ds, uid);
    show(ds.tab === "curse" ? "curse" : "mc");
  } else host.append(mcBox);

  // ================================================================ tab 1
  const st = { prob: startProb, is: startIS, seed: startSeed };
  let P: Problem = problem(st.prob, st.is);
  let run = new Running();
  let r = rng(st.seed);

  const row1 = htmlEl("div", "mc-row", "", mcBox);
  const pg = htmlEl("div", "mc-group", "", row1);
  htmlEl("span", "", "problem:", pg);
  const keys: ProbKey[] = ["auc", "motif", "rare"];
  const pBtns = presets(pg, [
    { label: "OD-curve area", apply: () => setProb("auc") },
    { label: "motif overlap P(K ≥ 18)", apply: () => setProb("motif") },
    { label: "rare event P(X > 5)", apply: () => setProb("rare") },
  ]);
  const isBox = checkbox(row1, "importance sampling, g = N(5, 1)", st.is, v => { st.is = v; reset(); });
  const isLab = isBox.parentElement as HTMLElement;

  const row2 = htmlEl("div", "mc-row", "", mcBox);
  const ag = htmlEl("div", "mc-group", "", row2);
  htmlEl("span", "", "add draws:", ag);
  const addBtns = [10, 100, 1000, 10000, 100000].map(k => button(ag, "+" + fmtInt(k), () => add(k)));
  button(ag, "reset", () => { st.prob = startProb; st.is = startIS; st.seed = startSeed; isBox.checked = st.is; seedS.set(st.seed); reset(); });
  const seedS = slider(row2, "seed", 1, 20, 1, st.seed, v => { st.seed = v; reset(); });

  const plots = htmlEl("div", "mc-plots", "", mcBox);
  const H = 360, WA = 550, WB = 550, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Running estimate with ±2 SE band" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Relative error against N, log-log" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro = htmlEl("div", "qb-readout mc-read", "", mcBox);

  const top = 40, bot = H - 52, L0 = 0, L1 = 6;
  const ax = scale(L0, L1, 84, WA - 14);
  const bx = scale(L0, L1, 70, WB - 14);
  const EY0 = -5, EY1 = 1;
  const by = scale(EY0, EY1, bot, top);

  function setProb(k: ProbKey): void { st.prob = k; reset(); }
  function reset(): void {
    P = problem(st.prob, st.is);
    run = new Running();
    r = rng(st.seed);
    draw();
  }
  function add(k: number): void {
    const m = Math.min(k, NMAXDRAWS - run.n);
    for (let i = 0; i < m; i++) {
      const y = P.draw(r);
      run.add(y, P.hit ? P.hit(y) : false);
    }
    run.mark();
    draw();
  }

  function xAxis(svg: SVGSVGElement, x: (v: number) => number, label: string): void {
    svgEl("line", { x1: x(L0), x2: x(L1), y1: bot, y2: bot, class: "axis" }, svg);
    for (let k = L0; k <= L1; k++) {
      svgEl("line", { x1: x(k), x2: x(k), y1: top, y2: bot, class: "grid" }, svg);
      svgEl("text", { x: x(k), y: bot + 22, "text-anchor": "middle", class: "sm" }, svg).textContent = pow10(k);
    }
    svgEl("text", { x: (x(L0) + x(L1)) / 2, y: H - 8, "text-anchor": "middle" }, svg).textContent = label;
  }
  function yLabel(svg: SVGSVGElement, txt: string): void {
    const t = svgEl("text", { x: 16, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(top + bot) / 2})` }, svg);
    t.textContent = txt;
  }
  function clipped(svg: SVGSVGElement, id: string, x0: number, x1: number): SVGGElement {
    const cp = svgEl("clipPath", { id }, svg);
    svgEl("rect", { x: x0, y: top - 2, width: x1 - x0, height: bot - top + 4 }, cp);
    return svgEl("g", { "clip-path": `url(#${id})` }, svg);
  }

  function draw(): void {
    pBtns.forEach((b, i) => b.classList.toggle("active", keys[i] === st.prob));
    isBox.disabled = st.prob !== "rare";
    isLab.classList.toggle("off", st.prob !== "rare");
    addBtns.forEach(b => { b.disabled = run.n >= NMAXDRAWS; });
    const hist = run.hist, ex = P.exact, n = run.n;

    // ---------- (a) running estimate ----------
    svgA.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent = "(a) running estimate ± 2 SE";
    // y range from the points with N >= 10, always including the exact value
    const use = hist.filter(h => h.n >= 10);
    let lo = ex, hi = ex;
    for (const h of use.length ? use : hist) { lo = Math.min(lo, h.est - 2 * h.se); hi = Math.max(hi, h.est + 2 * h.se); }
    if (P.key !== "auc") lo = Math.max(lo, 0);
    if (hi - lo < 0.02 * Math.abs(ex)) { lo = ex * 0.97; hi = ex * 1.03; }
    const pad = 0.06 * (hi - lo);
    lo -= pad; hi += pad;
    if (P.key !== "auc" && lo < 0) lo = Math.min(0, lo);
    const ay = scale(lo, hi, bot, top);
    // tiny probabilities are labelled in units of 10^k
    const ue = Math.abs(ex) < 1e-3 ? Math.floor(Math.log10(Math.abs(ex))) : 0, u = 10 ** ue;
    const tv = ticks(lo, hi, 5), stp = tv.length > 1 ? (tv[1] - tv[0]) / u : 1;
    const dec = Math.max(0, Math.ceil(-Math.log10(stp) - 1e-9) + (Math.abs(stp / 10 ** Math.floor(Math.log10(stp)) - 2.5) < 1e-9 ? 1 : 0));
    for (const v of tv) {
      svgEl("line", { x1: ax(L0), x2: ax(L1), y1: ay(v), y2: ay(v), class: "grid" }, svgA);
      svgEl("text", { x: ax(L0) - 6, y: ay(v) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = (v / u).toFixed(dec);
    }
    xAxis(svgA, ax, "number of draws N");
    yLabel(svgA, ue ? `estimate (× ${pow10(ue)})` : "estimate");
    const ga = clipped(svgA, "mc-clipA-" + uid, ax(L0), ax(L1));
    const lx = (k: number) => ax(Math.log10(k));
    const band = hist.filter(h => h.n >= 2);
    if (band.length > 1) {
      const up = band.map(h => `${lx(h.n)},${ay(h.est + 2 * h.se)}`);
      const dn = band.slice().reverse().map(h => `${lx(h.n)},${ay(h.est - 2 * h.se)}`);
      svgEl("polygon", { points: up.concat(dn).join(" "), fill: COL.a, "fill-opacity": 0.2, stroke: "none" }, ga);
    }
    svgEl("line", { x1: ax(L0), x2: ax(L1), y1: ay(ex), y2: ay(ex), stroke: COL.c, "stroke-width": 3, "stroke-dasharray": "10 6" }, ga);
    if (hist.length > 1) svgEl("polyline", { points: hist.map(h => `${lx(h.n)},${ay(h.est)}`).join(" "), fill: "none", stroke: COL.a, "stroke-width": 2.5 }, ga);
    if (n > 0) svgEl("circle", { cx: lx(n), cy: ay(run.mean), r: 7, fill: COL.a, stroke: "#fff", "stroke-width": 2 }, ga);
    const ext = svgEl("text", { x: ax(L1) - 6, y: ay(ex) - 8, "text-anchor": "end", "font-weight": 600 }, svgA);
    ext.style.fill = COL.c; ext.textContent = `exact ${sci(ex, 4)}`;
    if (n === 0) svgEl("text", { x: (ax(L0) + ax(L1)) / 2, y: (top + bot) / 2 + 30, "text-anchor": "middle", class: "sm" }, svgA).textContent = "click +10 to draw the first points";

    // ---------- (b) error vs N ----------
    svgB.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgB).textContent = "(b) relative error (log–log)";
    for (let k = EY0; k <= EY1; k++) {
      svgEl("line", { x1: bx(L0), x2: bx(L1), y1: by(k), y2: by(k), class: "grid" }, svgB);
      svgEl("text", { x: bx(L0) - 6, y: by(k) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = pow10(k);
    }
    xAxis(svgB, bx, "number of draws N");
    yLabel(svgB, "|estimate − exact| / exact");
    const gb = clipped(svgB, "mc-clipB-" + uid, bx(L0), bx(L1));
    const clampE = (e: number) => Math.max(EY0 - 0.3, Math.log10(Math.max(e, 1e-30)));
    const relSd = P.sd / ex;
    const ref = Array.from({ length: 61 }, (_, i) => L0 + i * (L1 - L0) / 60).map(l => `${bx(l)},${by(clampE(relSd / Math.sqrt(10 ** l)))}`).join(" ");
    svgEl("polyline", { points: ref, fill: "none", stroke: COL.ref, "stroke-width": 2.5, "stroke-dasharray": "9 6" }, gb);
    const pts = hist.map(h => [bx(Math.log10(h.n)), by(clampE(Math.abs(h.est - ex) / ex))]);
    if (pts.length > 1) svgEl("polyline", { points: pts.map(p => p.join(",")).join(" "), fill: "none", stroke: COL.b, "stroke-width": 1.5, "stroke-opacity": 0.6 }, gb);
    pts.forEach(([x, y]) => svgEl("circle", { cx: x, cy: y, r: 3.5, fill: COL.b }, gb));
    const lr = svgEl("text", { x: bx(L0) + 10, y: by(EY0) - 38, "font-weight": 600 }, svgB);
    lr.style.fill = "#555"; lr.textContent = "- - true SE σ/√N: slope −½";
    const lg = svgEl("text", { x: bx(L0) + 10, y: by(EY0) - 12, "font-weight": 600 }, svgB);
    lg.style.fill = COL.b; lg.textContent = "● actual error of this run";

    // ---------- readout ----------
    ro.replaceChildren();
    const item = (k: string, v: string, cls = "") => {
      const s = htmlEl("span", cls, "", ro);
      if (k) htmlEl("b", "", k + " ", s);
      s.append(v);
    };
    item("", P.label);
    item("N =", fmtInt(n));
    if (n > 0) {
      const se = run.se, err = run.mean - ex;
      item("estimate", `${sci(run.mean, 4)} ± ${sci(2 * se, 2)}`);
      item("SE", sci(se, 3));
      item("exact", sci(ex, 4));
      item("error", `${sci(err, 2)} (${pct(Math.abs(err) / ex)})`);
      if (P.hit) item(P.key === "motif" ? "placements with K ≥ 18:" : "draws above 5:", fmtInt(run.hits));
      if (n > 1) {
        const inside = Math.abs(err) <= 2 * se;
        item("", inside ? "exact value inside ± 2 SE" : se === 0 ? "SE = 0: the error bar claims certainty" : "exact value outside ± 2 SE",
          inside ? "mc-ok" : "mc-warn");
      }
    }
    if (n >= NMAXDRAWS) item("", "maximum of 10⁶ draws reached: reset or change the seed");
  }
  reset();
}

// ================================================================ tab 2
function buildCurse(box: HTMLElement, ds: DOMStringMap, uid: number): () => void {
  const D1 = 20;
  const st = {
    d: Math.max(1, Math.min(D1, Math.round(Number(ds.d) || 10))),
    N: [1e3, 1e4, 1e5].includes(Number(ds.budget)) ? Number(ds.budget) : 1e5,
    seed: 1,
  };
  const start = { ...st };
  const TOL = 0.01, US = 1e-6;   // target relative error; 1 µs per evaluation

  const row = htmlEl("div", "mc-row", "", box);
  const dS = slider(row, "dimension d", 1, D1, 1, st.d, v => { st.d = v; draw(); });
  const bg = htmlEl("div", "mc-group", "", row);
  htmlEl("span", "", "budget N =", bg);
  const budgets = [1e3, 1e4, 1e5];
  const bBtns = presets(bg, budgets.map(b => ({ label: pow10(Math.log10(b)), apply: () => { st.N = b; recompute(); } })));
  const sS = slider(row, "seed", 1, 20, 1, st.seed, v => { st.seed = v; recompute(); });
  button(row, "reset", () => { Object.assign(st, start); dS.set(st.d); sS.set(st.seed); recompute(); });

  const plots = htmlEl("div", "mc-plots", "", box);
  const H = 360, WA = 550, WB = 550, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Relative error against dimension at a fixed budget" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Evaluations needed for 1% error against dimension" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro = htmlEl("div", "qb-readout mc-read", "", box);

  const top = 40, bot = H - 52;
  const ax = scale(0.5, D1 + 0.5, 70, WA - 14);
  const EY0 = -8, EY1 = 1;
  const ay = scale(EY0, EY1, bot, top);
  const bx = scale(0.5, D1 + 0.5, 70, WB - 14);
  const NY1 = 42;
  const byy = scale(0, NY1, bot, top);
  const clampE = (e: number) => Math.max(EY0, Math.min(EY1, Math.log10(Math.max(e, 1e-30))));

  // seeded Monte Carlo runs, one per d, computed a few at a time so the page stays responsive
  let runs: ({ err: number; relse: number } | undefined)[] = [];
  let token = 0;
  const needs: (ReturnType<typeof gridNeed> | null)[] = Array(D1 + 1).fill(null);

  function recompute(): void {
    runs = Array(D1 + 1).fill(undefined);
    const my = ++token;
    const order = [st.d, ...Array.from({ length: D1 }, (_, i) => i + 1).filter(d => d !== st.d)];
    const step = () => {
      if (my !== token) return;
      const t0 = performance.now();
      while (order.length && performance.now() - t0 < 40) {
        const d = order.shift() as number;
        runs[d] = mcRun(st.N, d, st.seed);
      }
      draw();
      if (order.length) setTimeout(step, 0);
    };
    draw();
    setTimeout(step, 0);
  }

  function draw(): void {
    bBtns.forEach((b, i) => b.classList.toggle("active", budgets[i] === st.N));
    const d = st.d, N = st.N;

    // ---------- (c) error at a fixed budget ----------
    svgA.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent = `(c) relative error with N = ${pow10(Math.log10(N))} evaluations`;
    for (let k = EY0; k <= EY1; k++) {
      svgEl("line", { x1: ax(0.5), x2: ax(D1 + 0.5), y1: ay(k), y2: ay(k), class: "grid" }, svgA);
      if ((k - EY0) % 2 === 0 || k === EY1) svgEl("text", { x: ax(0.5) - 6, y: ay(k) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = pow10(k);
    }
    svgEl("line", { x1: ax(0.5), x2: ax(D1 + 0.5), y1: bot, y2: bot, class: "axis" }, svgA);
    for (let k = 1; k <= D1; k++) if (k === 1 || k % 2 === 0) svgEl("text", { x: ax(k), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgA).textContent = String(k);
    svgEl("text", { x: (70 + WA) / 2, y: H - 8, "text-anchor": "middle" }, svgA).textContent = "dimension d";
    const yl = svgEl("text", { x: 16, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(top + bot) / 2})` }, svgA);
    yl.textContent = "|Ẑ − Z| / Z";
    svgEl("rect", { x: ax(d - 0.5), y: top, width: ax(1) - ax(0), height: bot - top, fill: "#000", "fill-opacity": 0.07 }, svgA);
    const ds_ = Array.from({ length: D1 }, (_, i) => i + 1);
    const gridPts = ds_.map(k => [ax(k), ay(clampE(gridErr(perAxis(N, k), k)))]);
    svgEl("polyline", { points: gridPts.map(p => p.join(",")).join(" "), fill: "none", stroke: COL.b, "stroke-width": 2.5 }, svgA);
    gridPts.forEach(([x, y]) => svgEl("path", { d: `M${x},${y - 7} l6.5,11 h-13 z`, fill: COL.b }, svgA));
    svgEl("polyline", { points: ds_.map(k => `${ax(k)},${ay(clampE(mcRelSE(N, k)))}`).join(" "), fill: "none", stroke: COL.a, "stroke-width": 2.5, "stroke-dasharray": "9 6" }, svgA);
    ds_.forEach(k => { const rr = runs[k]; if (rr) svgEl("circle", { cx: ax(k), cy: ay(clampE(rr.err)), r: 6, fill: COL.a, stroke: "#fff", "stroke-width": 1.5 }, svgA); });
    const lg = (txt: string, y: number, c: string) => { const t = svgEl("text", { x: ax(6.5), y, "font-weight": 600 }, svgA); t.style.fill = c; t.textContent = txt; };
    lg("▲ grid, n = ⌊N^(1/d)⌋ per axis", ay(EY0) - 52, COL.b);
    lg("- - Monte Carlo, expected SE", ay(EY0) - 30, COL.a);
    lg("● Monte Carlo, one run", ay(EY0) - 8, COL.a);

    // ---------- (d) evaluations for 1 % ----------
    svgB.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgB).textContent = "(d) evaluations needed for 1 % error";
    for (let k = 0; k <= 40; k += 10) {
      svgEl("line", { x1: bx(0.5), x2: bx(D1 + 0.5), y1: byy(k), y2: byy(k), class: "grid" }, svgB);
      svgEl("text", { x: bx(0.5) - 6, y: byy(k) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = pow10(k);
    }
    svgEl("line", { x1: bx(0.5), x2: bx(D1 + 0.5), y1: bot, y2: bot, class: "axis" }, svgB);
    for (let k = 1; k <= D1; k++) if (k === 1 || k % 2 === 0) svgEl("text", { x: bx(k), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgB).textContent = String(k);
    svgEl("text", { x: (70 + WB) / 2, y: H - 8, "text-anchor": "middle" }, svgB).textContent = "dimension d";
    const yl2 = svgEl("text", { x: 16, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(top + bot) / 2})` }, svgB);
    yl2.textContent = "evaluations";
    svgEl("rect", { x: bx(d - 0.5), y: top, width: bx(1) - bx(0), height: bot - top, fill: "#000", "fill-opacity": 0.07 }, svgB);
    // times at 1 µs per evaluation
    const tl: [number, string][] = [[Math.log10(1 / US), "1 s"], [Math.log10(3600 / US), "1 h"], [Math.log10(3.156e7 / US), "1 year"], [Math.log10(4.35e17 / US), "age of the universe"]];
    for (const [v, s] of tl) {
      svgEl("line", { x1: bx(0.5), x2: bx(D1 + 0.5), y1: byy(v), y2: byy(v), stroke: COL.ref, "stroke-width": 1.5, "stroke-dasharray": "3 5" }, svgB);
      svgEl("text", { x: bx(0.5) + 6, y: byy(v) - 5, class: "sm" }, svgB).textContent = s + " at 1 µs";
    }
    for (let k = 1; k <= D1; k++) if (!needs[k]) needs[k] = gridNeed(k, TOL);
    const gp = ds_.map(k => [bx(k), byy(Math.min(NY1, (needs[k] as { log10N: number }).log10N))]);
    svgEl("polyline", { points: gp.map(p => p.join(",")).join(" "), fill: "none", stroke: COL.b, "stroke-width": 2.5 }, svgB);
    gp.forEach(([x, y]) => svgEl("path", { d: `M${x},${y - 7} l6.5,11 h-13 z`, fill: COL.b }, svgB));
    const mp = ds_.map(k => [bx(k), byy(Math.max(0, mcNeed(k, TOL)))]);
    svgEl("polyline", { points: mp.map(p => p.join(",")).join(" "), fill: "none", stroke: COL.a, "stroke-width": 2.5 }, svgB);
    mp.forEach(([x, y]) => svgEl("circle", { cx: x, cy: y, r: 5, fill: COL.a }, svgB));
    const lg2 = (txt: string, y: number, c: string) => { const t = svgEl("text", { x: bx(0.5) + 10, y, "font-weight": 600 }, svgB); t.style.fill = c; t.textContent = txt; };
    lg2("▲ grid (prior quantiles)", byy(NY1) + 22, COL.b);
    lg2("● Monte Carlo (SE = 1 %)", byy(NY1) + 44, COL.a);

    // ---------- readout ----------
    const n = perAxis(N, d), nd = Math.pow(n, d);
    const ge = gridErr(n, d), rr = runs[d], nd1 = needs[d] as { n: number; log10N: number };
    const mcN = Math.pow(10, mcNeed(d, TOL));
    ro.replaceChildren();
    const item = (k: string, v: string) => { const s = htmlEl("span", "", "", ro); if (k) htmlEl("b", "", k + " ", s); s.append(v); };
    item(`d = ${d}:`, `grid ${n} per axis = ${fmtInt(nd)} points,`);
    item("grid error", ge < 1e-8 ? "< 1e-8" : pct(ge));
    item("MC expected SE", pct(mcRelSE(N, d)));
    item("MC run", rr ? `error ${pct(rr.err)}, SE ${pct(rr.relse)}` : "computing …");
    const gN = nd1.log10N < 7 ? fmtInt(Math.pow(10, nd1.log10N)) : "≈ " + pow10(Math.round(nd1.log10N));
    item("for 1 %:", `grid ${nd1.n} per axis = ${gN} evaluations (${fmtTime(Math.pow(10, nd1.log10N) * US)});`);
    item("Monte Carlo", `${fmtInt(Math.ceil(mcN))} draws (${fmtTime(mcN * US)})`);
  }

  let started = false;
  return () => { if (!started) { started = true; recompute(); } else draw(); };
}

autoMount(NAME, mount);

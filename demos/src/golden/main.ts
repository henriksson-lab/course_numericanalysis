// 1D minimiser step-through for the optimisation-1d unit. Three objective functions of the
// lecture: the mRNA decay SSE(k), the circadian RSS(T) (two local minima) and the
// zero-truncated Poisson minus log-likelihood. Three methods: golden-section search,
// successive parabolic interpolation and Brent's method (R's optimize, step for step).
//  (a) f with the bracket [a, b] (drag a and b), the interior points, the parabola of a
//      parabolic step, every evaluated point; "plot the curve" hides/shows f, "multistart"
//      runs Brent on six sub-intervals; "flat bottom" zooms in on the minimum, where f values
//      are identical doubles over a band of width about sqrt(eps).
//  (b) bracket width and error of the current estimate against the number of evaluations
//      (log scale), with 0.618^n and the flat-bottom floor; a table of every step.
// Auto-mounts on <div data-demo="golden">. Optional attributes:
//   data-problem="decay|circadian|ztp"  data-method="golden|parabolic|brent"
//   data-tol="default|1e-10|1e-15"  data-a, data-b (initial bracket)  data-curve="1|0"
//   data-zoom="1" (start in the flat-bottom view)
import { svgEl, htmlEl, button, presets, scale, svgPoint, draggable, autoMount } from "../lib/ui";
import {
  PROBLEMS, ProbKey, Problem, Method, Run, Pt, makeRun, runAll, brent, parabola, vertex, flatWidth,
  fmtSci, fx, EPS, G,
} from "./math";

const NAME = "golden";
// Okabe-Ito, as in make-optimisation-1d-figs.R
const COL = {
  curve: "#0072B2", gold: "#E69F00", par: "#D55E00", brent: "#0072B2", grn: "#009E73",
  pur: "#CC79A7", sky: "#56B4E9", grey: "#999999", ink: "#222222", width: "#555555",
};
const MCOL: Record<Method, string> = { golden: COL.gold, parabolic: COL.pur, brent: COL.brent };
const MNAME: Record<Method, string> = { golden: "golden section", parabolic: "parabolic interpolation", brent: "Brent (optimize)" };
const TOLS: { label: string; v: number }[] = [
  { label: "ε^¼ ≈ 1.2e-4 (R default)", v: Math.pow(EPS, 0.25) },
  { label: "1e-10", v: 1e-10 },
  { label: "1e-15", v: 1e-15 },
];

const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .gs-row { display: flex; flex-wrap: wrap; gap: .35em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .gs-group { display: flex; flex-wrap: wrap; gap: .35em; align-items: center; }
[data-demo="${NAME}"] .gs-plots { display: flex; flex-wrap: wrap; gap: 10px; width: 100%; align-items: flex-start; }
[data-demo="${NAME}"] .gs-main { flex: 1 1 520px; min-width: 0; max-width: 680px; }
[data-demo="${NAME}"] .gs-side { flex: 1 1 400px; min-width: 0; max-width: 560px; display: flex; flex-direction: column; gap: 6px; }
[data-demo="${NAME}"] svg { width: 100%; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #444; }
[data-demo="${NAME}"] .gs-table { max-height: 190px; overflow-y: auto; border: 1px solid #ddd; border-radius: 6px; }
[data-demo="${NAME}"] .gs-table table { border-collapse: collapse; width: 100%; font-family: ui-monospace, monospace; font-size: 17px; }
[data-demo="${NAME}"] .gs-table th { position: sticky; top: 0; background: #f0f1f4; font-weight: 600; text-align: right; padding: 1px 6px; }
[data-demo="${NAME}"] .gs-table td { text-align: right; padding: 0 6px; white-space: nowrap; }
[data-demo="${NAME}"] .gs-table tr.last td { background: #fff3d6; }
[data-demo="${NAME}"] .gs-table td.k { text-align: left; }
[data-demo="${NAME}"] .gs-legend { display: flex; flex-wrap: wrap; gap: .1em 1em; font-size: 17px; }
[data-demo="${NAME}"] .gs-sw { display: inline-block; width: 22px; height: 0; margin-right: 5px; vertical-align: middle; }
[data-demo="${NAME}"] .gs-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .gs-ok { color: #00704f; font-weight: 600; }
[data-demo="${NAME}"] .qb-readout { line-height: 1.4; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .gs-lab { font-weight: 600; }
[data-demo="${NAME}"] .qb-btn.on { outline: 2px solid ${COL.grn}; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; gap: .25em; }
.reveal [data-demo="${NAME}"] .gs-plots { flex-wrap: nowrap; }
.reveal [data-demo="${NAME}"] .gs-main { flex: 0 0 600px; max-width: 600px; }
.reveal [data-demo="${NAME}"] .gs-side { flex: 1 1 auto; }
.reveal [data-demo="${NAME}"] .gs-table { max-height: 150px; }
.reveal [data-demo="${NAME}"] .gs-table table { font-size: 16px; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; padding: .3em .5em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

const sup: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const pow10 = (k: number): string => k === 0 ? "1" : "10" + String(k).split("").map(c => sup[c]).join("");
const ulp = (v: number): number => Math.pow(2, Math.floor(Math.log2(Math.abs(v))) - 52);

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const uid = ++mountCount;
  const pk = (s: string | undefined): ProbKey => (s === "circadian" || s === "ztp" ? s : "decay");
  const mk = (s: string | undefined): Method => (s === "parabolic" || s === "brent" ? s : "golden");
  const tolIdx = (s: string | undefined): number => s === "1e-10" ? 1 : s === "1e-15" ? 2 : 0;
  const num = (s: string | undefined): number => (s !== undefined && s !== "" && Number.isFinite(Number(s)) ? Number(s) : NaN);

  const init = { prob: pk(ds.problem), method: mk(ds.method), tol: tolIdx(ds.tol), zoom: ds.zoom === "1" };
  const st = {
    prob: init.prob, method: init.method, tol: init.tol,
    a0: 0, b0: 1, curve: true, zoom: init.zoom, lw: -7.5,   // lw: log10 half-width of the zoom
    run: null as Run | null, multi: null as { pts: Pt[]; evals: number } | null,
  };
  let timer: number | null = null;

  function setProblem(k: ProbKey, first = false): void {
    st.prob = k;
    const p = PROBLEMS[k];
    st.a0 = p.bracket[0]; st.b0 = p.bracket[1];
    if (first) {
      const a = num(ds.a), b = num(ds.b);
      if (Number.isFinite(a) && Number.isFinite(b) && a < b) { st.a0 = Math.max(p.view[0], a); st.b0 = Math.min(p.view[1], b); }
    }
    st.curve = first && ds.curve !== undefined ? ds.curve !== "0" : k !== "circadian";
    st.multi = null;
    st.lw = k === "decay" ? -8.5 : -7.5;
    newRun();
  }
  function newRun(): void {
    stop();
    const p = PROBLEMS[st.prob];
    st.run = makeRun(st.method, p.f, st.a0, st.b0, TOLS[st.tol].v);
  }
  function stop(): void {
    if (timer !== null) { clearInterval(timer); timer = null; }
    runBtn.textContent = "Run ▶";
  }

  // ---- controls ----
  const row1 = htmlEl("div", "gs-row", "", host);
  const g1 = htmlEl("div", "gs-group", "", row1);
  htmlEl("span", "gs-lab", "problem:", g1);
  const keys: ProbKey[] = ["decay", "circadian", "ztp"];
  const pBtns = presets(g1, keys.map(k => ({ label: PROBLEMS[k].label, apply: () => { setProblem(k); draw(); } })));
  const g2 = htmlEl("div", "gs-group", "", row1);
  htmlEl("span", "gs-lab", "method:", g2);
  const methods: Method[] = ["golden", "brent", "parabolic"];
  const mBtns = presets(g2, methods.map(m => ({
    label: m === "golden" ? "golden section" : m === "brent" ? "Brent" : "parabolic",
    apply: () => { st.method = m; newRun(); draw(); },
  })));

  const row2 = htmlEl("div", "gs-row", "", host);
  const g3 = htmlEl("div", "gs-group", "", row2);
  button(g3, "Step", () => { stop(); st.run?.step(); draw(); });
  const runBtn = button(g3, "Run ▶", () => {
    if (timer !== null) { stop(); return; }
    if (st.run?.done) newRun();
    runBtn.textContent = "Pause ❚❚";
    timer = window.setInterval(() => {
      if (!st.run || st.run.done) { stop(); return; }
      st.run.step(); draw();
      if (st.run.done) stop();
    }, 320);
  });
  button(g3, "To end ⇥", () => { stop(); if (st.run) runAll(st.run); draw(); });
  button(g3, "Reset", () => {
    st.method = init.method; st.tol = init.tol; st.zoom = init.zoom;
    setProblem(init.prob, true); draw();
  });
  const g4 = htmlEl("div", "gs-group", "", row2);
  htmlEl("span", "gs-lab", "tol:", g4);
  const tBtns = presets(g4, TOLS.map((t, i) => ({ label: t.label, apply: () => { st.tol = i; newRun(); draw(); } })));
  const g5 = htmlEl("div", "gs-group", "", row2);
  const curveBtn = button(g5, "plot the curve", () => { st.curve = !st.curve; draw(); });
  const multiBtn = button(g5, "multistart", () => {
    if (st.multi) { st.multi = null; draw(); return; }
    const p = PROBLEMS[st.prob];
    const L = st.a0, R = st.b0, n = 5, pts: Pt[] = [];
    let evals = 0;
    for (let i = 0; i < n; i++) {
      const lo = L + i * (R - L) / n, hi = L + (i + 1) * (R - L) / n;
      const r = runAll(brent(p.f, lo, hi, 1e-10));
      evals += r.evals.length;
      const x = r.rows[r.rows.length - 1].best, f = p.f(x);
      const atEdge = Math.min(x - lo, hi - x) < 1e-3 * (hi - lo);    // not a minimum inside
      if (!atEdge && !pts.some(q => Math.abs(q.x - x) < 1e-4 * (R - L))) pts.push({ x, f });
    }
    st.multi = { pts: pts.sort((u, v) => u.x - v.x), evals };
    st.curve = true;
    draw();
  });
  const zoomBtn = button(g5, "zoom: flat bottom", () => { st.zoom = !st.zoom; draw(); });
  const zg = htmlEl("span", "gs-group", "", g5);
  const zIn = button(zg, "zoom in ×10", () => { st.lw = Math.max(-11, st.lw - 1); draw(); });
  const zOut = button(zg, "zoom out ×10", () => { st.lw = Math.min(-2, st.lw + 1); draw(); });

  // ---- plots ----
  const plots = htmlEl("div", "gs-plots", "", host);
  const mainBox = htmlEl("div", "gs-main", "", plots);
  const side = htmlEl("div", "gs-side", "", plots);
  const W = 620, H = 400;
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Objective function with the bracket and the points of the minimiser" }, mainBox);
  const WL = 520, HL = 250;
  const svgL = svgEl("svg", { viewBox: `0 0 ${WL} ${HL}`, role: "img", "aria-label": "Bracket width and error against the number of evaluations, log scale" }, side);
  const legend = htmlEl("div", "gs-legend", "", side);
  const tbox = htmlEl("div", "gs-table", "", side);
  const ro = htmlEl("div", "qb-readout", "", host);

  const L = 74, Rm = W - 14, top = 36, bot = H - 66;

  // drag handles for a and b (created once, repositioned on draw)
  const handles: SVGCircleElement[] = [];
  const hLabels: SVGTextElement[] = [];
  const hLayer = svgEl("g", {});
  for (let i = 0; i < 2; i++) {
    const c = svgEl("circle", { r: 12, class: "handle", stroke: COL.ink }, hLayer);
    handles.push(c);
    const t = svgEl("text", { "text-anchor": "middle", "font-weight": 700 }, hLayer);
    t.style.pointerEvents = "none";
    hLabels.push(t);
    draggable(c, ev => svgPoint(svg, ev), (px) => {
      if (st.zoom) return;
      const p = PROBLEMS[st.prob];
      const sx = scale(p.view[0], p.view[1], L, Rm);
      const minW = 0.04 * (p.view[1] - p.view[0]);
      let v = sx.invert(px);
      v = Math.max(p.view[0], Math.min(p.view[1], v));
      const dec = p.key === "circadian" ? 1 : 2;
      v = Number(v.toFixed(dec));
      if (i === 0) st.a0 = Math.min(v, st.b0 - minW); else st.b0 = Math.max(v, st.a0 + minW);
      st.multi = null;
      newRun(); draw();
    });
  }

  function draw(): void {
    const p = PROBLEMS[st.prob];
    const run = st.run!;
    pBtns.forEach((b, i) => b.classList.toggle("active", keys[i] === st.prob));
    mBtns.forEach((b, i) => b.classList.toggle("active", methods[i] === st.method));
    tBtns.forEach((b, i) => b.classList.toggle("active", i === st.tol));
    curveBtn.textContent = st.curve ? "hide the curve" : "plot the curve";
    curveBtn.classList.toggle("on", st.curve);
    multiBtn.classList.toggle("on", !!st.multi);
    zoomBtn.textContent = st.zoom ? "full curve" : "zoom: flat bottom";
    zoomBtn.classList.toggle("on", st.zoom);
    zg.style.display = st.zoom ? "flex" : "none";
    curveBtn.style.display = multiBtn.style.display = st.zoom ? "none" : "";
    zIn.disabled = st.lw <= -11; zOut.disabled = st.lw >= -2;
    if (st.zoom) drawZoom(p, run); else drawMain(p, run);
    drawLog(p, run);
    drawTable(p, run);
    drawReadout(p, run);
  }

  // ---------------------------------------------------------------- (a) the function
  function drawMain(p: Problem, run: Run): void {
    svg.replaceChildren();
    const sx = scale(p.view[0], p.view[1], L, Rm), sy = scale(p.ylim[0], p.ylim[1], bot, top);
    const cy = (v: number) => sy(Math.max(p.ylim[0] - 0.08 * (p.ylim[1] - p.ylim[0]), Math.min(p.ylim[1] + 0.05 * (p.ylim[1] - p.ylim[0]), v)));
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svg).textContent = `(a) ${p.fname}: ${MNAME[run.method]}`;
    // axes
    for (const v of niceTicks(p.ylim[0], p.ylim[1], 5)) {
      svgEl("line", { x1: L, x2: Rm, y1: sy(v), y2: sy(v), class: "grid" }, svg);
      svgEl("text", { x: L - 6, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = tickFmt(v);
    }
    svgEl("line", { x1: L, x2: Rm, y1: bot, y2: bot, class: "axis" }, svg);
    for (const v of niceTicks(p.view[0], p.view[1], 7)) {
      svgEl("line", { x1: sx(v), x2: sx(v), y1: bot, y2: bot + 5, class: "axis" }, svg);
      svgEl("text", { x: sx(v), y: bot + 40, "text-anchor": "middle", class: "sm" }, svg).textContent = tickFmt(v);
    }
    svgEl("text", { x: (L + Rm) / 2, y: H - 4, "text-anchor": "middle" }, svg).textContent = p.xname;
    const yl = svgEl("text", { x: 18, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 18 ${(top + bot) / 2})` }, svg);
    yl.textContent = p.fname;

    const clip = svgEl("clipPath", { id: `gs-clip-${uid}` }, svg);
    svgEl("rect", { x: L, y: top - 6, width: Rm - L, height: bot - top + 6 }, clip);
    const pg = svgEl("g", { "clip-path": `url(#gs-clip-${uid})` }, svg);

    // bracket
    const br = run.bracket();
    if (br) {
      svgEl("rect", { x: sx(br[0]), y: top, width: Math.max(1.5, sx(br[1]) - sx(br[0])), height: bot - top, fill: COL.sky, "fill-opacity": 0.2 }, pg);
    }
    // the part golden section will throw away next
    if (run.method === "golden" && run.evals.length && !run.done && br) {
      const [q1, q2] = run.inner();
      const [d0, d1] = q1.f < q2.f ? [q2.x, br[1]] : [br[0], q1.x];
      svgEl("rect", { x: sx(d0), y: top, width: Math.max(1, sx(d1) - sx(d0)), height: bot - top, fill: COL.par, "fill-opacity": 0.12 }, pg);
      svgEl("text", { x: (sx(d0) + sx(d1)) / 2, y: top + 20, "text-anchor": "middle", class: "sm" }, pg).textContent =
        sx(d1) - sx(d0) > 70 ? "discard next" : "";
    }
    // the curve
    if (st.curve) {
      const n = 600, pts: string[] = [];
      for (let i = 0; i <= n; i++) {
        const x = p.view[0] + (p.view[1] - p.view[0]) * i / n;
        pts.push(`${sx(x).toFixed(1)},${cy(p.f(x)).toFixed(1)}`);
      }
      svgEl("polyline", { points: pts.join(" "), fill: "none", stroke: COL.curve, "stroke-width": 3 }, pg);
    } else {
      svgEl("text", { x: (L + Rm) / 2, y: bot - 14, "text-anchor": "middle", class: "sm" }, svg).textContent =
        "curve hidden: the method sees only its own points";
    }
    // the parabola: the one used by the last step (Brent) or the next one (parabolic)
    let par: Pt[] | undefined;
    let parLabel = "";
    if (run.method === "parabolic" && run.evals.length && !run.done) { par = run.inner(); parLabel = "next: vertex"; }
    else if (run.rows.length && run.rows[run.rows.length - 1].parab) { par = run.rows[run.rows.length - 1].parab; parLabel = "vertex"; }
    if (par && new Set(par.map(q => q.x)).size === 3) {
      const P = parabola(par), n = 300, pts: string[] = [];
      for (let i = 0; i <= n; i++) {
        const x = p.view[0] + (p.view[1] - p.view[0]) * i / n;
        pts.push(`${sx(x).toFixed(1)},${cy(P(x)).toFixed(1)}`);
      }
      svgEl("polyline", { points: pts.join(" "), fill: "none", stroke: COL.par, "stroke-width": 2.5, "stroke-dasharray": "8 6" }, pg);
      const v = vertex(par);
      if (Number.isFinite(v) && v > p.view[0] && v < p.view[1]) {
        svgEl("line", { x1: sx(v), x2: sx(v), y1: top, y2: bot, stroke: COL.par, "stroke-width": 2, "stroke-dasharray": "3 4" }, pg);
        const t = svgEl("text", { x: sx(v) + 5, y: bot - 10, class: "sm" }, pg);
        t.style.fill = COL.par; t.textContent = parLabel;
      }
      for (const q of par) svgEl("circle", { cx: sx(q.x), cy: cy(q.f), r: 9, fill: "none", stroke: COL.par, "stroke-width": 2.5 }, pg);
    }
    // all evaluated points
    for (const e of run.evals) {
      if (!Number.isFinite(e.f)) continue;
      svgEl("circle", { cx: sx(e.x), cy: cy(e.f), r: 4.5, fill: COL.ink, "fill-opacity": 0.55 }, pg);
    }
    // interior points
    const inner = run.inner();
    const names = run.method === "golden" ? ["x₁", "x₂"] : run.method === "brent" ? ["x", "w", "v"] : [];
    inner.forEach((q, i) => {
      if (!Number.isFinite(q.f) || run.method === "parabolic") return;
      const c = run.method === "golden" ? COL.gold : i === 0 ? COL.par : COL.gold;
      svgEl("line", { x1: sx(q.x), x2: sx(q.x), y1: cy(q.f), y2: bot, stroke: c, "stroke-width": 1.5, "stroke-dasharray": "2 4" }, pg);
      svgEl("circle", { cx: sx(q.x), cy: cy(q.f), r: i === 0 && run.method === "brent" ? 8 : 7, fill: c, stroke: "#fff", "stroke-width": 1.5 }, pg);
      if (run.method === "brent" && i > 0 && (q.x === inner[0].x || (i === 2 && q.x === inner[1].x))) return;
      const t = svgEl("text", { x: sx(q.x), y: cy(q.f) - 14 - (run.method === "brent" ? 18 * i : 0), "text-anchor": "middle", "font-weight": 700 }, pg);
      t.style.fill = c; t.textContent = names[i];
    });
    // the newest point
    if (run.evals.length) {
      const e = run.evals[run.evals.length - 1];
      if (Number.isFinite(e.f) && e.x >= p.view[0] && e.x <= p.view[1])
        svgEl("circle", { cx: sx(e.x), cy: cy(e.f), r: 12, fill: "none", stroke: COL.ink, "stroke-width": 2 }, pg);
      else if (e.x < p.view[0] || e.x > p.view[1]) {
        const t = svgEl("text", { x: e.x < p.view[0] ? L + 6 : Rm - 6, y: top + 70, "text-anchor": e.x < p.view[0] ? "start" : "end", "font-weight": 700 }, svg);
        t.style.fill = COL.par; t.textContent = `${e.x < p.view[0] ? "◀ " : ""}new point at ${fx(e.x, 2)}${e.x > p.view[1] ? " ▶" : ""}`;
      }
    }
    // multistart results
    if (st.multi && st.multi.pts.length) {
      const best = st.multi.pts.reduce((u, v) => (v.f < u.f ? v : u));
      for (const q of st.multi.pts) {
        const X = sx(q.x), Y = cy(q.f), s = 10;
        svgEl("path", { d: `M${X},${Y - s} L${X + s},${Y} L${X},${Y + s} L${X - s},${Y} Z`, fill: q === best ? COL.grn : "#fff", stroke: COL.grn, "stroke-width": 2.5 }, pg);
      }
      const t = svgEl("text", { x: sx(best.x), y: cy(best.f) + 32, "text-anchor": "middle", "font-weight": 700 }, pg);
      t.style.fill = COL.grn; t.textContent = `best ${fx(best.x, 2)}`;
    }
    // drag handles for the initial bracket
    svg.appendChild(hLayer);
    const hx = [st.a0, st.b0];
    for (let i = 0; i < 2; i++) {
      handles[i].setAttribute("cx", String(sx(hx[i]))); handles[i].setAttribute("cy", String(bot + 17));
      handles[i].style.display = "";
      hLabels[i].setAttribute("x", String(sx(hx[i]))); hLabels[i].setAttribute("y", String(bot + 23));
      hLabels[i].textContent = i === 0 ? "a" : "b";
      hLabels[i].style.fontSize = "16px";
    }
  }

  // ---------------------------------------------------------------- (a') the flat bottom
  function drawZoom(p: Problem, run: Run): void {
    svg.replaceChildren();
    handles.forEach(h => (h.style.display = "none"));
    const Wd = Math.pow(10, st.lw);
    const xs = p.xstar, f0 = p.f(xs);
    const u = ulp(f0);
    const ymax = Math.max(0.5 * p.f2 * Wd * Wd, 6 * u) * 1.15;
    const ymin = -0.25 * ymax;
    const topZ = top + 26;
    const sx = scale(-Wd, Wd, L, Rm), sy = scale(ymin, ymax, bot, topZ);
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svg).textContent = `(a) the bottom of ${p.fname}, zoomed in`;
    const yT = niceTicks(ymin, ymax, 5);
    for (const v of yT) {
      svgEl("line", { x1: L, x2: Rm, y1: sy(v), y2: sy(v), class: v === 0 ? "axis" : "grid" }, svg);
      svgEl("text", { x: L - 6, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = sciTick(v);
    }
    svgEl("line", { x1: L, x2: Rm, y1: bot, y2: bot, class: "axis" }, svg);
    for (const v of [-1, -0.5, 0, 0.5, 1]) {
      svgEl("line", { x1: sx(v * Wd), x2: sx(v * Wd), y1: bot, y2: bot + 5, class: "axis" }, svg);
      svgEl("text", { x: sx(v * Wd), y: bot + 24, "text-anchor": v === -1 ? "start" : v === 1 ? "end" : "middle", class: "sm" }, svg).textContent = sciTick(v * Wd);
    }
    svgEl("text", { x: (L + Rm) / 2, y: H - 14, "text-anchor": "middle" }, svg).textContent =
      `distance from the minimiser, ${p.xname.split(" (")[0].replace(/^.* /, "")} − ${fx(xs, 6)}`;
    const yl = svgEl("text", { x: 16, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(top + bot) / 2})` }, svg);
    yl.textContent = "f(x) − f(x̂)";

    const clip = svgEl("clipPath", { id: `gs-clipz-${uid}` }, svg);
    svgEl("rect", { x: L, y: topZ - 6, width: Rm - L, height: bot - topZ + 6 }, clip);
    const pg = svgEl("g", { "clip-path": `url(#gs-clipz-${uid})` }, svg);
    // predicted flat band
    const fw = flatWidth(p);
    if (fw < Wd * 1.2) {
      svgEl("rect", { x: sx(-fw), y: topZ, width: sx(fw) - sx(-fw), height: bot - topZ, fill: COL.grn, "fill-opacity": 0.13 }, pg);
      const t = svgEl("text", { x: sx(0), y: topZ + 22, "text-anchor": "middle", class: "sm" }, pg);
      t.style.fill = COL.grn; t.textContent = `±√(2ε|f|/f″) = ±${fmtSci(fw)}`;
    }
    // parabola 1/2 f'' d^2
    const n = 300, pp: string[] = [];
    for (let i = 0; i <= n; i++) { const d = -Wd + 2 * Wd * i / n; pp.push(`${sx(d).toFixed(1)},${sy(0.5 * p.f2 * d * d).toFixed(1)}`); }
    svgEl("polyline", { points: pp.join(" "), fill: "none", stroke: COL.par, "stroke-width": 2.5 }, pg);
    // computed values
    const m = 401;
    for (let i = 0; i < m; i++) {
      const d = -Wd + 2 * Wd * i / (m - 1);
      const v = p.f(xs + d) - f0;
      svgEl("circle", { cx: sx(d), cy: sy(Math.max(ymin, Math.min(ymax, v))), r: 2.6, fill: COL.curve }, pg);
    }
    // the current estimate of the run
    if (run.rows.length) {
      const e = run.rows[run.rows.length - 1].best - xs;
      if (Math.abs(e) <= Wd) {
        svgEl("line", { x1: sx(e), x2: sx(e), y1: topZ + 30, y2: bot, stroke: MCOL[run.method], "stroke-width": 3 }, pg);
        const t = svgEl("text", { x: sx(e) + (e > 0 ? -6 : 6), y: topZ + 50, "text-anchor": e > 0 ? "end" : "start", "font-weight": 700 }, pg);
        t.style.fill = MCOL[run.method]; t.textContent = "current estimate";
      }
    }
    // legend, under the title
    svgEl("circle", { cx: L + 8, cy: top + 6, r: 4.5, fill: COL.curve }, svg);
    svgEl("text", { x: L + 20, y: top + 12, class: "sm" }, svg).textContent = "f computed in doubles";
    svgEl("line", { x1: L + 230, x2: L + 256, y1: top + 6, y2: top + 6, stroke: COL.par, "stroke-width": 3 }, svg);
    svgEl("text", { x: L + 262, y: top + 12, class: "sm" }, svg).textContent = "½ f″(x̂) δ² (exact)";
  }

  // ---------------------------------------------------------------- (b) convergence
  function drawLog(p: Problem, run: Run): void {
    svgL.replaceChildren();
    const lL = 66, lR = WL - 12, lT = 34, lB = HL - 44;
    const nmax = Math.max(30, Math.min(100, Math.ceil((run.evals.length + 2) / 10) * 10));
    const Y0 = -16, Y1 = 2;
    const sx = scale(0, nmax, lL, lR), sy = scale(Y0, Y1, lB, lT);
    const cy = (v: number) => sy(Math.max(Y0, Math.min(Y1, Math.log10(Math.max(v, 1e-300)))));
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgL).textContent = "(b) width and error (log scale)";
    for (let k = Y0; k <= Y1; k += 4) {
      svgEl("line", { x1: lL, x2: lR, y1: sy(k), y2: sy(k), class: "grid" }, svgL);
      svgEl("text", { x: lL - 6, y: sy(k) + 6, "text-anchor": "end", class: "sm" }, svgL).textContent = pow10(k);
    }
    svgEl("line", { x1: lL, x2: lR, y1: lB, y2: lB, class: "axis" }, svgL);
    const step = nmax <= 40 ? 10 : 20;
    for (let k = 0; k <= nmax; k += step) {
      svgEl("text", { x: sx(k), y: lB + 22, "text-anchor": "middle", class: "sm" }, svgL).textContent = String(k);
    }
    svgEl("text", { x: (lL + lR) / 2, y: HL - 4, "text-anchor": "middle", class: "sm" }, svgL).textContent = "evaluations of f";
    const clip = svgEl("clipPath", { id: `gs-clipl-${uid}` }, svgL);
    svgEl("rect", { x: lL, y: lT - 4, width: lR - lL, height: lB - lT + 8 }, clip);
    const g = svgEl("g", { "clip-path": `url(#gs-clipl-${uid})` }, svgL);
    // reference 0.618^n and the flat-bottom floor
    const L0 = st.b0 - st.a0;
    svgEl("line", { x1: sx(1), y1: cy(L0), x2: sx(nmax), y2: cy(L0 * Math.pow(G, nmax - 1)), stroke: COL.grey, "stroke-width": 2, "stroke-dasharray": "8 6" }, g);
    const fw = flatWidth(p);
    svgEl("line", { x1: lL, x2: lR, y1: cy(fw), y2: cy(fw), stroke: COL.grn, "stroke-width": 2, "stroke-dasharray": "3 4" }, g);
    // series
    const series = (vals: [number, number][], c: string, dash: string) => {
      const v = vals.filter(q => Number.isFinite(q[1]) && q[1] >= 0);
      if (!v.length) return;
      svgEl("polyline", { points: v.map(([n, y]) => `${sx(n)},${cy(y)}`).join(" "), fill: "none", stroke: c, "stroke-width": 2.5, "stroke-dasharray": dash }, g);
      for (const [n, y] of v) svgEl("circle", { cx: sx(n), cy: cy(y), r: 3.5, fill: c }, g);
    };
    const rows = run.rows;
    if (run.method !== "parabolic") series(rows.map(r => [r.n, r.b - r.a]), COL.width, "");
    series(rows.map(r => [r.n, Math.abs(r.best - p.xstar)]), MCOL[run.method], "");
    // legend (HTML, below the plot)
    legend.replaceChildren();
    const lg = (txt: string, c: string, dash: boolean) => {
      const s = htmlEl("span", "", "", legend);
      const sw = htmlEl("span", "gs-sw", "", s);
      sw.style.borderTop = `3px ${dash ? "dashed" : "solid"} ${c}`;
      s.append(txt);
    };
    if (run.method !== "parabolic") lg("width b − a", COL.width, false);
    lg("error |x − x̂|", MCOL[run.method], false);
    lg("0.618ⁿ(b₀ − a₀)", COL.grey, true);
    lg("flat bottom", COL.grn, true);
  }

  // ---------------------------------------------------------------- table
  function drawTable(p: Problem, run: Run): void {
    const dec = p.key === "circadian" ? 6 : 8;
    tbox.replaceChildren();
    const tb = htmlEl("table", "", "", tbox);
    const hr = htmlEl("tr", "", "", htmlEl("thead", "", "", tb));
    const heads = run.method === "parabolic" ? ["n", "step", "new x", "f(new x)", "note"] : ["n", "step", "new x", "f(new x)", "width b − a"];
    for (const h of heads) htmlEl("th", "", h, hr);
    const body = htmlEl("tbody", "", "", tb);
    if (!run.rows.length) {
      const tr = htmlEl("tr", "", "", body);
      const td = htmlEl("td", "k", "press Step", tr); td.colSpan = 5;
    }
    run.rows.forEach((r, i) => {
      const tr = htmlEl("tr", i === run.rows.length - 1 ? "last" : "", "", body);
      htmlEl("td", "", String(r.n), tr);
      htmlEl("td", "k", r.kind === "parabolic" ? "parab." : r.kind, tr);
      htmlEl("td", "", fx(r.x, dec), tr);
      htmlEl("td", "", Number.isFinite(r.f) ? (Math.abs(r.f) < 1e-3 ? r.f.toExponential(5) : r.f.toPrecision(7)) : String(r.f), tr);
      htmlEl("td", run.method === "parabolic" ? "k" : "", run.method === "parabolic" ? (r.note ?? "") : fmtSci(r.b - r.a, 3), tr);
    });
    tbox.scrollTop = tbox.scrollHeight;
  }

  // ---------------------------------------------------------------- readout
  function drawReadout(p: Problem, run: Run): void {
    ro.replaceChildren();
    const line = (parts: [string, string, string?][]) => {
      const d = htmlEl("div", "", "", ro);
      parts.forEach(([k, v, cls], i) => {
        if (i) d.append("   ");
        if (k) htmlEl("b", "", k + " ", d);
        const s = htmlEl("span", cls ?? "", v, d);
        void s;
      });
    };
    if (st.zoom) {
      const f0 = p.f(p.xstar), Wd = Math.pow(10, st.lw);
      let same = 0, lo = Infinity, hi = -Infinity;
      for (let i = 0; i < 401; i++) {
        const d = -Wd + 2 * Wd * i / 400;
        if (p.f(p.xstar + d) === f0) { same++; lo = Math.min(lo, d); hi = Math.max(hi, d); }
      }
      line([["x̂ =", p.xstar.toPrecision(12)], ["f(x̂) =", f0.toPrecision(10)], ["f″(x̂) =", p.f2.toPrecision(4)],
        ["gap between doubles near f(x̂):", fmtSci(ulp(f0), 2)], ["predicted flat half-width:", fmtSci(flatWidth(p), 3)]]);
      line([["zoom ±", fmtSci(Wd, 1)], ["points with f(x) == f(x̂):", `${same} of 401` + (same > 1 ? `, from ${fmtSci(lo)} to ${fmtSci(hi)}` : ""), same > 1 ? "gs-warn" : ""]]);
      return;
    }
    const nE = run.evals.length;
    const br = run.bracket();
    const last = run.rows[run.rows.length - 1];
    const est = last ? last.best : NaN;
    line([
      ["evaluations:", String(nE)],
      ["bracket:", br ? `[${fx(br[0], p.digits)}, ${fx(br[1], p.digits)}], width ${fmtSci(br[1] - br[0], 3)}` : "none (parabolic keeps no bracket)"],
      ["estimate x =", Number.isFinite(est) ? est.toPrecision(10) : "–"],
      ["error |x − x̂| =", Number.isFinite(est) ? fmtSci(Math.abs(est - p.xstar), 2) : "–"],
    ]);
    const parts: [string, string, string?][] = [];
    if (run.done) parts.push(["stopped:", run.why, /domain|away|line/.test(run.why) ? "gs-warn" : ""]);
    else parts.push(["next:", run.next()]);
    if (p.other !== undefined && Number.isFinite(est) && Math.abs(est - p.other) < 1 && nE > 3)
      parts.push(["", `a local minimum: RSS ${p.f(p.other).toFixed(2)}, but RSS(${fx(p.xstar, 2)}) = ${p.fstar.toFixed(2)}`, "gs-warn"]);
    if (st.multi && st.multi.pts.length) {
      const best = st.multi.pts.reduce((u, v) => (v.f < u.f ? v : u));
      parts.push([`multistart (Brent on 5 sub-intervals, ${st.multi.evals} evals), interior minima:`,
        st.multi.pts.map(q => `${fx(q.x, 2)} (f ${q.f.toPrecision(4)})`).join(", ") + ` → best ${fx(best.x, 3)}`, "gs-ok"]);
    }
    line(parts);
  }

  setProblem(init.prob, true);
  draw();
}

// ---- axis helpers ----
function niceTicks(lo: number, hi: number, n: number): number[] {
  const span = hi - lo, raw = span / n;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const st = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => s >= raw) ?? 10 * mag;
  const out: number[] = [];
  for (let v = Math.ceil(lo / st - 1e-9) * st; v <= hi + 1e-9 * span; v += st) out.push(Math.abs(v) < 1e-12 * span ? 0 : v);
  return out;
}
function tickFmt(v: number): string {
  if (v === 0) return "0";
  return String(+v.toPrecision(4));
}
function sciTick(v: number): string {
  if (v === 0) return "0";
  const [m, e] = v.toExponential(1).split("e");
  return `${+m}e${Number(e)}`;
}

autoMount(NAME, mount);

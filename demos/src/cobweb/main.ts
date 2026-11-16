// Cobweb explorer for the bisection-fixedpoint unit. A discrete-time population model
// N_{t+1} = g(N_t) (Beverton–Holt, Ricker, or the lecture's harvested Ricker stock,
// K = 100) is a fixed-point iteration.
//  Tab "Cobweb": (a) the cobweb plot of g with the diagonal; click or drag in it to move the
//    starting population N0; steady states marked filled (attracting, |g'| < 1) or open
//    (repelling). (b) the same iterates as a time series. Readout: g'(N*) at every steady
//    state, the ratio of the last two steps, and the long-run behaviour (steady state,
//    k-cycle, collapse, or no repeat = chaos).
//  Tab "Bisection": (a) f(N) = g(N) - N with a draggable starting bracket [a, b]; step
//    through bisection. (b) error bound and actual error of the midpoint against the step,
//    log scale. With H near 14.26 the two roots merge into a double root: no sign change.
// Auto-mounts on <div data-demo="cobweb">. Optional attributes:
//   data-model="harvest|ricker|bh"  data-r="0.5"  data-h="10"  data-n0="30"  data-years="30"
//   data-mode="cobweb|bisect"  data-bisect="1|0" (0 hides the bisection tab)
//   data-a="50" data-b="100" (starting bracket)
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, svgPoint, autoMount, Slider } from "../lib/ui";
import { Model, Params, g, f, range, steadyStates, maxCatch, orbit, longRun, bisect, Steady } from "./math";

const NAME = "cobweb";
// Okabe–Ito, as in make-bisection-fixedpoint-figs.R
const COL = {
  map: "#0072B2", path: "#D55E00", stable: "#009E73", unstable: "#E69F00", gg: "#CC79A7",
  diag: "#999999", sky: "#56B4E9", ink: "#222222", grey: "#777777", root: "#D55E00", brk: "#E69F00",
};
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .cw-row { display: flex; flex-wrap: wrap; gap: .35em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .cw-group { display: flex; flex-wrap: wrap; gap: .35em; align-items: center; }
[data-demo="${NAME}"] .cw-row .qb-slider { flex: 1 1 300px; min-width: 260px; max-width: 420px; grid-template-columns: auto 1fr 3.4em; }
[data-demo="${NAME}"] .cw-plots { display: flex; flex-wrap: wrap; gap: 8px; width: 100%; }
[data-demo="${NAME}"] .cw-plots svg { min-width: 0; height: auto; touch-action: none; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 19px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 17px; fill: #555; }
[data-demo="${NAME}"] .cw-read { display: flex; flex-wrap: wrap; gap: .15em 1.3em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .cw-read b { font-weight: 600; }
[data-demo="${NAME}"] .cw-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .cw-sep { width: 1px; align-self: stretch; background: #ccc; }
[data-demo="${NAME}"] .qb-btn.cw-tab { font-size: .95em; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; gap: .3em; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; padding: .3em .5em; }
.reveal [data-demo="${NAME}"] .qb-btn { font-size: .78em; padding: .18em .55em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- formatting ----
const sub: Record<string, string> = { "0": "₀", "1": "₁", "2": "₂", "3": "₃", "4": "₄", "5": "₅", "6": "₆", "7": "₇", "8": "₈", "9": "₉" };
const subs = (k: number): string => String(k).split("").map(c => sub[c] ?? c).join("");
const supd: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const pow10 = (k: number): string => k === 0 ? "1" : "10" + String(k).split("").map(c => supd[c]).join("");
const minus = (s: string): string => s.replace(/^-/, "−");
const fx = (v: number, d = 2): string => minus(v.toFixed(d));
function sci(v: number, d = 2): string {
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 0.01 && a < 1e5) return minus(v.toPrecision(d + 1));
  const [m, e] = v.toExponential(d).split("e");
  return minus(`${m}e${Number(e)}`);
}

const MODELS: { key: Model; label: string }[] = [
  { key: "bh", label: "Beverton–Holt" },
  { key: "ricker", label: "Ricker" },
  { key: "harvest", label: "harvested Ricker" },
];
const R_RANGE: Record<Model, [number, number, number]> = {
  bh: [0.5, 5, 0.1], ricker: [0.1, 3.5, 0.01], harvest: [0.1, 3, 0.01],
};

type Mode = "cobweb" | "bisect";
interface State { mode: Mode; model: Model; r: number; H: number; N0: number; n: number; gg: boolean; a: number; b: number; k: number }

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const uid = ++mountCount;
  const num = (s: string | undefined, def: number) => (s !== undefined && s !== "" && Number.isFinite(Number(s)) ? Number(s) : def);
  const model0: Model = ds.model === "bh" || ds.model === "ricker" ? ds.model : "harvest";
  const withBisect = ds.bisect !== "0";
  const init: State = {
    mode: withBisect && ds.mode === "bisect" ? "bisect" : "cobweb",
    model: model0,
    r: num(ds.r, model0 === "bh" ? 2 : 0.5),
    H: num(ds.h, model0 === "harvest" ? 10 : 0),
    N0: num(ds.n0, model0 === "harvest" ? 30 : 20),
    n: Math.round(num(ds.years, 30)),
    gg: false,
    a: num(ds.a, 50), b: num(ds.b, 100), k: 0,
  };
  const st: State = { ...init };
  const P = (): Params => ({ model: st.model, r: st.r, H: st.model === "harvest" ? st.H : 0 });

  // ---------------- controls ----------------
  const row1 = htmlEl("div", "cw-row", "", host);
  let tabBtns: HTMLButtonElement[] = [];
  if (withBisect) {
    const tg = htmlEl("div", "cw-group", "", row1);
    tabBtns = presets(tg, [
      { label: "Iterate: cobweb", apply: () => { st.mode = "cobweb"; rebuild(); } },
      { label: "Bisection on f = g − N", apply: () => { st.mode = "bisect"; rebuild(); } },
    ]);
    tabBtns.forEach(b => b.classList.add("cw-tab"));
    htmlEl("span", "cw-sep", "", row1);
  }
  const mg = htmlEl("div", "cw-group", "", row1);
  htmlEl("span", "", "model:", mg);
  const modelBtns = presets(mg, MODELS.map(m => ({ label: m.label, apply: () => setModel(m.key) })));
  button(row1, "reset", () => { Object.assign(st, init); rebuild(); });

  const row2 = htmlEl("div", "cw-row", "", host);
  const row3 = htmlEl("div", "cw-row", "", host);
  const row4 = htmlEl("div", "cw-row", "", host);

  const plots = htmlEl("div", "cw-plots", "", host);
  const svgA = svgEl("svg", { role: "img" }, plots);
  const svgB = svgEl("svg", { role: "img" }, plots);
  const ro = htmlEl("div", "qb-readout cw-read", "", host);

  let rS: Slider | null = null, hS: Slider | null = null, n0S: Slider | null = null, nS: Slider | null = null;
  let presetBtns: HTMLButtonElement[] = [];

  function setModel(m: Model): void {
    st.model = m;
    st.r = m === "bh" ? 2 : 0.5;
    st.H = m === "harvest" ? 10 : 0;
    st.N0 = m === "harvest" ? 30 : 20;
    st.k = 0;
    if (m === "harvest") { st.a = 50; st.b = 100; } else { st.a = 50; st.b = 150; }
    rebuild();
  }

  /** Rebuild the model- and mode-dependent controls, then draw. */
  function rebuild(): void {
    tabBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.mode === "cobweb")));
    modelBtns.forEach((b, i) => b.classList.toggle("active", MODELS[i].key === st.model));
    row2.replaceChildren(); row3.replaceChildren(); row4.replaceChildren();
    const [r0, r1, rstep] = R_RANGE[st.model];
    rS = slider(row2, st.model === "bh" ? "growth R" : "growth r", r0, r1, rstep, st.r, v => { st.r = v; st.k = 0; draw(); });
    hS = null;
    if (st.model === "harvest") hS = slider(row2, "catch H", 0, 16, 0.01, st.H, v => { st.H = v; st.k = 0; draw(); });
    n0S = nS = null;
    if (st.mode === "cobweb") {
      n0S = slider(row3, "start N₀", 0, 150, 0.1, st.N0, v => { st.N0 = v; draw(); }, v => v.toFixed(1));
      nS = slider(row3, "years", 1, 200, 1, st.n, v => { st.n = v; draw(); });
      const sg = htmlEl("div", "cw-group", "", row3);
      button(sg, "−1 year", () => { st.n = Math.max(1, st.n - 1); nS?.set(st.n); draw(); });
      button(sg, "+1 year", () => { st.n = Math.min(200, st.n + 1); nS?.set(st.n); draw(); });
      htmlEl("span", "", "presets:", row4);
      const pg = htmlEl("div", "cw-group", "", row4);
      const set = (o: Partial<State>) => () => { Object.assign(st, o); syncSliders(); draw(); };
      const items: { label: string; apply: () => void }[] =
        st.model === "bh" ? [
          { label: "R = 2, N₀ = 20", apply: set({ r: 2, N0: 20, n: 10 }) },
          { label: "R = 3, N₀ = 50", apply: set({ r: 3, N0: 50, n: 10 }) },
          { label: "R = 0.8", apply: set({ r: 0.8, N0: 50, n: 30 }) },
        ] : st.model === "ricker" ? [0.5, 1.5, 1.9, 2.3, 2.6, 3.0, 3.12].map(r => ({
          label: `r = ${r}`, apply: set({ r, N0: 20, n: r >= 2.6 ? 60 : 30 }),
        })) : [
          { label: "N₀ = 30", apply: set({ r: 0.5, H: 10, N0: 30, n: 30 }) },
          { label: "N₀ = 21", apply: set({ r: 0.5, H: 10, N0: 21, n: 40 }) },
          { label: "N₀ = 20", apply: set({ r: 0.5, H: 10, N0: 20, n: 40 }) },
          { label: "H = 14", apply: set({ r: 0.5, H: 14, N0: 60, n: 60 }) },
          { label: "H = 14.26", apply: set({ r: 0.5, H: 14.26, N0: 60, n: 100 }) },
          { label: "H = 15", apply: set({ r: 0.5, H: 15, N0: 60, n: 40 }) },
          { label: "r = 2.3", apply: set({ r: 2.3, H: 10, N0: 30, n: 40 }) },
        ];
      presetBtns = presets(pg, items);
      checkbox(row4, "show g(g(N))", st.gg, v => { st.gg = v; draw(); });
    } else {
      const bg = htmlEl("div", "cw-group", "", row3);
      htmlEl("span", "", "bisection:", bg);
      button(bg, "step", () => { st.k = Math.min(60, st.k + 1); draw(); });
      button(bg, "+10 steps", () => { st.k = Math.min(60, st.k + 10); draw(); });
      button(bg, "restart", () => { st.k = 0; draw(); });
      htmlEl("span", "qb-note", "drag in the plot to move the bracket ends a, b", row3);
      htmlEl("span", "", "bracket:", row4);
      const pg = htmlEl("div", "cw-group", "", row4);
      const br = (a: number, b: number) => () => { st.a = a; st.b = b; st.k = 0; draw(); };
      const hi = range(P());
      const items = st.model === "harvest" ? [
        { label: "[0, 50]", apply: br(0, 50) }, { label: "[50, 100]", apply: br(50, 100) },
        { label: "[0, 100]", apply: br(0, 100) }, { label: "[40, 47]", apply: br(40, 47) },
      ] : [
        { label: "[50, 150]", apply: br(50, Math.min(150, hi)) }, { label: "[10, 200]", apply: br(10, 200) },
        { label: "[0, 150]", apply: br(0, 150) },
      ];
      presetBtns = presets(pg, items);
      if (st.model === "harvest") {
        const hg = htmlEl("div", "cw-group", "", row4);
        htmlEl("span", "", "catch:", hg);
        presets(hg, [10, 14, 14.26, 15].map(H => ({ label: `H = ${H}`, apply: () => { st.H = H; st.r = 0.5; st.k = 0; syncSliders(); draw(); } })));
      }
    }
    draw();
  }
  function syncSliders(): void {
    rS?.set(st.r); hS?.set(st.H); n0S?.set(st.N0); nS?.set(st.n);
  }

  // ---------------- geometry ----------------
  const H = 340, top = 34, bot = H - 48;
  let WA = 360, WB = 740;
  let ax = scale(0, 150, 62, 62 + (bot - top)), ay = scale(0, 150, bot, top);
  let fy = scale(-1, 1, bot, top);

  function setSizes(wa: number, wb: number): void {
    WA = wa; WB = wb;
    const WT = WA + WB;
    svgA.setAttribute("viewBox", `0 0 ${WA} ${H}`);
    svgB.setAttribute("viewBox", `0 0 ${WB} ${H}`);
    svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 4px)`;
    svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 4px)`;
  }

  // ---------------- pointer ----------------
  let dragEnd: "a" | "b" | null = null;
  function pointer(ev: PointerEvent, start: boolean): void {
    const [px] = svgPoint(svgA, ev);
    const x = ax.invert(px), hi = range(P());
    if (st.mode === "cobweb") {
      st.N0 = Math.round(Math.min(hi, Math.max(0, x)) * 10) / 10;
      if (st.N0 > 150) st.N0 = 150;
      n0S?.set(st.N0);
    } else {
      const v = Math.round(Math.min(hi, Math.max(0, x)) * 2) / 2;
      if (start) dragEnd = Math.abs(v - st.a) <= Math.abs(v - st.b) ? "a" : "b";
      if (dragEnd === "a") st.a = v; else st.b = v;
      if (st.a > st.b) { [st.a, st.b] = [st.b, st.a]; dragEnd = dragEnd === "a" ? "b" : "a"; }
      st.k = 0;
    }
    presetBtns.forEach(b => b.classList.remove("active"));
    draw();
  }
  svgA.addEventListener("pointerdown", ev => {
    const [, py] = svgPoint(svgA, ev);
    if (py < top - 10 || py > bot + 30) return;
    ev.preventDefault();
    svgA.setPointerCapture(ev.pointerId);
    pointer(ev, true);
    const move = (e: PointerEvent) => pointer(e, false);
    const up = () => {
      svgA.removeEventListener("pointermove", move);
      svgA.removeEventListener("pointerup", up);
      svgA.removeEventListener("pointercancel", up);
      dragEnd = null;
    };
    svgA.addEventListener("pointermove", move);
    svgA.addEventListener("pointerup", up);
    svgA.addEventListener("pointercancel", up);
  });

  // ---------------- drawing helpers ----------------
  const T = (parent: Element, x: number, y: number, s: string, attrs: Record<string, string | number> = {}, fill?: string) => {
    const t = svgEl("text", { x, y, ...attrs }, parent);
    t.textContent = s;
    if (fill) t.style.fill = fill;
    return t;
  };
  function clipRect(svg: SVGSVGElement, id: string, x0: number, x1: number): SVGGElement {
    const cp = svgEl("clipPath", { id }, svg);
    svgEl("rect", { x: x0, y: top - 2, width: x1 - x0, height: bot - top + 4 }, cp);
    return svgEl("g", { "clip-path": `url(#${id})` }, svg);
  }
  function ticksFor(hi: number): number[] {
    const step = hi <= 150 ? 50 : hi <= 300 ? 100 : 200;
    const t: number[] = [];
    for (let v = 0; v <= hi + 1e-9; v += step) t.push(v);
    return t;
  }
  const yLabel = (svg: SVGSVGElement, x: number, s: string) =>
    T(svg, x, (top + bot) / 2, s, { "text-anchor": "middle", transform: `rotate(-90 ${x} ${(top + bot) / 2})` });
  const classify = (d: number): string =>
    Math.abs(d) < 1 ? `attracting, ${d >= 0 ? "monotone" : "alternating"}` : `repelling${d < 0 ? ", alternating" : ""}`;
  function marker(parent: Element, x: number, y: number, s: Steady, rad = 7): void {
    const stable = Math.abs(s.d) < 1;
    svgEl("circle", { cx: x, cy: y, r: rad, fill: stable ? COL.stable : "#fff", stroke: stable ? "#fff" : COL.unstable, "stroke-width": stable ? 1.5 : 3 }, parent);
  }

  // ---------------- draw ----------------
  function draw(): void {
    if (st.mode === "cobweb") drawCobweb(); else drawBisect();
  }

  function drawCobweb(): void {
    setSizes(360, 740);
    const p = P(), hi = range(p);
    const side = bot - top;
    ax = scale(0, hi, 62, 62 + side); ay = scale(0, hi, bot, top);
    if (st.N0 > hi) st.N0 = hi;
    const o = orbit(p, st.N0, st.n);
    const ss = steadyStates(p, hi);
    svgA.setAttribute("aria-label", "Cobweb plot of the population map");
    svgB.setAttribute("aria-label", "Population against year");

    // ---- (a) cobweb ----
    svgA.replaceChildren();
    T(svgA, 8, 22, "(a) cobweb: drag to set N₀", { class: "ttl" });
    for (const v of ticksFor(hi)) {
      svgEl("line", { x1: ax(0), x2: ax(hi), y1: ay(v), y2: ay(v), class: v === 0 ? "axis" : "grid" }, svgA);
      svgEl("line", { x1: ax(v), x2: ax(v), y1: top, y2: bot, class: v === 0 ? "axis" : "grid" }, svgA);
      T(svgA, ax(0) - 6, ay(v) + 6, String(v), { "text-anchor": "end", class: "sm" });
      T(svgA, ax(v), bot + 20, String(v), { "text-anchor": "middle", class: "sm" });
    }
    T(svgA, ax(hi / 2), H - 6, "Nₜ", { "text-anchor": "middle" });
    yLabel(svgA, 16, "Nₜ₊₁");
    const pg = clipRect(svgA, `cw-clipA-${uid}`, ax(0), ax(hi));
    svgEl("line", { x1: ax(0), y1: ay(0), x2: ax(hi), y2: ay(hi), stroke: COL.diag, "stroke-width": 2 }, pg);
    const curve = (fn: (N: number) => number, stroke: string, dash = "", w = 3) => {
      const pts: string[] = [];
      for (let i = 0; i <= 400; i++) {
        const N = hi * i / 400, y = Math.max(-hi * 0.1, Math.min(hi * 1.2, fn(N)));
        pts.push(`${ax(N).toFixed(1)},${ay(y).toFixed(1)}`);
      }
      svgEl("polyline", { points: pts.join(" "), fill: "none", stroke, "stroke-width": w, "stroke-dasharray": dash }, pg);
    };
    if (st.gg) curve(N => { const v = g(p, N); return v > 0 ? g(p, v) : -1; }, COL.gg, "7 5", 2.5);
    curve(N => g(p, N), COL.map, "", 3);
    // cobweb path
    const xs = o.x, path: string[] = [`${ax(xs[0])},${ay(0)}`];
    for (let t = 0; t < xs.length - 1; t++) {
      path.push(`${ax(xs[t])},${ay(xs[t + 1])}`);
      if (o.collapse === t + 1) break;
      path.push(`${ax(xs[t + 1])},${ay(xs[t + 1])}`);
    }
    svgEl("polyline", { points: path.join(" "), fill: "none", stroke: COL.path, "stroke-width": 2, "stroke-opacity": 0.85, "stroke-linejoin": "round" }, pg);
    for (const s of ss) marker(pg, ax(s.x), ay(s.x), s);
    if (o.collapse >= 0) {
      const xc = ax(xs[o.collapse - 1]);
      T(svgA, xc + 4, ay(0) - 8, "✕ collapse", { "font-weight": 600 }, COL.path);
    }
    // N0 handle
    svgEl("circle", { cx: ax(st.N0), cy: ay(0), r: 10, fill: "#fff", stroke: COL.path, "stroke-width": 3, class: "handle" }, svgA);
    T(svgA, ax(hi) - 4, ay(hi * 0.86), "y = x", { "text-anchor": "end", class: "sm" });

    // ---- (b) time series ----
    svgB.replaceChildren();
    T(svgB, 8, 22, `(b) population Nₜ, years 0–${st.n}`, { class: "ttl" });
    const L = 62, R = WB - 14;
    const tx = scale(0, st.n, L, R), ty = ay;
    for (const v of ticksFor(hi)) {
      svgEl("line", { x1: L, x2: R, y1: ty(v), y2: ty(v), class: v === 0 ? "axis" : "grid" }, svgB);
      T(svgB, L - 6, ty(v) + 6, String(v), { "text-anchor": "end", class: "sm" });
    }
    const tstep = st.n <= 10 ? 1 : st.n <= 20 ? 2 : st.n <= 50 ? 5 : st.n <= 100 ? 10 : 20;
    for (let t = 0; t <= st.n; t += tstep) {
      svgEl("line", { x1: tx(t), x2: tx(t), y1: bot, y2: bot + 5, class: "axis" }, svgB);
      T(svgB, tx(t), bot + 20, String(t), { "text-anchor": "middle", class: "sm" });
    }
    T(svgB, (L + R) / 2, H - 6, "year t", { "text-anchor": "middle" });
    yLabel(svgB, 16, "Nₜ");
    for (const s of ss) {
      if (s.x === 0) continue;
      const c = Math.abs(s.d) < 1 ? COL.stable : COL.unstable;
      svgEl("line", { x1: L, x2: R, y1: ty(s.x), y2: ty(s.x), stroke: c, "stroke-width": 2, "stroke-dasharray": "7 5" }, svgB);
    }
    const pgB = clipRect(svgB, `cw-clipB-${uid}`, L - 8, R + 8);
    svgEl("polyline", { points: xs.map((v, t) => `${tx(t).toFixed(1)},${ty(v).toFixed(1)}`).join(" "), fill: "none", stroke: COL.path, "stroke-width": 2 }, pgB);
    const rad = st.n <= 60 ? 4.5 : st.n <= 120 ? 3 : 2;
    xs.forEach((v, t) => svgEl("circle", { cx: tx(t), cy: ty(v), r: rad, fill: COL.path }, pgB));
    // legend, on the title line
    const lg = (x: number, c: string, s: string, filled: boolean) => {
      svgEl("circle", { cx: x, cy: 16, r: 6, fill: filled ? c : "#fff", stroke: c, "stroke-width": filled ? 0 : 3 }, svgB);
      T(svgB, x + 11, 22, s, { class: "sm" });
    };
    lg(R - 228, COL.stable, "attracting", true);
    lg(R - 118, COL.unstable, "repelling", false);

    // ---- readout ----
    ro.replaceChildren();
    const item = (k: string, v: string, cls = "") => { const s = htmlEl("span", cls, "", ro); if (k) htmlEl("b", "", k + " ", s); s.append(v); return s; };
    const last = xs[xs.length - 1];
    item(`N₀ = ${st.N0.toFixed(1)},`, `N${subs(st.n)} = ${last.toFixed(3)}`);
    const n = xs.length - 1;
    if (o.collapse < 0 && n >= 2) {
      const d1 = xs[n] - xs[n - 1], d0 = xs[n - 1] - xs[n - 2];
      item(`step ratio d${subs(n - 1)}/d${subs(n - 2)} =`, Math.abs(d1) > 1e-12 && d0 !== 0 ? fx(d1 / d0, 3) : "– (steps below 10⁻¹²)");
    }
    const lr = longRun(p, st.N0);
    const lrTxt = lr.kind === "collapse" ? `collapses (N ≤ 0) in year ${lr.year}`
      : lr.kind === "none" ? "no repeat within 64 years: chaos"
      : lr.kind === "slow" ? "still creeping towards a steady state (|g′| ≈ 1)"
      : lr.values.length === 1 ? `settles at ${lr.values[0].toFixed(3)}`
      : `${lr.values.length}-cycle ${lr.values.map(v => v.toFixed(2)).join(", ")}`;
    item("long run (after 3000 years):", lrTxt, lr.kind === "collapse" ? "cw-warn" : "");
    const brk = htmlEl("span", "", "", ro); brk.style.flexBasis = "100%"; brk.style.height = "0";
    if (ss.length === 0) {
      const [hm, nm] = maxCatch(st.r);
      item("no steady state:", `H = ${st.H.toFixed(2)} is above the maximum sustainable catch ${hm.toFixed(2)} (at N = ${nm.toFixed(1)})`, "cw-warn");
    } else {
      for (const s of ss) {
        const sp = item("", "");
        const dot = htmlEl("span", "", Math.abs(s.d) < 1 ? "●" : "○", sp);
        dot.style.color = Math.abs(s.d) < 1 ? COL.stable : COL.unstable; dot.style.fontWeight = "700";
        sp.append(" ");
        htmlEl("b", "", `N* = ${s.x.toFixed(2)}:`, sp);
        sp.append(` g′(N*) = ${fx(s.d, 3)}, ${classify(s.d)}`);
      }
      if (st.model === "harvest") { const [hm] = maxCatch(st.r); item("max. catch", `${hm.toFixed(2)}`); }
    }
  }

  function drawBisect(): void {
    setSizes(640, 460);
    const p = P(), hi = range(p);
    const L = 70, R = WA - 14;
    ax = scale(0, hi, L, R);
    st.a = Math.min(st.a, hi); st.b = Math.min(st.b, hi);
    // y range: show the hump of f clearly
    let fmax = -Infinity, fmin = Infinity;
    for (let i = 0; i <= 600; i++) { const v = f(p, hi * i / 600); fmax = Math.max(fmax, v); fmin = Math.min(fmin, v); }
    const pos = Math.max(0, fmax);
    const y1 = Math.max(3, 1.2 * pos), y0 = Math.max(fmin, -Math.max(16, 2.5 * pos));
    fy = scale(y0, y1, bot, top);
    svgA.setAttribute("aria-label", "f(N) = g(N) − N with the bisection bracket");
    svgB.setAttribute("aria-label", "Bisection error against step, log scale");

    const steps = bisect(p, st.a, st.b, st.k + 1);
    if (steps && steps.length && st.k > steps.length - 1) st.k = steps.length - 1;   // f(m) = 0 exactly: done
    const fa = f(p, st.a), fb = f(p, st.b);
    const ok = steps !== null;
    let ca = st.a, cb = st.b;
    let lastM: number | null = null;
    if (ok && steps.length) {
      const k = Math.min(st.k, steps.length - 1);
      const cur = steps[k];
      ca = cur.a; cb = cur.b;
      if (k > 0) lastM = steps[k - 1].m;
    }
    // reference root: full-precision bisection of the starting bracket
    const root = ok ? (steps.length ? bisect(p, st.a, st.b, 200)!.slice(-1)[0].m : (fa === 0 ? st.a : st.b)) : NaN;

    // ---- (a) f with bracket ----
    svgA.replaceChildren();
    T(svgA, 8, 22, `(a) f(N) = g(N) − N, bracket after ${st.k} step${st.k === 1 ? "" : "s"}`, { class: "ttl" });
    const yStep = niceStep((y1 - y0) / 5);
    for (let v = Math.ceil(y0 / yStep) * yStep; v <= y1 + 1e-9; v += yStep) {
      svgEl("line", { x1: L, x2: R, y1: fy(v), y2: fy(v), class: "grid" }, svgA);
      T(svgA, L - 6, fy(v) + 6, fx(v, yStep < 1 ? 1 : 0), { "text-anchor": "end", class: "sm" });
    }
    for (const v of ticksFor(hi)) {
      svgEl("line", { x1: ax(v), x2: ax(v), y1: bot, y2: bot + 5, class: "axis" }, svgA);
      T(svgA, ax(v), bot + 20, String(v), { "text-anchor": "middle", class: "sm" });
    }
    svgEl("line", { x1: L, x2: R, y1: bot, y2: bot, class: "axis" }, svgA);
    T(svgA, (L + R) / 2, H - 6, "N", { "text-anchor": "middle" });
    yLabel(svgA, 16, "f(N)");
    const pg = clipRect(svgA, `cw-clipC-${uid}`, L, R);
    // current bracket shading
    svgEl("rect", { x: ax(ca), y: top, width: Math.max(2, ax(cb) - ax(ca)), height: bot - top, fill: COL.sky, "fill-opacity": ok ? 0.25 : 0.12 }, pg);
    svgEl("line", { x1: L, x2: R, y1: fy(0), y2: fy(0), stroke: COL.ink, "stroke-width": 1.5 }, pg);
    const pts: string[] = [];
    for (let i = 0; i <= 600; i++) { const N = hi * i / 600; pts.push(`${ax(N).toFixed(1)},${fy(Math.max(y0 - 5, Math.min(y1 + 5, f(p, N)))).toFixed(1)}`); }
    svgEl("polyline", { points: pts.join(" "), fill: "none", stroke: COL.map, "stroke-width": 3 }, pg);
    if (lastM !== null) {
      svgEl("line", { x1: ax(lastM), x2: ax(lastM), y1: fy(0), y2: fy(f(p, lastM)), stroke: COL.root, "stroke-width": 2, "stroke-dasharray": "4 3" }, pg);
      svgEl("circle", { cx: ax(lastM), cy: fy(f(p, lastM)), r: 6, fill: COL.root }, pg);
    }
    // starting bracket handles
    for (const [x, fv, lab] of [[st.a, fa, "a"], [st.b, fb, "b"]] as [number, number, string][]) {
      svgEl("line", { x1: ax(x), x2: ax(x), y1: fy(0), y2: fy(Math.max(y0, Math.min(y1, fv))), stroke: COL.brk, "stroke-width": 2 }, pg);
      svgEl("rect", { x: ax(x) - 8, y: fy(Math.max(y0, Math.min(y1, fv))) - 8, width: 16, height: 16, fill: COL.brk }, pg);
      svgEl("circle", { cx: ax(x), cy: fy(0), r: 10, fill: "#fff", stroke: COL.brk, "stroke-width": 3, class: "handle" }, svgA);
      T(svgA, ax(x), fy(0) + (fv < 0 ? -16 : 30), lab, { "text-anchor": "middle", "font-weight": 700 }, COL.brk);
    }
    if (!ok) T(svgA, R - 4, top + 26, "f(a), f(b) same sign: no bracket", { "text-anchor": "end", "font-weight": 700 }, "#b00020");

    // ---- (b) error against step ----
    svgB.replaceChildren();
    T(svgB, 8, 22, "(b) error of the midpoint", { class: "ttl" });
    const BL = 66, BR = WB - 12, KMAX = 60, E0 = -16, E1 = 2;
    const bx = scale(0, KMAX, BL, BR), by = scale(E0, E1, bot, top);
    for (let e = E0; e <= E1; e += 2) {
      svgEl("line", { x1: BL, x2: BR, y1: by(e), y2: by(e), class: e === E0 ? "axis" : "grid" }, svgB);
      T(svgB, BL - 6, by(e) + 6, pow10(e), { "text-anchor": "end", class: "sm" });
    }
    for (let k = 0; k <= KMAX; k += 10) {
      svgEl("line", { x1: bx(k), x2: bx(k), y1: bot, y2: bot + 5, class: "axis" }, svgB);
      T(svgB, bx(k), bot + 20, String(k), { "text-anchor": "middle", class: "sm" });
    }
    T(svgB, (BL + BR) / 2, H - 6, "step k", { "text-anchor": "middle" });
    yLabel(svgB, 14, "error");
    const pgB = clipRect(svgB, `cw-clipD-${uid}`, BL - 6, BR + 6);
    const ly = (v: number) => by(Math.max(E0, Math.log10(Math.max(v, 1e-30))));
    if (ok && steps.length) {
      const all = bisect(p, st.a, st.b, KMAX + 1)!;
      const bound = all.map(s => (s.b - s.a) / 2);
      svgEl("polyline", { points: bound.map((v, k) => `${bx(k)},${ly(v)}`).join(" "), fill: "none", stroke: COL.brk, "stroke-width": 2, "stroke-opacity": 0.35 }, pgB);
      const kk = Math.min(st.k, all.length - 1);
      svgEl("polyline", { points: bound.slice(0, kk + 1).map((v, k) => `${bx(k)},${ly(v)}`).join(" "), fill: "none", stroke: COL.brk, "stroke-width": 3 }, pgB);
      for (let k = 0; k <= kk; k++) {
        const e = Math.abs(all[k].m - root);
        svgEl("circle", { cx: bx(k), cy: ly(e === 0 ? 1e-30 : e), r: 4.5, fill: "#fff", stroke: COL.ink, "stroke-width": 2 }, pgB);
      }
      svgEl("line", { x1: bx(st.k), x2: bx(st.k), y1: top, y2: bot, stroke: COL.ink, "stroke-dasharray": "2 3" }, svgB);
      // legend
      svgEl("line", { x1: BR - 200, x2: BR - 176, y1: top + 8, y2: top + 8, stroke: COL.brk, "stroke-width": 3 }, svgB);
      T(svgB, BR - 170, top + 14, "bound (b − a)/2", {}, "#333");
      svgEl("circle", { cx: BR - 188, cy: top + 34, r: 5, fill: "#fff", stroke: COL.ink, "stroke-width": 2 }, svgB);
      T(svgB, BR - 170, top + 40, "actual |m − N*|", {}, "#333");
    } else {
      T(svgB, (BL + BR) / 2, (top + bot) / 2, "bisection cannot start", { "text-anchor": "middle", "font-weight": 700 }, "#b00020");
    }

    // ---- readout ----
    ro.replaceChildren();
    const item = (k: string, v: string, cls = "") => { const s = htmlEl("span", cls, "", ro); if (k) htmlEl("b", "", k + " ", s); s.append(v); return s; };
    item("start", `[${st.a}, ${st.b}]: f(a) = ${fx(fa, 3)}, f(b) = ${fx(fb, 3)}`);
    if (!ok) {
      item("", "no sign change: bisection cannot start", "cw-warn");
      const ss = steadyStates(p, hi);
      if (st.model === "harvest") {
        const [hm, nm] = maxCatch(st.r);
        item("max. catch", `${hm.toFixed(4)} at N = ${nm.toFixed(2)}; roots on [0, ${hi}]: ${ss.length ? ss.map(s => s.x.toFixed(2)).join(", ") : "none"}`);
      }
    } else if (steps.length) {
      const w = cb - ca, m = 0.5 * (ca + cb);
      item(`after ${st.k} step${st.k === 1 ? "" : "s"}:`, `[${ca.toPrecision(10)}, ${cb.toPrecision(10)}]`);
      item("width", sci(w, 2));
      item("midpoint", `${m.toPrecision(12)}, f(m) = ${sci(f(p, m), 2)}`);
      item("error", `${sci(Math.abs(m - root), 2)} (bound ${sci(w / 2, 2)})`);
      const need = Math.ceil(Math.log2((st.b - st.a) / 1e-8));
      item("steps to width 10⁻⁸:", String(need));
    } else item("", "f is exactly 0 at a bracket end");
  }

  rebuild();
}

function niceStep(raw: number): number {
  const e = Math.pow(10, Math.floor(Math.log10(raw)));
  const m = raw / e;
  return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * e;
}

autoMount(NAME, mount);

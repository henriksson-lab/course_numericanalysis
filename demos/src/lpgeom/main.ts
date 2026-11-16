// LP geometry explorer for the constrained-lp unit.
// Tab 1, "LP in 2 variables": the lecture's ATP LP, max c1 x1 + c2 x2 (respiration,
// fermentation) subject to carbon x1 + x2 <= 6, proteome 4 x1 + x2 <= 12, and optionally
// oxygen x1 <= o and a demand x1 + x2 >= d. The feasible polygon, the objective as parallel
// level lines with a draggable gradient arrow, the optimal corner (or optimal edge), a Step
// button that walks the simplex path corner by corner, shadow prices measured by re-solving
// (left and right of a kink), and a scan of the best ATP against one limit.
// Tab 2, "FBA toy network": the 9-reaction network of the lecture, solved by a small dense
// simplex in TS; knock-outs, uptake limits, FVA ranges, shadow prices, oxygen scan.
// Auto-mounts on <div data-demo="lpgeom">. Optional attributes:
//   data-preset="lecture|oxygen|degenerate|tie|unbounded|infeasible"
//   data-fba="1|0" (show the FBA tab)   data-tab="lp|fba"
import { svgEl, htmlEl, slider, button, presets, checkbox, svgPoint, draggable, autoMount, Slider } from "../lib/ui";
import { Vec, HalfPlane, solve2, edges, isFar, BIG, active, fba, fva, maxResidual, RXN, FBAOpts, Edge } from "./lp";

const NAME = "lpgeom";
// Okabe-Ito, as in make-constrained-lp-figs.R
const OI = { orange: "#E69F00", sky: "#56B4E9", green: "#009E73", blue: "#0072B2", red: "#D55E00", purple: "#CC79A7", grey: "#999999", ink: "#222" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .45em; }
[data-demo="${NAME}"] .lg-main { display: flex; flex-wrap: wrap; gap: 12px; width: 100%; align-items: flex-start; }
[data-demo="${NAME}"] .lg-plot { flex: 1 1 380px; max-width: 470px; min-width: 0; }
[data-demo="${NAME}"] .lg-plot svg { width: 100%; height: auto; touch-action: none; }
[data-demo="${NAME}"] .lg-col { flex: 1 1 300px; display: flex; flex-direction: column; gap: .4em; min-width: 0; }
[data-demo="${NAME}"] .lg-row { display: flex; flex-wrap: wrap; gap: .35em .8em; align-items: center; }
[data-demo="${NAME}"] .lg-con { display: grid; grid-template-columns: 1.4em 1fr; align-items: center; gap: .2em; }
[data-demo="${NAME}"] .lg-con .qb-slider { grid-template-columns: 7.2em 1fr 2.8em; }
[data-demo="${NAME}"] .lg-con input[type=checkbox] { width: 1.1em; height: 1.1em; margin: 0; }
[data-demo="${NAME}"] .lg-obj .qb-slider { grid-template-columns: 7.2em 1fr 2.8em; }
[data-demo="${NAME}"] .lg-h { font-weight: 600; margin-top: .15em; }
[data-demo="${NAME}"] .qb-readout { white-space: normal; line-height: 1.4; }
[data-demo="${NAME}"] .qb-readout div { margin: .08em 0; }
[data-demo="${NAME}"] .lg-bad { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .lg-ok { color: #8a3000; font-weight: 600; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 19px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #444; }
[data-demo="${NAME}"] .lg-ko { display: grid; grid-template-columns: repeat(3, auto); gap: .1em .8em; font-size: .9em; }
[data-demo="${NAME}"] .lg-fplot { flex: 1 1 420px; max-width: 600px; min-width: 0; display: flex; flex-direction: column; gap: 4px; }
[data-demo="${NAME}"] .lg-fplot svg { width: 100%; height: auto; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; }
.reveal [data-demo="${NAME}"] .lg-main { flex-wrap: nowrap; }
.reveal [data-demo="${NAME}"] .lg-plot { flex: 0 0 420px; }
.reveal [data-demo="${NAME}"] .lg-col { flex: 0 0 345px; }
.reveal [data-demo="${NAME}"] .lg-col.lg-right { flex: 0 0 350px; }
.reveal [data-demo="${NAME}"] .lg-fplot { flex: 0 0 600px; }
.reveal [data-demo="${NAME}"] .lg-col.lg-fside { flex: 0 0 520px; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; }
.reveal [data-demo="${NAME}"] .qb-btn { font-size: .78em; padding: .2em .5em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- number formatting: 2 -> "2", 1.5 -> "1.5", 0.83333 -> "0.833" ----
const fmt = (v: number, d = 3): string => {
  if (!Number.isFinite(v)) return v > 0 ? "∞" : "–";
  const s = (Math.abs(v) < 5e-10 ? 0 : v).toFixed(d);
  return s.includes(".") ? s.replace(/0+$/, "").replace(/\.$/, "") : s;
};
const pt = (v: Vec): string => `(${fmt(v[0], 2)}, ${fmt(v[1], 2)})`;

// ---- the limits of the ATP LP ----
interface Limit { key: string; name: string; label: string; a: Vec; ge: boolean; min: number; max: number; step: number; def: number; col: string; dash: string }
const LIMITS: Limit[] = [
  { key: "carbon", name: "carbon", label: "carbon x₁+x₂ ≤", a: [1, 1], ge: false, min: 0, max: 10, step: 0.5, def: 6, col: OI.blue, dash: "" },
  { key: "proteome", name: "proteome", label: "proteome 4x₁+x₂ ≤", a: [4, 1], ge: false, min: 0, max: 24, step: 0.5, def: 12, col: OI.green, dash: "" },
  { key: "oxygen", name: "oxygen", label: "oxygen x₁ ≤", a: [1, 0], ge: false, min: 0, max: 8, step: 0.5, def: 1.5, col: OI.purple, dash: "" },
  { key: "demand", name: "demand", label: "demand x₁+x₂ ≥", a: [1, 1], ge: true, min: 0, max: 12, step: 0.5, def: 8, col: "#555555", dash: "9 6" },
];
const hp = (L: Limit, b: number): HalfPlane => L.ge ? { a: [-L.a[0], -L.a[1]], b: -b } : { a: L.a, b };
const NONNEG: HalfPlane[] = [{ a: [-1, 0], b: 0 }, { a: [0, -1], b: 0 }];

interface Preset { key: string; label: string; on: boolean[]; b: number[]; c: Vec; focus: number }
const PRESETS: Preset[] = [
  { key: "lecture", label: "ATP LP", on: [true, true, false, false], b: [6, 12, 1.5, 8], c: [6, 3], focus: 0 },
  { key: "oxygen", label: "oxygen ≤ 1.5", on: [true, true, true, false], b: [6, 12, 1.5, 8], c: [6, 3], focus: 2 },
  { key: "degenerate", label: "degenerate corner", on: [true, true, true, false], b: [6, 12, 2, 8], c: [6, 3], focus: 2 },
  { key: "tie", label: "tie: 8x₁+2x₂", on: [true, true, false, false], b: [6, 12, 1.5, 8], c: [8, 2], focus: 1 },
  { key: "unbounded", label: "unbounded", on: [false, false, true, false], b: [6, 12, 3, 8], c: [6, 3], focus: 2 },
  { key: "infeasible", label: "infeasible", on: [true, false, false, true], b: [6, 12, 1.5, 8], c: [6, 3], focus: 3 },
];

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const uid = ++mountCount;
  const withFba = ds.fba !== "0";
  const lpBox = htmlEl("div", "", ""), fbaBox = htmlEl("div", "", "");
  for (const b of [lpBox, fbaBox]) Object.assign(b.style, { display: "flex", flexDirection: "column", gap: ".45em", width: "100%" });
  if (withFba) {
    const tabs = htmlEl("div", "lg-row", "", host);
    const tb = presets(tabs, [
      { label: "LP in 2 variables", apply: () => show("lp") },
      { label: "FBA toy network", apply: () => show("fba") },
    ]);
    const show = (t: string) => {
      lpBox.style.display = t === "lp" ? "flex" : "none";
      fbaBox.style.display = t === "fba" ? "flex" : "none";
      tb.forEach((b, i) => b.classList.toggle("active", (i === 0) === (t === "lp")));
    };
    host.append(lpBox, fbaBox);
    show(ds.tab === "fba" ? "fba" : "lp");
    buildFba(fbaBox, uid);
  } else host.append(lpBox);
  buildLp(lpBox, uid, ds.preset ?? "lecture");
}

// ================================================================ tab 1: 2-variable LP
function buildLp(box: HTMLElement, uid: number, presetKey: string): void {
  const p0 = PRESETS.find(p => p.key === presetKey) ?? PRESETS[0];
  const st = { on: p0.on.slice(), b: p0.b.slice(), c: [p0.c[0], p0.c[1]] as Vec, focus: p0.focus };
  let path: number[] = [];          // vertex indices visited by the simplex walk
  let pathMsg: string[] = [];
  let walkDone = false;

  const main = htmlEl("div", "lg-main", "", box);
  const plotDiv = htmlEl("div", "lg-plot", "", main);
  const W = 430, H = 500;
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Feasible polygon and objective lines" }, plotDiv);
  const ctl = htmlEl("div", "lg-col", "", main);
  const right = htmlEl("div", "lg-col lg-right", "", main);

  // ---- controls ----
  const pBtns = presets(ctl, PRESETS.map(p => ({ label: p.label, apply: () => applyPreset(p) })));
  htmlEl("div", "lg-h", "limits (tick to switch on)", ctl);
  const chks: HTMLInputElement[] = [], sls: Slider[] = [];
  LIMITS.forEach((L, k) => {
    const row = htmlEl("div", "lg-con", "", ctl);
    const cb = htmlEl("input", "", "", row);
    cb.type = "checkbox"; cb.checked = st.on[k]; cb.style.accentColor = L.col;
    cb.addEventListener("change", () => { st.on[k] = cb.checked; if (cb.checked) st.focus = k; cb.blur(); changed(); });
    const s = slider(row, L.label, L.min, L.max, L.step, st.b[k], v => { st.b[k] = v; st.focus = k; changed(); }, v => fmt(v, 1));
    const nm = s.input.parentElement?.querySelector<HTMLElement>(".qb-slider-name");
    if (nm) { nm.style.color = L.col; nm.style.fontWeight = "600"; nm.style.whiteSpace = "normal"; nm.style.lineHeight = "1.1"; }
    s.input.style.accentColor = L.col;
    chks.push(cb); sls.push(s);
  });
  htmlEl("div", "lg-h", "objective: max ATP = c₁x₁ + c₂x₂ (or drag the arrow)", ctl);
  const obj = htmlEl("div", "lg-obj", "", ctl);
  const c1S = slider(obj, "c₁ (resp.)", 0, 10, 0.5, st.c[0], v => { st.c[0] = v; changed(); }, v => fmt(v, 1));
  const c2S = slider(obj, "c₂ (ferm.)", 0, 10, 0.5, st.c[1], v => { st.c[1] = v; changed(); }, v => fmt(v, 1));
  for (const s of [c1S, c2S]) s.input.style.accentColor = OI.orange;
  const sr = htmlEl("div", "lg-row", "", ctl);
  htmlEl("span", "lg-h", "simplex:", sr);
  const stepB = button(sr, "step ▶", () => step());
  button(sr, "all ⏭", () => { for (let i = 0; i < 20 && !walkDone; i++) step(); });
  button(sr, "restart path", () => { resetPath(); draw(); });

  // ---- right column: readout and scan ----
  const ro = htmlEl("div", "qb-readout", "", right);
  const scanRow = htmlEl("div", "lg-row", "", right);
  htmlEl("span", "", "scan:", scanRow);
  const fBtns = presets(scanRow, LIMITS.map((L, k) => ({ label: L.name, apply: () => { st.focus = k; draw(); } })));
  const SW = 350, SH = 230;
  const scan = svgEl("svg", { viewBox: `0 0 ${SW} ${SH}`, role: "img", "aria-label": "Best ATP against one limit" }, right);
  scan.style.width = "100%"; scan.style.height = "auto";

  // ---- geometry of the main plot (equal scales, so the gradient is perpendicular) ----
  const X0 = 56, Y0 = 446, PX = 32, XM = 11, YM = 13;
  const sx = (v: number) => X0 + PX * v, sy = (v: number) => Y0 - PX * v;
  const CA: Vec = [258, 170], CK = 10;     // compass anchor (px) and px per unit of c

  const halfPlanes = (b = st.b, on = st.on, skip = -1): HalfPlane[] =>
    LIMITS.map((L, k) => (on[k] && k !== skip) ? hp(L, b[k]) : null).filter((h): h is HalfPlane => h !== null);
  const allCons = (): { h: HalfPlane; name: string }[] => [
    { h: NONNEG[0], name: "x₁ ≥ 0" }, { h: NONNEG[1], name: "x₂ ≥ 0" },
    ...LIMITS.map((L, k) => st.on[k] ? { h: hp(L, st.b[k]), name: L.name } : null).filter((x): x is { h: HalfPlane; name: string } => x !== null),
  ];

  function applyPreset(p: Preset): void {
    st.on = p.on.slice(); st.b = p.b.slice(); st.c = [p.c[0], p.c[1]]; st.focus = p.focus;
    chks.forEach((cb, k) => cb.checked = st.on[k]);
    sls.forEach((s, k) => s.set(st.b[k]));
    c1S.set(st.c[0]); c2S.set(st.c[1]);
    changed(true);
  }
  function changed(keepPreset = false): void {
    if (!keepPreset) pBtns.forEach(b => b.classList.remove("active"));
    resetPath(); draw();
  }
  function resetPath(): void { path = []; pathMsg = []; walkDone = false; }

  // a closed-form "phase 1": the origin if feasible, otherwise the worst real corner
  function startVertex(poly: Vec[]): number {
    const o = poly.findIndex(v => Math.abs(v[0]) + Math.abs(v[1]) < 1e-9);
    if (o >= 0) return o;
    let best = -1, bz = Infinity;
    poly.forEach((v, i) => { if (!isFar(v)) { const z = st.c[0] * v[0] + st.c[1] * v[1]; if (z < bz) { bz = z; best = i; } } });
    return best;
  }

  const atp = (v: Vec) => st.c[0] * v[0] + st.c[1] * v[1];
  const perUnit = (e: Edge, names: string[]): string => {
    const f = e.freed.map(k => names[k]);
    if (!f.length) return "per unit length";
    const n = f[0];
    return n === "x₁ ≥ 0" ? "per unit x₁" : n === "x₂ ≥ 0" ? "per unit x₂" : `per unit of ${n} freed`;
  };

  function step(): void {
    const sol = solve2(halfPlanes(), st.c);
    if (sol.status === "infeasible" || walkDone) return;
    const poly = sol.poly, cons = allCons(), names = cons.map(c => c.name), hs = cons.map(c => c.h);
    if (!path.length) {
      const s0 = startVertex(poly);
      if (s0 < 0) return;
      path.push(s0);
      pathMsg.push(`start at ${pt(poly[s0])}, ATP ${fmt(atp(poly[s0]))}` +
        (Math.abs(poly[s0][0]) + Math.abs(poly[s0][1]) > 1e-9 ? " (the origin is not feasible: phase 1 finds a corner)" : ""));
      draw(); return;
    }
    const i = path[path.length - 1], v = poly[i];
    const es = edges(poly, i, hs, st.c);
    const desc = es.map(e => `${e.rate >= 0 ? "+" : "−"}${fmt(Math.abs(e.rate))} ${perUnit(e, names)}`).join(", ");
    const best = es.reduce<Edge | null>((a, e) => e.rate > 1e-9 && (!a || e.rate > a.rate) ? e : a, null);
    if (!best) {
      walkDone = true;
      const zero = es.find(e => Math.abs(e.rate) <= 1e-9);
      pathMsg.push(`edges from ${pt(v)}: ${desc || "none"}. No improving edge: optimal.` +
        (zero ? ` Rate 0 along the ${zero.along.map(k => names[k]).join("/")} edge: a tie.` : ""));
    } else if (best.far) {
      walkDone = true; path.push(best.to);
      pathMsg.push(`edges from ${pt(v)}: ${desc}. The improving edge never ends: unbounded.`);
    } else {
      path.push(best.to);
      const w = poly[best.to];
      pathMsg.push(`${pt(v)} → ${pt(w)} along ${best.along.map(k => names[k]).join("/")}: ${desc.split(", ")[es.indexOf(best)]}; ATP ${fmt(atp(v))} → ${fmt(atp(w))}`);
    }
    draw();
  }

  // shadow price of limit k: change of the best ATP per unit of its right-hand side
  function price(k: number): [number, number] {
    const h = 1e-4;
    const z = (bk: number) => {
      const b = st.b.slice(); b[k] = bk;
      const s = solve2(halfPlanes(b), st.c);
      return s.status === "optimal" || s.status === "tie" ? s.z : NaN;
    };
    const z0 = z(st.b[k]);
    return [(z0 - z(st.b[k] - h)) / h, (z(st.b[k] + h) - z0) / h];
  }

  // drag the gradient arrow
  const toC = (ev: PointerEvent): [number, number] => {
    const [x, y] = svgPoint(svg, ev);
    const q = (v: number) => Math.max(0, Math.min(10, Math.round(v * 2) / 2));
    let a = q((x - CA[0]) / CK), b = q((CA[1] - y) / CK);
    if (a === 0 && b === 0) a = 0.5;
    return [a, b];
  };

  // segment of the line a.x = b inside the view, ordered by decreasing x2
  function lineSeg(a: Vec, b: number): [Vec, Vec] | null {
    const ps: Vec[] = [];
    const add = (p: Vec) => { if (p[0] >= -1e-9 && p[0] <= XM + 1e-9 && p[1] >= -1e-9 && p[1] <= YM + 1e-9) ps.push(p); };
    if (Math.abs(a[1]) > 1e-12) { add([0, b / a[1]]); add([XM, (b - a[0] * XM) / a[1]]); }
    if (Math.abs(a[0]) > 1e-12) { add([b / a[0], 0]); add([(b - a[1] * YM) / a[0], YM]); }
    if (ps.length < 2) return null;
    ps.sort((p, q) => q[1] - p[1] || p[0] - q[0]);
    return [ps[0], ps[ps.length - 1]];
  }
  const clipView = (v: Vec): Vec => [Math.min(v[0], XM + 1), Math.min(v[1], YM + 1)];

  function draw(): void {
    const hs = halfPlanes();
    const sol = solve2(hs, st.c);
    const poly = sol.poly;
    fBtns.forEach((b, k) => { b.classList.toggle("active", k === st.focus); b.disabled = !st.on[k]; b.style.opacity = st.on[k] ? "1" : ".4"; });
    stepB.disabled = sol.status === "infeasible" || walkDone;

    svg.replaceChildren();
    const defs = svgEl("defs", {}, svg);
    const mk = svgEl("marker", { id: `lg-arr-${uid}`, viewBox: "0 0 10 10", refX: 8, refY: 5, markerWidth: 4, markerHeight: 4, orient: "auto-start-reverse" }, defs);
    svgEl("path", { d: "M0,0 L10,5 L0,10 z", fill: OI.red }, mk);
    const mk2 = svgEl("marker", { id: `lg-arro-${uid}`, viewBox: "0 0 10 10", refX: 8, refY: 5, markerWidth: 4, markerHeight: 4, orient: "auto" }, defs);
    svgEl("path", { d: "M0,0 L10,5 L0,10 z", fill: OI.orange }, mk2);
    const cp = svgEl("clipPath", { id: `lg-clip-${uid}` }, defs);
    svgEl("rect", { x: sx(0), y: sy(YM), width: PX * XM, height: PX * YM }, cp);

    // axes and grid
    for (let v = 0; v <= XM; v += 2) {
      svgEl("line", { x1: sx(v), x2: sx(v), y1: sy(0), y2: sy(YM), class: "grid" }, svg);
      svgEl("text", { x: sx(v), y: sy(0) + 22, "text-anchor": "middle", class: "sm" }, svg).textContent = String(v);
    }
    for (let v = 0; v <= YM; v += 2) {
      svgEl("line", { x1: sx(0), x2: sx(XM), y1: sy(v), y2: sy(v), class: "grid" }, svg);
      svgEl("text", { x: sx(0) - 7, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = String(v);
    }
    svgEl("line", { x1: sx(0), x2: sx(XM), y1: sy(0), y2: sy(0), class: "axis", "stroke-width": 1.5 }, svg);
    svgEl("line", { x1: sx(0), x2: sx(0), y1: sy(0), y2: sy(YM), class: "axis", "stroke-width": 1.5 }, svg);
    svgEl("text", { x: sx(XM / 2), y: H - 10, "text-anchor": "middle" }, svg).textContent = "respiration flux x₁";
    svgEl("text", { x: 16, y: sy(YM / 2), "text-anchor": "middle", transform: `rotate(-90 16 ${sy(YM / 2)})` }, svg).textContent = "fermentation flux x₂";

    const g = svgEl("g", { "clip-path": `url(#lg-clip-${uid})` }, svg);
    const polyPts = (p: Vec[]) => p.map(v => `${sx(Math.min(v[0], 50))},${sy(Math.min(v[1], 50))}`).join(" ");
    // shaded regions
    if (sol.status === "infeasible") {
      const dk = LIMITS.findIndex(L => L.ge);
      if (st.on[dk]) {
        const rest = solve2(halfPlanes(st.b, st.on, dk), st.c).poly;
        if (rest.length > 2) svgEl("polygon", { points: polyPts(rest), fill: OI.sky, "fill-opacity": 0.3 }, g);
        const dem = solve2([hp(LIMITS[dk], st.b[dk])], st.c).poly;
        if (dem.length > 2) svgEl("polygon", { points: polyPts(dem), fill: "#777", "fill-opacity": 0.15 }, g);
      }
    } else if (poly.length > 2) svgEl("polygon", { points: polyPts(poly), fill: OI.sky, "fill-opacity": 0.38 }, g);

    // objective level lines
    const cn = st.c;
    if (cn[0] > 0 || cn[1] > 0) {
      const zv = cn[0] * XM + cn[1] * YM;
      const levels: number[] = (sol.status === "optimal" || sol.status === "tie") && sol.z > 0
        ? [0.25, 0.5, 0.75, 1.25, 1.5, 2].map(f => f * sol.z) : [1, 2, 3, 4, 5, 6, 7].map(k => k * zv / 8);
      for (const z of levels) {
        const s = lineSeg(cn, z);
        if (s) svgEl("line", { x1: sx(s[0][0]), y1: sy(s[0][1]), x2: sx(s[1][0]), y2: sy(s[1][1]), stroke: OI.orange, "stroke-width": 2, "stroke-dasharray": "7 6", "stroke-opacity": 0.85 }, g);
      }
    }
    // constraint lines
    LIMITS.forEach((L, k) => {
      if (!st.on[k]) return;
      const s = lineSeg(L.a, st.b[k]);
      if (!s) return;
      svgEl("line", { x1: sx(s[0][0]), y1: sy(s[0][1]), x2: sx(s[1][0]), y2: sy(s[1][1]), stroke: L.col, "stroke-width": 3, "stroke-dasharray": L.dash }, g);
    });
    // optimal level line, optimal edge, optimum
    if ((sol.status === "optimal" || sol.status === "tie") && (cn[0] > 0 || cn[1] > 0)) {
      const s = lineSeg(cn, sol.z);
      if (s) svgEl("line", { x1: sx(s[0][0]), y1: sy(s[0][1]), x2: sx(s[1][0]), y2: sy(s[1][1]), stroke: OI.orange, "stroke-width": 3.5 }, g);
      if (sol.status === "tie" && sol.opt.length >= 2) {
        const [a, b] = [clipView(poly[sol.opt[0]]), clipView(poly[sol.opt[1]])];
        svgEl("line", { x1: sx(a[0]), y1: sy(a[1]), x2: sx(b[0]), y2: sy(b[1]), stroke: OI.red, "stroke-width": 8, "stroke-linecap": "round" }, g);
      }
    }
    // simplex path
    for (let i = 1; i < path.length; i++) {
      const a = poly[path[i - 1]], b = clipView(poly[path[i]]);
      const len = Math.hypot(sx(b[0]) - sx(a[0]), sy(b[1]) - sy(a[1]));
      const sh = Math.min(12, len / 3) / Math.max(len, 1e-9);
      svgEl("line", {
        x1: sx(a[0]) + (sx(b[0]) - sx(a[0])) * sh, y1: sy(a[1]) + (sy(b[1]) - sy(a[1])) * sh,
        x2: sx(b[0]) - (sx(b[0]) - sx(a[0])) * sh, y2: sy(b[1]) - (sy(b[1]) - sy(a[1])) * sh,
        stroke: OI.red, "stroke-width": 5, "marker-end": `url(#lg-arr-${uid})`,
      }, g);
    }
    // corners
    let cx = 0, cy = 0, nr = 0;
    poly.forEach(v => { if (!isFar(v)) { cx += v[0]; cy += v[1]; nr++; } });
    if (nr) { cx /= nr; cy /= nr; }
    poly.forEach((v, i) => {
      if (isFar(v) || v[0] > XM || v[1] > YM) return;
      const isOpt = (sol.status === "optimal" || sol.status === "tie") && sol.opt.includes(i);
      const onPath = path.includes(i);
      svgEl("circle", { cx: sx(v[0]), cy: sy(v[1]), r: isOpt ? 10 : 6, fill: isOpt ? OI.red : onPath ? OI.red : OI.ink, stroke: "#fff", "stroke-width": 2 }, svg);
      let dx = v[0] - cx, dy = v[1] - cy;
      const dl = Math.hypot(dx, dy) || 1; dx /= dl; dy /= dl;
      const lx = sx(v[0]) + dx * 16, ly = sy(v[1]) - dy * 16 + 6;
      const t = svgEl("text", { x: Math.max(sx(0) + 4, Math.min(W - 4, lx)), y: Math.min(sy(0) - 6, Math.max(sy(YM) + 14, ly)), "text-anchor": dx > 0.3 ? "start" : dx < -0.3 ? "end" : "middle", "font-weight": isOpt ? 700 : 400 }, svg);
      if (dx < -0.3 && lx < sx(0) + 60) t.setAttribute("text-anchor", "start");
      t.style.fill = isOpt ? OI.red : OI.ink;
      t.style.paintOrder = "stroke"; t.style.stroke = "#fff"; t.style.strokeWidth = "4px";
      t.textContent = pt(v);
    });
    // constraint labels
    LIMITS.forEach((L, k) => {
      if (!st.on[k]) return;
      const b = st.b[k];
      const p: Vec = L.key === "carbon" ? [b * 0.72, b * 0.28] : L.key === "proteome" ? [b / 4 * 0.12, b * 0.88]
        : L.key === "oxygen" ? [b, 12.2] : [b * 0.88, b * 0.12];
      if (p[0] > XM - 0.5 || p[1] > YM) return;
      const t = svgEl("text", { x: sx(p[0]) + 7, y: sy(p[1]) - 4, "font-weight": 600 }, svg);
      t.style.fill = L.col; t.style.paintOrder = "stroke"; t.style.stroke = "#fff"; t.style.strokeWidth = "4px";
      t.textContent = L.name;
    });
    // compass: the gradient arrow
    svgEl("rect", { x: CA[0] - 26, y: 30, width: W - CA[0] + 22, height: CA[1] - 30 + 18, rx: 8, fill: "#fff", "fill-opacity": 0.88, stroke: "#ddd" }, svg);
    const ot = svgEl("text", { x: (CA[0] - 26 + W - 4) / 2, y: 52, "text-anchor": "middle", "font-weight": 600 }, svg);
    ot.style.fill = "#8a5a00";
    ot.textContent = `ATP = ${fmt(cn[0], 1)}x₁ + ${fmt(cn[1], 1)}x₂`;
    const tip: Vec = [CA[0] + CK * cn[0], CA[1] - CK * cn[1]];
    const al = Math.hypot(tip[0] - CA[0], tip[1] - CA[1]), sh = Math.max(0, al - 13) / Math.max(al, 1e-9);
    svgEl("line", { x1: CA[0], y1: CA[1], x2: CA[0] + (tip[0] - CA[0]) * sh, y2: CA[1] + (tip[1] - CA[1]) * sh, stroke: OI.orange, "stroke-width": 5, "marker-end": `url(#lg-arro-${uid})` }, svg);
    svgEl("circle", { cx: CA[0], cy: CA[1], r: 4, fill: OI.orange }, svg);
    const hd = svgEl("circle", { cx: tip[0], cy: tip[1], r: 12, class: "handle", stroke: OI.orange }, svg);
    draggable(hd, toC, (a, b) => {
      if (a === st.c[0] && b === st.c[1]) return;
      st.c = [a, b]; c1S.set(a); c2S.set(b); changed();
    });

    // ---------- readout ----------
    ro.replaceChildren();
    const line = (html: string, cls = "") => { const d = htmlEl("div", cls, "", ro); d.innerHTML = html; };
    const cons = allCons();
    if (sol.status === "infeasible") line("<b>infeasible</b>: no point satisfies every constraint", "lg-bad");
    else if (sol.status === "unbounded") line("<b>unbounded</b>: ATP grows without limit", "lg-bad");
    else if (sol.status === "all") line("objective 0: every feasible point is optimal");
    else {
      const v = poly[sol.opt[0]];
      if (sol.status === "tie") {
        const w = poly[sol.opt[1]];
        line(`<b>tie</b>: the whole edge ${pt(v)} – ${isFar(w) ? "∞" : pt(w)} is optimal, ATP ${fmt(sol.z)}`, "lg-ok");
      } else line(`<b>optimum</b> x = ${pt(v)}, <b>ATP ${fmt(sol.z)}</b>`, "lg-ok");
      const act = cons.filter(c => active(c.h, v)).map(c => c.name);
      line(`active at ${pt(v)}: ${act.join(", ")}`);
    }
    if (sol.status === "optimal" || sol.status === "tie") {
      const items: string[] = [];
      LIMITS.forEach((L, k) => {
        if (!st.on[k]) return;
        const [l, r] = price(k);
        const s = Math.abs(l - r) < 1e-6 ? fmt(r) : `${fmt(l)} below / ${fmt(r)} above (kink)`;
        items.push(`${L.name} <b>${s}</b>`);
      });
      line(`shadow prices (ATP per unit of limit): ${items.join("; ")}`);
    }
    if (pathMsg.length) {
      const d = htmlEl("div", "", "", ro);
      d.style.borderTop = "1px solid #ccc"; d.style.paddingTop = ".2em";
      pathMsg.forEach((m, i) => { const p = htmlEl("div", "", "", d); p.textContent = `${i}. ${m}`; });
    } else if (sol.status !== "infeasible") line("<i>step ▶ walks the simplex path from (0, 0)</i>");

    drawScan();
  }

  function drawScan(): void {
    const k = st.focus, L = LIMITS[k];
    scan.replaceChildren();
    svgEl("text", { x: 6, y: 20, class: "ttl" }, scan).textContent = `best ATP vs ${L.name} limit`;
    if (!st.on[k]) {
      svgEl("text", { x: SW / 2, y: SH / 2, "text-anchor": "middle", class: "sm" }, scan).textContent = `tick ${L.name} to scan it`;
      return;
    }
    const N = 241, xs: number[] = [], zs: number[] = [];
    let unb = false;
    for (let i = 0; i < N; i++) {
      const b = L.min + (L.max - L.min) * i / (N - 1);
      const bb = st.b.slice(); bb[k] = b;
      const s = solve2(halfPlanes(bb), st.c);
      xs.push(b);
      if (s.status === "unbounded") unb = true;
      zs.push(s.status === "optimal" || s.status === "tie" || s.status === "all" ? s.z : NaN);
    }
    const fin = zs.filter(Number.isFinite);
    const zmax = fin.length ? Math.max(1, ...fin) * 1.1 : 1;
    const l = 50, r = SW - 12, t = 32, bt = SH - 46;
    const px = (v: number) => l + (v - L.min) / (L.max - L.min) * (r - l);
    const py = (v: number) => bt - v / zmax * (bt - t);
    svgEl("line", { x1: l, x2: r, y1: bt, y2: bt, class: "axis" }, scan);
    svgEl("line", { x1: l, x2: l, y1: t, y2: bt, class: "axis" }, scan);
    const tk = L.max > 12 ? 4 : 2;
    for (let v = L.min; v <= L.max + 1e-9; v += tk) {
      svgEl("line", { x1: px(v), x2: px(v), y1: bt, y2: bt + 5, class: "axis" }, scan);
      svgEl("text", { x: px(v), y: bt + 21, "text-anchor": "middle", class: "sm" }, scan).textContent = String(v);
    }
    const ystep = zmax > 60 ? 20 : zmax > 24 ? 10 : 5;
    for (let v = 0; v <= zmax; v += ystep) {
      svgEl("text", { x: l - 6, y: py(v) + 6, "text-anchor": "end", class: "sm" }, scan).textContent = String(v);
      if (v > 0) svgEl("line", { x1: l, x2: r, y1: py(v), y2: py(v), class: "grid" }, scan);
    }
    svgEl("text", { x: (l + r) / 2, y: SH - 4, "text-anchor": "middle", class: "sm" }, scan).textContent = `${L.name} limit`;
    let d = "", pen = false;
    zs.forEach((z, i) => {
      if (!Number.isFinite(z)) { pen = false; return; }
      d += `${pen ? "L" : "M"}${px(xs[i]).toFixed(1)},${py(z).toFixed(1)} `; pen = true;
    });
    if (d) svgEl("path", { d, fill: "none", stroke: L.col, "stroke-width": 3 }, scan);
    if (!fin.length) svgEl("text", { x: (l + r) / 2, y: (t + bt) / 2, "text-anchor": "middle", class: "sm" }, scan).textContent = unb ? "unbounded for every value" : "infeasible for every value";
    // current value
    const cur = solve2(halfPlanes(), st.c);
    svgEl("line", { x1: px(st.b[k]), x2: px(st.b[k]), y1: t, y2: bt, stroke: OI.ink, "stroke-dasharray": "3 4" }, scan);
    if (cur.status === "optimal" || cur.status === "tie") {
      svgEl("circle", { cx: px(st.b[k]), cy: py(cur.z), r: 7, fill: OI.red, stroke: "#fff", "stroke-width": 2 }, scan);
      const [pl, pr] = price(k);
      const txt = Math.abs(pl - pr) < 1e-6 ? `slope ${fmt(pr)}` : `slope ${fmt(pl)} | ${fmt(pr)}`;
      const right = px(st.b[k]) < (l + r) / 2;
      const tt = svgEl("text", { x: px(st.b[k]) + (right ? 10 : -10), y: Math.max(t + 30, Math.min(bt - 10, py(cur.z) + (cur.z > zmax * 0.6 ? 26 : -12))), "text-anchor": right ? "start" : "end", "font-weight": 600 }, scan);
      tt.style.fill = OI.red; tt.style.paintOrder = "stroke"; tt.style.stroke = "#fff"; tt.style.strokeWidth = "4px";
      tt.textContent = txt;
    }
  }

  pBtns.forEach(b => b.classList.toggle("active", b.textContent === p0.label));
  draw();
}

// ================================================================ tab 2: FBA
function buildFba(box: HTMLElement, uid: number): void {
  void uid;
  const DEF: FBAOpts = { glc: 10, o2: 15, ko: RXN.map(() => false) };
  const st: FBAOpts & { fva: boolean } = { glc: DEF.glc, o2: DEF.o2, ko: DEF.ko.slice(), fva: true };

  const main = htmlEl("div", "lg-main", "", box);
  const left = htmlEl("div", "lg-fplot", "", main);
  const BW = 600, BH = 356;
  const bars = svgEl("svg", { viewBox: `0 0 ${BW} ${BH}`, role: "img", "aria-label": "Fluxes of the toy network" }, left);
  const OW = 600, OH = 200;
  const o2p = svgEl("svg", { viewBox: `0 0 ${OW} ${OH}`, role: "img", "aria-label": "Growth against the oxygen limit" }, left);
  const side = htmlEl("div", "lg-col lg-fside", "", main);

  const pB = presets(side, [
    { label: "wild type", apply: () => set(10, 15, []) },
    { label: "oxygen ≤ 4", apply: () => set(10, 4, []) },
    { label: "KO RES", apply: () => set(10, 15, [5]) },
    { label: "KO GLY + BYP1", apply: () => set(10, 15, [2, 3]) },
  ]);
  const gS = slider(side, "glucose uptake ≤", 0, 20, 0.5, st.glc, v => { st.glc = v; changed(); }, v => fmt(v, 1));
  const oS = slider(side, "oxygen uptake ≤", 0, 20, 0.5, st.o2, v => { st.o2 = v; changed(); }, v => fmt(v, 1));
  htmlEl("div", "lg-h", "knock out (flux forced to 0):", side);
  const koBox = htmlEl("div", "lg-ko", "", side);
  const kos = RXN.map((r, j) => checkbox(koBox, r, st.ko[j], v => { st.ko[j] = v; changed(); }));
  checkbox(side, "show flux variability (FVA) ranges", st.fva, v => { st.fva = v; draw(); });
  const ro = htmlEl("div", "qb-readout", "", side);

  function set(g: number, o: number, ko: number[]): void {
    st.glc = g; st.o2 = o; st.ko = RXN.map((_, j) => ko.includes(j));
    gS.set(g); oS.set(o); kos.forEach((c, j) => c.checked = st.ko[j]);
    draw();
  }
  function changed(): void { pB.forEach(b => b.classList.remove("active")); draw(); }

  const price = (which: "glc" | "o2"): [number, number] => {
    const h = 1e-4, z = (v: number) => { const o = { ...st, [which]: v }; const r = fba(o); return r.status === "optimal" ? r.z : NaN; };
    const z0 = z(st[which]);
    return [(z0 - z(st[which] - h)) / h, (z(st[which] + h) - z0) / h];
  };
  const pstr = ([l, r]: [number, number]) => Number.isNaN(l) || Math.abs(l - r) < 1e-6 ? fmt(r) : `${fmt(l)} below / ${fmt(r)} above`;

  function draw(): void {
    const sol = fba(st);
    const ranges = sol.status === "optimal" && st.fva ? fva(st, sol.z) : null;
    // ---- bars ----
    bars.replaceChildren();
    svgEl("text", { x: 6, y: 22, class: "ttl" }, bars).textContent = "fluxes at maximal growth" + (ranges ? " (band: FVA range)" : "");
    const l = 96, r = BW - 70, t = 40, rh = 30;
    const vmax = Math.max(20, ...(sol.x.length ? sol.x : [0]), ...(ranges ? ranges.map(x => Number.isFinite(x[1]) ? x[1] : 0) : [0]));
    const xm = Math.ceil(vmax / 5) * 5;
    const px = (v: number) => l + v / xm * (r - l);
    const bt = t + RXN.length * rh;
    for (let v = 0; v <= xm; v += 5) {
      svgEl("line", { x1: px(v), x2: px(v), y1: t - 4, y2: bt, class: v === 0 ? "axis" : "grid" }, bars);
      svgEl("text", { x: px(v), y: bt + 20, "text-anchor": "middle", class: "sm" }, bars).textContent = String(v);
    }
    RXN.forEach((name, j) => {
      const y = t + j * rh;
      const nm = svgEl("text", { x: l - 8, y: y + rh / 2 + 6, "text-anchor": "end", "font-weight": j === 8 ? 700 : 400 }, bars);
      nm.textContent = name;
      if (st.ko[j]) { nm.style.fill = OI.grey; nm.style.textDecoration = "line-through"; }
      if (ranges) {
        const [lo, hi] = ranges[j];
        const w = Math.max(3, px(Math.min(hi, xm)) - px(lo));
        svgEl("rect", { x: px(lo) - (hi - lo < 1e-6 ? 1.5 : 0), y: y + 3, width: w, height: rh - 6, fill: OI.sky, "fill-opacity": 0.45, rx: 3 }, bars);
      }
      if (sol.status === "optimal") {
        const v = sol.x[j];
        svgEl("rect", { x: px(0), y: y + rh / 2 - 6, width: Math.max(0, px(v) - px(0)), height: 12, fill: j === 8 ? OI.red : OI.blue }, bars);
        const tx = svgEl("text", { x: r + 6, y: y + rh / 2 + 6, "font-weight": j === 8 ? 700 : 400 }, bars);
        tx.textContent = ranges && ranges[j][1] - ranges[j][0] > 1e-6 ? `${fmt(v, 2)}*` : fmt(v, 2);
        if (st.ko[j]) tx.textContent = "KO";
      }
    });
    svgEl("text", { x: (l + r) / 2, y: BH - 4, "text-anchor": "middle", class: "sm" }, bars).textContent = "flux (toy units; * = not determined)";

    // ---- oxygen scan ----
    o2p.replaceChildren();
    svgEl("text", { x: 6, y: 20, class: "ttl" }, o2p).textContent = "growth (red), acetate (orange) vs oxygen limit";
    const ol = 60, or = OW - 14, ot = 32, ob = OH - 42;
    const N = 81, os: number[] = [], gs: number[] = [], as: number[] = [];
    for (let i = 0; i < N; i++) {
      const o = 20 * i / (N - 1), s = fba({ ...st, o2: o });
      os.push(o); gs.push(s.status === "optimal" ? s.z : NaN); as.push(s.status === "optimal" ? s.x[7] : NaN);
    }
    // acetate at minimal-flux is not unique in general; here it is determined by growth
    const ym = Math.max(5, ...gs.filter(Number.isFinite), ...as.filter(Number.isFinite)) * 1.08;
    const qx = (v: number) => ol + v / 20 * (or - ol), qy = (v: number) => ob - v / ym * (ob - ot);
    svgEl("line", { x1: ol, x2: or, y1: ob, y2: ob, class: "axis" }, o2p);
    svgEl("line", { x1: ol, x2: ol, y1: ot, y2: ob, class: "axis" }, o2p);
    for (let v = 0; v <= 20; v += 2) {
      svgEl("text", { x: qx(v), y: ob + 20, "text-anchor": "middle", class: "sm" }, o2p).textContent = String(v);
    }
    const ys = ym > 40 ? 10 : 5;
    for (let v = 0; v <= ym; v += ys) {
      svgEl("text", { x: ol - 6, y: qy(v) + 6, "text-anchor": "end", class: "sm" }, o2p).textContent = String(v);
      if (v > 0) svgEl("line", { x1: ol, x2: or, y1: qy(v), y2: qy(v), class: "grid" }, o2p);
    }
    svgEl("text", { x: (ol + or) / 2, y: OH - 4, "text-anchor": "middle", class: "sm" }, o2p).textContent = "oxygen uptake limit";
    const path = (ys2: number[], col: string) => {
      let d = "", pen = false;
      ys2.forEach((y, i) => { if (!Number.isFinite(y)) { pen = false; return; } d += `${pen ? "L" : "M"}${qx(os[i]).toFixed(1)},${qy(y).toFixed(1)} `; pen = true; });
      if (d) svgEl("path", { d, fill: "none", stroke: col, "stroke-width": 3 }, o2p);
    };
    path(as, OI.orange); path(gs, OI.red);
    svgEl("line", { x1: qx(st.o2), x2: qx(st.o2), y1: ot, y2: ob, stroke: OI.ink, "stroke-dasharray": "3 4" }, o2p);
    if (sol.status === "optimal") {
      svgEl("circle", { cx: qx(st.o2), cy: qy(sol.z), r: 7, fill: OI.red, stroke: "#fff", "stroke-width": 2 }, o2p);
      svgEl("circle", { cx: qx(st.o2), cy: qy(sol.x[7]), r: 6, fill: OI.orange, stroke: "#fff", "stroke-width": 2 }, o2p);
    }

    // ---- readout ----
    ro.replaceChildren();
    const line = (html: string, cls = "") => { const d = htmlEl("div", cls, "", ro); d.innerHTML = html; };
    if (sol.status !== "optimal") { line(`<b>${sol.status}</b>`, "lg-bad"); return; }
    line(`<b>growth (BIO) ${fmt(sol.z)}</b>${sol.z < 1e-9 ? " — lethal" : ""}`, "lg-ok");
    line(`glucose used ${fmt(sol.x[0])}, oxygen used ${fmt(sol.x[1])}, acetate secreted ${fmt(sol.x[7])}`);
    line(`shadow prices: glucose <b>${pstr(price("glc"))}</b>, oxygen <b>${pstr(price("o2"))}</b>`);
    line(`max |S v| = ${sol.x.length ? maxResidual(sol.x).toExponential(1) : "–"} (steady state)`);
    if (ranges) {
      const free = RXN.filter((_, j) => ranges[j][1] - ranges[j][0] > 1e-6);
      line(free.length ? `not determined at the optimum: ${free.map(n => `${n} [${fmt(ranges[RXN.indexOf(n)][0], 2)}, ${fmt(ranges[RXN.indexOf(n)][1], 2)}]`).join(", ")}` : "every flux is determined");
    }
  }
  pB[0].classList.add("active");
  draw();
}

autoMount(NAME, mount);

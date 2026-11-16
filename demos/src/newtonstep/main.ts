// Newton / secant / bisection step-through for the newton unit.
// Three equations from the lecture: the two-site dose-response f(c) = R(c) - 0.5 (root
// c50 = 3.2613), the zero-truncated Poisson score S(lambda) (root lambda-hat = 1.2185),
// and Michaelis-Menten in log substrate sigma(u) - 1/2 (root u = 0, Newton step u - sinh u).
// Drag (or click) in the plot to move the start value(s); Step takes one step and draws the
// tangent (Newton), the chord (secant) or the bracket (bisection) and where it hits zero.
// The table lists the iterates with f, the error against the wanted root, the number of
// correct digits -log10(|e| / max(|root|, 1)) and the factor by which the digits grew:
// about 2 for Newton, 1.6 for secant, while bisection adds 0.3 digits per step.
// Presets reproduce the lecture's runs and its three failures (flat tangent -> negative
// dose; bad start for the MLE runs off; cycling at S0 = 8.82 Km).
// Auto-mounts on <div data-demo="newtonstep">. Optional attributes:
//   data-problem="dose|score|mm"  data-method="newton|secant|bisection"
//   data-x0="1"  data-x1="2"  (start values; x1 only for secant / bisection)
//   data-steps="0"  (steps already taken when the demo appears)
import { svgEl, htmlEl, button, presets, scale, svgPoint, autoMount } from "../lib/ui";
import { PROBLEMS, Problem, ProbKey, Method, Run, start, step, steps, errDigits, sub, MAXIT } from "./math";

const NAME = "newtonstep";
// Okabe-Ito, as in make-newton-figs.R: curve blue, Newton red, secant green, bisection grey
const COL = { curve: "#0072B2", newton: "#D55E00", secant: "#009E73", bisection: "#6b6b6b", root: "#E69F00", ink: "#222", grey: "#999999" };
const METHOD_LABEL: Record<Method, string> = { newton: "Newton", secant: "secant", bisection: "bisection" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .ns-row { display: flex; flex-wrap: wrap; gap: .35em 1.2em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .ns-group { display: flex; flex-wrap: wrap; gap: .35em; align-items: center; }
[data-demo="${NAME}"] .ns-main { display: flex; flex-wrap: wrap; gap: 12px; width: 100%; align-items: flex-start; }
[data-demo="${NAME}"] .ns-main svg { flex: 1 1 520px; min-width: 0; max-width: 680px; height: auto; touch-action: none; cursor: crosshair; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] svg .hl { font-size: 18px; font-weight: 600; }
[data-demo="${NAME}"] .ns-side { flex: 1 1 400px; min-width: 0; display: flex; flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .ns-status { min-height: 2.6em; }
[data-demo="${NAME}"] .ns-status.warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .ns-status.good { color: #006b4f; font-weight: 600; }
[data-demo="${NAME}"] .ns-tab { max-height: 330px; overflow-y: auto; border: 1px solid #ddd; border-radius: 6px; }
[data-demo="${NAME}"] table { border-collapse: collapse; width: 100%; font-family: ui-monospace, monospace; font-size: 15px; }
[data-demo="${NAME}"] th, [data-demo="${NAME}"] td { padding: 2px 6px; text-align: right; white-space: nowrap; }
[data-demo="${NAME}"] th { position: sticky; top: 0; background: #f0f1f4; font-family: system-ui, sans-serif; font-weight: 600; }
[data-demo="${NAME}"] tr.last td { background: #fff3e0; }
[data-demo="${NAME}"] td .ok { color: #006b4f; font-weight: 700; }
[data-demo="${NAME}"] td .rest { color: #999; }
[data-demo="${NAME}"] .qb-btn.go { font-weight: 600; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; gap: .3em; }
.reveal [data-demo="${NAME}"] .ns-main { flex-wrap: nowrap; }
.reveal [data-demo="${NAME}"] .ns-main svg { flex: 0 0 600px; max-width: 600px; }
.reveal [data-demo="${NAME}"] .ns-side { flex: 1 1 auto; }
.reveal [data-demo="${NAME}"] .ns-tab { max-height: 300px; }
.reveal [data-demo="${NAME}"] table { font-size: 18px; }
.reveal [data-demo="${NAME}"] th, .reveal [data-demo="${NAME}"] td { padding: 1px 5px; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 18px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- formatting ----
function full(v: number): string {
  if (!Number.isFinite(v)) return String(v);
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a < 1e-3 || a >= 1e6) return v.toExponential(10).replace("e+", "e");
  return v.toPrecision(16);
}
function sci(v: number, d = 1): string {
  if (!Number.isFinite(v)) return String(v);
  if (v === 0) return "0";
  return v.toExponential(d).replace("e+", "e");
}
const short = (v: number): string => !Number.isFinite(v) ? String(v)
  : Math.abs(v) >= 1e4 ? sci(v, 2) : String(+v.toPrecision(4));

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const uid = ++mountCount;
  const ds = host.dataset;
  const pk0: ProbKey = ds.problem === "score" || ds.problem === "mm" ? ds.problem : "dose";
  const m0: Method = ds.method === "secant" || ds.method === "bisection" ? ds.method : "newton";
  const num = (s: string | undefined): number | undefined => s !== undefined && s.trim() !== "" && Number.isFinite(Number(s)) ? Number(s) : undefined;

  const st = {
    pk: pk0, m: m0,
    pts: {} as Record<Method, number[]>,
    run: null as unknown as Run,
    timer: 0,
  };
  const loadDefaults = (pk: ProbKey) => {
    const d = PROBLEMS[pk].defaults;
    st.pts = { newton: [...d.newton], secant: [...d.secant], bisection: [...d.bisection] };
  };
  loadDefaults(pk0);
  const x0a = num(ds.x0), x1a = num(ds.x1);
  if (x0a !== undefined) st.pts[m0][0] = x0a;
  if (x1a !== undefined && m0 !== "newton") st.pts[m0][1] = x1a;
  const P = () => PROBLEMS[st.pk];

  // ---- controls ----
  const row1 = htmlEl("div", "ns-row", "", host);
  const g1 = htmlEl("div", "ns-group", "", row1);
  htmlEl("span", "", "equation:", g1);
  const pKeys: ProbKey[] = ["dose", "score", "mm"];
  const pBtns = presets(g1, pKeys.map(k => ({ label: PROBLEMS[k].label, apply: () => { st.pk = k; loadDefaults(k); restart(); } })));
  const g2 = htmlEl("div", "ns-group", "", row1);
  htmlEl("span", "", "method:", g2);
  const mKeys: Method[] = ["newton", "secant", "bisection"];
  const mBtns = presets(g2, mKeys.map(k => ({ label: METHOD_LABEL[k], apply: () => { st.m = k; restart(); } })));

  const row2 = htmlEl("div", "ns-row", "", host);
  const g3 = htmlEl("div", "ns-group", "", row2);
  button(g3, "Step ▸", () => { stop(); step(st.run); draw(); }, "go");
  const runBtn = button(g3, "Run", () => { if (st.timer) stop(); else play(); }, "go");
  button(g3, "Back", () => { stop(); back(); });
  button(g3, "Reset", () => { stop(); st.run = start(P(), st.m, st.pts[st.m]); draw(); });
  const g4 = htmlEl("div", "ns-group", "", row2);
  htmlEl("span", "", "lecture:", g4);
  const pre = (label: string, pk: ProbKey, m: Method, pts: number[]) => ({
    label, apply: () => { st.pk = pk; loadDefaults(pk); st.m = m; st.pts[m] = [...pts]; restart(); },
  });
  const preBtns = presets(g4, [
    pre("Newton c₀ = 1", "dose", "newton", [1]),
    pre("secant 1, 2", "dose", "secant", [1, 2]),
    pre("bisection [1, 10]", "dose", "bisection", [1, 10]),
    pre("MLE from x̄ = 1.73", "score", "newton", [1.73]),
  ]);
  const g5 = htmlEl("div", "ns-group", "", row2);
  htmlEl("span", "", "failures:", g5);
  const failBtns = presets(g5, [
    pre("1 flat tangent", "dose", "newton", [10]),
    pre("2 bad start", "score", "newton", [3]),
    pre("3 cycling", "mm", "newton", [2.1773189849653]),
  ]);
  const allPre = [...preBtns, ...failBtns];
  allPre.forEach(b => b.addEventListener("click", () => allPre.forEach(x => x.classList.toggle("active", x === b))));
  const clearPre = () => allPre.forEach(x => x.classList.remove("active"));
  pBtns.forEach(b => b.addEventListener("click", clearPre));
  mBtns.forEach(b => b.addEventListener("click", clearPre));

  // ---- plot and table ----
  const main = htmlEl("div", "ns-main", "", host);
  const W = 640, H = 430;
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "The function, the start values and the Newton, secant or bisection steps" }, main);
  const side = htmlEl("div", "ns-side", "", main);
  const info = htmlEl("div", "qb-readout", "", side);
  const statusEl = htmlEl("div", "ns-status", "", side);
  const tabBox = htmlEl("div", "ns-tab", "", side);
  const table = htmlEl("table", "", "", tabBox);

  const L = 74, Rm = 14, T = 16, B = H - 52;

  // ---- dragging: the pointer moves the nearest start value ----
  let dragIdx = -1;
  const toX = (ev: PointerEvent) => scale(P().xr[0], P().xr[1], L, W - Rm).invert(svgPoint(svg, ev)[0]);
  const snap = (x: number) => {
    const p = P(), r = p.res;
    let v = Math.round(x / r) * r;
    v = Math.max(p.xr[0], Math.min(p.xr[1], v));
    return +v.toFixed(6);
  };
  const setPt = (i: number, x: number) => {
    const pts = st.pts[st.m];
    let v = snap(x);
    if (st.m !== "newton" && v === pts[1 - i]) v = +(v + P().res * (i === 0 ? -1 : 1)).toFixed(6);
    pts[i] = v;
    st.run = start(P(), st.m, pts);
    clearPre();
    draw();
  };
  svg.addEventListener("pointerdown", ev => {
    ev.preventDefault();
    stop();
    const x = toX(ev), pts = st.pts[st.m];
    dragIdx = pts.length === 1 ? 0 : (Math.abs(x - pts[0]) <= Math.abs(x - pts[1]) ? 0 : 1);
    svg.setPointerCapture(ev.pointerId);
    setPt(dragIdx, x);
  });
  svg.addEventListener("pointermove", ev => { if (dragIdx >= 0) setPt(dragIdx, toX(ev)); });
  const endDrag = () => { dragIdx = -1; };
  svg.addEventListener("pointerup", endDrag);
  svg.addEventListener("pointercancel", endDrag);

  // ---- run control ----
  function stop(): void {
    if (st.timer) { clearInterval(st.timer); st.timer = 0; }
    runBtn.textContent = "Run";
  }
  function play(): void {
    if (st.run.done) return;
    runBtn.textContent = "Pause";
    const wasCycling = st.run.cycling;
    let cyc = wasCycling;
    const tick = () => {
      step(st.run);
      // pause the first time a cycle is detected, so it can be looked at
      const stopNow = st.run.done || (st.run.cycling && !cyc);
      cyc = st.run.cycling;
      draw();
      if (stopNow) stop();
    };
    tick();
    if (st.timer === 0 && runBtn.textContent === "Run") return;
    st.timer = window.setInterval(tick, st.m === "bisection" ? 220 : 650);
  }
  function back(): void {
    const k = steps(st.run);
    if (k === 0) return;
    const r = start(P(), st.m, st.pts[st.m]);
    for (let i = 0; i < k - 1; i++) step(r);
    st.run = r;
    draw();
  }
  function restart(): void {
    stop();
    st.run = start(P(), st.m, st.pts[st.m]);
    draw();
  }

  // ---- drawing ----
  function draw(): void {
    const p = P(), run = st.run, m = st.m;
    pBtns.forEach((b, i) => b.classList.toggle("active", pKeys[i] === st.pk));
    mBtns.forEach((b, i) => b.classList.toggle("active", mKeys[i] === m));
    const sx = scale(p.xr[0], p.xr[1], L, W - Rm), sy = scale(p.yr[0], p.yr[1], B, T);
    const ySpan = p.yr[1] - p.yr[0];
    const cy = (y: number) => sy(Math.max(p.yr[0] - ySpan, Math.min(p.yr[1] + ySpan, y)));
    const col = COL[m];
    svg.replaceChildren();

    // axes and grid
    for (const v of p.yticks) {
      svgEl("line", { x1: L, x2: W - Rm, y1: sy(v), y2: sy(v), class: "grid" }, svg);
      svgEl("text", { x: L - 6, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = String(v);
    }
    for (const v of p.xticks) {
      svgEl("line", { x1: sx(v), x2: sx(v), y1: B, y2: B + 5, class: "axis" }, svg);
      svgEl("text", { x: sx(v), y: B + 23, "text-anchor": "middle", class: "sm" }, svg).textContent = String(v);
    }
    svgEl("line", { x1: L, x2: W - Rm, y1: B, y2: B, class: "axis" }, svg);
    svgEl("text", { x: (L + W - Rm) / 2, y: H - 6, "text-anchor": "middle" }, svg).textContent = p.xlab;
    const yl = svgEl("text", { x: 18, y: (T + B) / 2, "text-anchor": "middle", transform: `rotate(-90 18 ${(T + B) / 2})` }, svg);
    yl.textContent = p.fname;

    const clipId = `ns-clip-${uid}`;
    const cp = svgEl("clipPath", { id: clipId }, svg);
    svgEl("rect", { x: L, y: T, width: W - Rm - L, height: B - T }, cp);
    const g = svgEl("g", { "clip-path": `url(#${clipId})` }, svg);

    if (p.positive) {
      svgEl("rect", { x: L, y: T, width: sx(0) - L, height: B - T, fill: COL.grey, "fill-opacity": 0.18 }, g);
      svgEl("text", { x: L + 8, y: T + 22, class: "sm" }, g).textContent = `negative ${p.positive}`;
    }
    // zero line
    svgEl("line", { x1: L, x2: W - Rm, y1: sy(0), y2: sy(0), stroke: COL.ink, "stroke-width": 1.3 }, g);

    // the curve, broken at poles
    const n = 900;
    let seg: string[] = [], prevY = NaN;
    const flush = () => { if (seg.length > 1) svgEl("polyline", { points: seg.join(" "), fill: "none", stroke: COL.curve, "stroke-width": 3 }, g); seg = []; };
    for (let i = 0; i <= n; i++) {
      const x = p.xr[0] + (p.xr[1] - p.xr[0]) * i / n, y = p.f(x);
      if (!Number.isFinite(y) || (Number.isFinite(prevY) && Math.sign(y) !== Math.sign(prevY) && Math.abs(y - prevY) > ySpan)) flush();
      if (Number.isFinite(y)) seg.push(`${sx(x).toFixed(1)},${cy(y).toFixed(1)}`);
      prevY = y;
    }
    flush();

    // the wanted root
    svgEl("circle", { cx: sx(p.root), cy: sy(0), r: 8, fill: COL.root, stroke: COL.ink, "stroke-width": 1.5 }, g);
    const rl = svgEl("text", { x: sx(p.root) + 10, y: sy(0) + 24, class: "hl" }, g);
    rl.style.fill = "#9a6700"; rl.textContent = p.rootName;

    // steps: earlier ones faded, the last one strong
    const stepRows = run.rows.filter(r => r.k > 0);
    stepRows.forEach((r, i) => {
      const last = i === stepRows.length - 1;
      const op = last ? 1 : 0.4, w = last ? 2.8 : 2;
      if (!Number.isFinite(r.x)) return;
      if (m === "newton") {
        const x0 = r.from[0], f0 = p.f(x0);
        svgEl("line", { x1: sx(x0), x2: sx(x0), y1: sy(0), y2: cy(f0), stroke: col, "stroke-width": 2, "stroke-dasharray": "3 4", "stroke-opacity": op }, g);
        const dx = r.x - x0, xa = x0 - 0.25 * dx, xb = r.x + 0.25 * dx;
        const slope = (0 - f0) / dx;
        svgEl("line", { x1: sx(xa), y1: cy(f0 + slope * (xa - x0)), x2: sx(xb), y2: cy(slope * (xb - r.x)), stroke: col, "stroke-width": w, "stroke-opacity": op }, g);
        svgEl("circle", { cx: sx(x0), cy: cy(f0), r: last ? 6 : 4.5, fill: col, "fill-opacity": op }, g);
      } else if (m === "secant") {
        const [xa0, xb0] = r.from, fa = p.f(xa0), fb = p.f(xb0);
        const xs = [xa0, xb0, r.x], lo = Math.min(...xs), hi = Math.max(...xs), ext = 0.12 * (hi - lo);
        const slope = (fb - fa) / (xb0 - xa0);
        const yAt = (x: number) => fb + slope * (x - xb0);
        svgEl("line", { x1: sx(lo - ext), y1: cy(yAt(lo - ext)), x2: sx(hi + ext), y2: cy(yAt(hi + ext)), stroke: col, "stroke-width": w, "stroke-opacity": op }, g);
        for (const [x, y] of [[xa0, fa], [xb0, fb]]) svgEl("circle", { cx: sx(x), cy: cy(y), r: last ? 6 : 4.5, fill: col, "fill-opacity": op }, g);
        svgEl("line", { x1: sx(r.x), x2: sx(r.x), y1: sy(0), y2: cy(r.fx), stroke: col, "stroke-width": 2, "stroke-dasharray": "3 4", "stroke-opacity": op }, g);
      } else if (last) {
        // bracket before this step, the midpoint, and the kept half
        const [a, b] = r.from, yb = B - 18;
        svgEl("line", { x1: sx(a), x2: sx(b), y1: yb, y2: yb, stroke: COL.grey, "stroke-width": 6, "stroke-opacity": 0.6 }, g);
        svgEl("line", { x1: sx(r.a!), x2: sx(r.b!), y1: yb, y2: yb, stroke: col, "stroke-width": 6 }, g);
        for (const x of [r.a!, r.b!]) svgEl("line", { x1: sx(x), x2: sx(x), y1: yb - 12, y2: yb + 12, stroke: col, "stroke-width": 2.5 }, g);
        svgEl("line", { x1: sx(r.x), x2: sx(r.x), y1: sy(0), y2: cy(r.fx), stroke: col, "stroke-width": 2, "stroke-dasharray": "3 4" }, g);
        svgEl("circle", { cx: sx(r.x), cy: cy(r.fx), r: 6, fill: col }, g);
        svgEl("text", { x: sx((r.a! + r.b!) / 2), y: yb - 16, "text-anchor": "middle", class: "sm" }, g).textContent = "bracket";
      }
    });
    // iterates on the zero line, labelled when there is room
    const labelled: number[] = [];
    const lab = (x: number, txt: string, c: string, below: boolean) => {
      const px = sx(x);
      if (labelled.some(q => Math.abs(q - px) < 34)) return;
      labelled.push(px);
      const t = svgEl("text", { x: px, y: sy(0) + (below ? 26 : -12), "text-anchor": "middle", class: "hl" }, g);
      t.style.fill = c; t.textContent = txt;
    };
    for (let i = stepRows.length - 1; i >= 0; i--) {
      const r = stepRows[i];
      if (!Number.isFinite(r.x)) continue;
      svgEl("circle", { cx: sx(r.x), cy: sy(0), r: i === stepRows.length - 1 ? 6.5 : 4.5, fill: col, stroke: "#fff", "stroke-width": 1.2 }, g);
      if (m !== "bisection" || i === stepRows.length - 1) lab(r.x, `${p.sym}${sub(r.k + (m === "secant" ? 1 : 0))}`, col, i % 2 === 1);
    }

    // start handles (drag or click in the plot)
    const pts = st.pts[m];
    const names = m === "bisection" ? ["a", "b"] : pts.map((_, i) => `${p.sym}${sub(i)}`);
    pts.forEach((x, i) => {
      const hx = sx(x), hy = sy(0);
      svgEl("line", { x1: hx, x2: hx, y1: hy, y2: cy(p.f(x)), stroke: COL.ink, "stroke-width": 1.2, "stroke-dasharray": "2 3" }, g);
      svgEl("circle", { cx: hx, cy: hy, r: 12, fill: "#fff", stroke: col, "stroke-width": 3.5 }, g);
      const t = svgEl("text", { x: hx, y: hy + 6, "text-anchor": "middle", "font-weight": 700 }, g);
      t.style.fill = col; t.style.fontSize = "15px"; t.textContent = names[i];
    });

    // iterates outside the window: an arrow at the edge
    const lastRow = stepRows[stepRows.length - 1];
    if (lastRow && (lastRow.x < p.xr[0] || lastRow.x > p.xr[1] || !Number.isFinite(lastRow.x))) {
      const left = !(lastRow.x > p.xr[1]);
      const t = svgEl("text", { x: left ? L + 8 : W - Rm - 8, y: T + 50, "text-anchor": left ? "start" : "end", class: "hl" }, svg);
      t.style.fill = col;
      t.textContent = `${left ? "◀ " : ""}${p.sym}${sub(lastRow.k + (m === "secant" ? 1 : 0))} = ${short(lastRow.x)}${left ? "" : " ▶"}`;
    }
    svgEl("text", { x: W - Rm - 8, y: T + 22, "text-anchor": "end", class: "sm" }, svg).textContent = "drag or click to move the start";

    // ---- info, status, table ----
    const k = stepRows.length;
    const startTxt = m === "newton" ? `${p.sym}₀ = ${pts[0]}` : m === "secant" ? `${p.sym}₀ = ${pts[0]}, ${p.sym}₁ = ${pts[1]}` : `[a, b] = [${Math.min(...pts)}, ${Math.max(...pts)}]`;
    const mmTxt = p.key === "mm" && m === "newton" ? `  (S₀ = ${Math.exp(pts[0]).toPrecision(3)} Km)` : "";
    info.textContent = `${METHOD_LABEL[m]} from ${startTxt}${mmTxt}\nsteps ${k}   evaluations ${run.evals}` +
      (m === "newton" ? " (f and f′)" : " (f only)");
    statusEl.className = "ns-status" + (run.level ? " " + run.level : "");
    statusEl.textContent = run.status || (k === 0 ? (run.done ? "" : "Press Step or Run.") :
      (m === "bisection" ? "Each step halves the bracket: about 0.3 digits per step." :
        m === "newton" ? "Near the root the digits roughly double each step." : "Near the root the digits grow by about 1.6× per step."));

    table.replaceChildren();
    const hr = htmlEl("tr", "", "", htmlEl("thead", "", "", table));
    for (const h of ["step", p.sym + "ₖ", "f(" + p.sym + "ₖ)", "|error|", "digits", m === "bisection" ? "+ digits" : "× digits"]) htmlEl("th", "", h, hr);
    const tb = htmlEl("tbody", "", "", table);
    let prevD = NaN;
    const rootStr = full(p.root);
    run.rows.forEach((r, i) => {
      const tr = htmlEl("tr", i === run.rows.length - 1 && r.k > 0 ? "last" : "", "", tb);
      const kTxt = r.k === 0 ? (m === "bisection" ? (i === 0 ? "a" : "b") : `${p.sym}${sub(i)}`) : String(r.k);
      htmlEl("td", "", kTxt, tr);
      // x with the digits that agree with the root in bold
      const td = htmlEl("td", "", "", tr);
      const xs = full(r.x);
      let nm = 0;
      if (p.root !== 0 && Number.isFinite(r.x)) { while (nm < xs.length && nm < rootStr.length && xs[nm] === rootStr[nm]) nm++; }
      if (nm > 0 && xs.length === rootStr.length && nm < xs.length && /[.\d]/.test(xs[nm - 1])) {
        htmlEl("span", "ok", xs.slice(0, nm), td); htmlEl("span", "rest", xs.slice(nm), td);
      } else if (nm === xs.length && nm > 0) htmlEl("span", "ok", xs, td);
      else td.textContent = xs;
      htmlEl("td", "", sci(r.fx, 2), tr);
      const { e, d } = errDigits(p, r.x);
      htmlEl("td", "", Number.isFinite(e) ? sci(e, 1) : "–", tr);
      const dTxt = !Number.isFinite(d) ? (e === 0 ? "all" : "–") : d >= 15.5 ? "all" : Math.max(0, d).toFixed(1);
      htmlEl("td", "", dTxt, tr);
      const ratio = r.k === 0 || (m === "bisection" && r.k === 1) || !Number.isFinite(d) || d >= 15.5 || !Number.isFinite(prevD) ? ""
        : m === "bisection" ? ((d - prevD) >= 0 ? "+" : "−") + Math.abs(d - prevD).toFixed(1)
          : prevD >= 0.5 ? (d / prevD).toFixed(1) : "";
      htmlEl("td", "", ratio, tr);
      prevD = d;
    });
    tabBox.scrollTop = tabBox.scrollHeight;
  }

  st.run = start(P(), st.m, st.pts[st.m]);
  const pre0 = Math.max(0, Math.min(MAXIT, Math.round(num(ds.steps) ?? 0)));
  for (let i = 0; i < pre0; i++) step(st.run);
  draw();
}

autoMount(NAME, mount);

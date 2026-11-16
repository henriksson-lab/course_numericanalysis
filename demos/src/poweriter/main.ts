// Power iteration explorer for the power-iteration unit. x <- A x / ||x|| for a matrix from
// the lecture (2 x 2 plants, the fish Leslie matrix, its semelparous variant, the cell-state
// Markov chain) or a 2 x 2 matrix set with sliders.
//  (a) 2 x 2: the plane with the eigen-directions (dashed), the unit circle, the iterates
//      x_0 .. x_k swinging towards v_1, and y = A x_k; the start x_0 is a drag handle.
//      n > 2: bars of x_k (proportions when v_1 is positive) against v_1 (dashed outlines).
//  (b) error in x and in the eigenvalue estimate against step k on a log axis, with the
//      predicted line C |lambda_2/lambda_1|^k (or the inverse-iteration rate).
// Inverse iteration with a shift sigma is a toggle. All numbers are computed live in doubles.
// Auto-mounts on <div data-demo="poweriter">. Optional attributes:
//   data-preset="plant|plantb1|biennial|early|sym|fish|semel|cells|slow"
//   data-method="power|inverse"  data-sigma="1.4"  data-k="0"  data-kmax="60"
//   data-edit="1|0" (show the 2 x 2 entry sliders)
import { svgEl, htmlEl, slider, button, presets, scale, svgPoint, draggable, autoMount } from "../lib/ui";
import { Mat, C, analyse, iterate, vecErr, matvec, norm2, sum, cabs } from "./math";

const NAME = "poweriter";
// Okabe-Ito, as in make-power-iteration-figs.R
const OI = { orange: "#E69F00", sky: "#56B4E9", green: "#009E73", blue: "#0072B2", red: "#D55E00", purple: "#CC79A7", grey: "#999999" };
const COL = { x: OI.blue, lam: OI.orange, res: OI.green, pred: "#737373", v1: OI.blue, v2: OI.red, y: OI.green, x0: OI.purple, ink: "#222" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .pi-row { display: flex; flex-wrap: wrap; gap: .35em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .pi-group { display: flex; gap: .35em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .pi-row .qb-slider { min-width: 210px; flex: 1 1 210px; max-width: 360px; grid-template-columns: auto 1fr 3.4em; }
[data-demo="${NAME}"] .pi-ent .qb-slider { min-width: 170px; flex: 1 1 170px; max-width: 260px; }
[data-demo="${NAME}"] .pi-plots { display: flex; flex-wrap: wrap; gap: 10px; width: 100%; }
[data-demo="${NAME}"] .pi-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .pi-read { display: flex; flex-wrap: wrap; gap: .15em 1.3em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .pi-read b { font-weight: 600; }
[data-demo="${NAME}"] .pi-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .pi-off { opacity: .35; pointer-events: none; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 18px; gap: .25em; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; padding: .3em .5em; }
.reveal [data-demo="${NAME}"] .qb-btn { font-size: .85em; padding: .15em .55em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- the lecture's matrices ----
interface Preset { id: string; label: string; A: Mat; x0: number[]; names: string[]; sigma: number; sigmas: number[]; cols?: string[] }
const fishL = (F: number[]): Mat => [F, [0.3, 0, 0, 0], [0, 0.6, 0, 0], [0, 0, 0.5, 0]];
const P: Mat = [[0.80, 0.15, 0.05], [0.05, 0.90, 0.05], [0.05, 0.10, 0.85]];
const T = (M: Mat): Mat => M[0].map((_, j) => M.map(r => r[j]));
const Pslow: Mat = P.map((r, i) => r.map((v, j) => (i === j ? 1 : 0) + 0.01 * (v - (i === j ? 1 : 0))));
const AGE = ["age 0", "age 1", "age 2", "age 3"], AGECOL = [OI.blue, OI.orange, OI.green, OI.purple];
const CELL = ["S", "B", "L"], CELLCOL = [OI.red, OI.blue, OI.green];
const PRESETS: Preset[] = [
  { id: "plant", label: "plant (1.5, −0.5)", A: [[1, 1], [0.75, 0]], x0: [1, 0], names: ["class 1", "class 2"], sigma: 1.4, sigmas: [1.4, 0, -0.4, 0.5] },
  { id: "plantb1", label: "plant B1", A: [[0.6, 0.9], [0.8, 0]], x0: [0, 1], names: ["class 1", "class 2"], sigma: 1, sigmas: [1, 0, -0.5] },
  { id: "biennial", label: "biennial", A: [[0, 4], [0.25, 0]], x0: [1, 0], names: ["rosettes", "flowering"], sigma: 0.5, sigmas: [0.5, 1.2] },
  { id: "early", label: "early flowering", A: [[0.4, 4], [0.25, 0]], x0: [1, 0], names: ["rosettes", "flowering"], sigma: 1, sigmas: [1, -0.8] },
  { id: "sym", label: "symmetric (6, 1)", A: [[5, 2], [2, 2]], x0: [1, 0], names: ["trait 1", "trait 2"], sigma: 5, sigmas: [5, 1.2, 3.5] },
  { id: "fish", label: "fish", A: fishL([0, 1.5, 4, 5]), x0: [0, 1, 0, 0], names: AGE, sigma: 1.1, sigmas: [1.1, -0.3, -0.7], cols: AGECOL },
  { id: "semel", label: "semelparous fish", A: fishL([0, 0, 0, 12]), x0: [0, 1, 0, 0], names: AGE, sigma: 1.1, sigmas: [1.1, -1], cols: AGECOL },
  { id: "cells", label: "cell states Pᵀ", A: T(P), x0: [1, 0, 0], names: CELL, sigma: 0.82, sigmas: [1, 0.82, 0.77, 0.9], cols: CELLCOL },
  { id: "slow", label: "slow cells", A: T(Pslow), x0: [1, 0, 0], names: CELL, sigma: 1.005, sigmas: [1.005, 1], cols: CELLCOL },
];

// ---- formatting ----
const sup: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const pow10 = (k: number): string => k === 0 ? "1" : "10" + String(k).split("").map(c => sup[c]).join("");
function sci(v: number, d = 2): string {
  if (!Number.isFinite(v)) return "–";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e-3 && a < 1e5) return v.toPrecision(d);
  const [m, e] = v.toExponential(d - 1).split("e");
  return `${m}e${Number(e)}`;
}
const minus = (s: string): string => s.replace(/-/g, "−");
const fnum = (v: number, d = 4): string => minus((Math.abs(v) < 0.5 * 10 ** -d ? 0 : v).toFixed(d));
const fC = (z: C): string => z.im === 0 ? fnum(z.re) : `${fnum(z.re)} ${z.im > 0 ? "+" : "−"} ${fnum(Math.abs(z.im))}i`;
function eigList(e: C[]): string {
  const out: string[] = [];
  for (let i = 0; i < e.length; i++) {
    const z = e[i], w = e[i + 1];
    if (z.im !== 0 && w && Math.abs(w.re - z.re) < 1e-12 && Math.abs(w.im + z.im) < 1e-12) {
      out.push(`${fnum(z.re)} ± ${fnum(Math.abs(z.im))}i`); i++;
    } else out.push(fC(z));
  }
  return out.join(",  ");
}

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const uid = ++mountCount;
  const KMAX = Math.max(5, Math.min(200, Math.round(Number(ds.kmax) || 60)));
  const start = {
    preset: PRESETS.find(p => p.id === ds.preset) ?? PRESETS[0],
    inv: ds.method === "inverse",
    sigma: Number.isFinite(Number(ds.sigma)) && ds.sigma !== undefined ? Number(ds.sigma) : NaN,
    k: Math.max(0, Math.min(KMAX, Math.round(Number(ds.k) || 0))),
  };
  const st = {
    A: start.preset.A.map(r => [...r]) as Mat, x0: [...start.preset.x0], names: start.preset.names,
    cols: start.preset.cols, inv: start.inv,
    sigma: Number.isFinite(start.sigma) ? start.sigma : start.preset.sigma, k: start.k,
  };
  const editable = ds.edit !== "0";

  // ---- controls ----
  const row1 = htmlEl("div", "pi-row", "", host);
  const pg = htmlEl("div", "pi-group", "", row1);
  htmlEl("span", "", "matrix:", pg);
  const pBtns = presets(pg, PRESETS.map(p => ({ label: p.label, apply: () => loadPreset(p, true) })));

  const row2 = htmlEl("div", "pi-row", "", host);
  const mg = htmlEl("div", "pi-group", "", row2);
  const mBtns = presets(mg, [
    { label: "power iteration", apply: () => { st.inv = false; update(); } },
    { label: "inverse iteration", apply: () => { st.inv = true; update(); } },
  ]);
  const sS = slider(row2, "shift σ", -1.5, 6, 0.005, st.sigma, v => { st.sigma = v; update(); }, v => minus(v.toFixed(3)));
  const sigBox = htmlEl("div", "pi-group", "", row2);

  const row3 = htmlEl("div", "pi-row", "", host);
  const kS = slider(row3, "step k", 0, KMAX, 1, st.k, v => { stop(); st.k = v; draw(); });
  const kg = htmlEl("div", "pi-group", "", row3);
  button(kg, "k = 0", () => { stop(); setK(0); });
  button(kg, "step +1", () => { stop(); setK(Math.min(KMAX, st.k + 1)); });
  button(kg, "+10", () => { stop(); setK(Math.min(KMAX, st.k + 10)); });
  const playB = button(kg, "play", () => { if (timer) stop(); else play(); });
  const sg = htmlEl("div", "pi-group", "", row3);
  htmlEl("span", "", "start:", sg);
  const startBox = htmlEl("span", "pi-group", "", sg);
  button(kg, "reset", () => {
    stop(); st.inv = start.inv; loadPreset(start.preset, false);
    st.sigma = Number.isFinite(start.sigma) ? start.sigma : start.preset.sigma; sS.set(st.sigma);
    st.k = start.k; kS.set(st.k); update();
  });

  const row4 = htmlEl("div", "pi-row pi-ent", "", host);
  const entS = ["a₁₁", "a₁₂", "a₂₁", "a₂₂"].map((lab, i) =>
    slider(row4, lab, -1, 5, 0.05, 0, v => {
      st.A[i >> 1][i & 1] = v; pBtns.forEach(b => b.classList.remove("active")); update();
    }, v => minus(v.toFixed(2))));

  const plots = htmlEl("div", "pi-plots", "", host);
  const H = 345, WA = 440, WB = 660, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Iterates and eigenvectors" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Error against step, log scale" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro = htmlEl("div", "qb-readout pi-read", "", host);

  let timer = 0;
  let sigButtons: { vals: number[]; btns: HTMLButtonElement[] } = { vals: [], btns: [] };
  function stop(): void { if (timer) { clearInterval(timer); timer = 0; } playB.textContent = "play"; }
  function play(): void {
    if (st.k >= KMAX) setK(0);
    playB.textContent = "pause";
    timer = window.setInterval(() => { if (st.k >= KMAX || !host.isConnected) { stop(); return; } setK(st.k + 1); }, 450);
  }
  function setK(k: number): void { st.k = k; kS.set(k); draw(); }

  function loadPreset(p: Preset, keepK: boolean): void {
    st.A = p.A.map(r => [...r]); st.x0 = [...p.x0]; st.names = p.names; st.cols = p.cols;
    st.sigma = p.sigma; sS.set(st.sigma);
    sigBox.replaceChildren();
    const sb = presets(sigBox, p.sigmas.map(v => ({ label: "σ = " + minus(String(v)), apply: () => { st.sigma = v; sS.set(v); st.inv = true; update(); } })));
    sigButtons = { vals: p.sigmas, btns: sb };
    pBtns.forEach((b, i) => b.classList.toggle("active", PRESETS[i] === p));
    if (!keepK) st.k = 0;
    buildStart();
    update();
  }
  function buildStart(): void {
    startBox.replaceChildren();
    const n = st.A.length;
    const items = st.names.map((nm, i) => ({ label: nm, apply: () => { st.x0 = st.A.map((_, j) => j === i ? 1 : 0); update(); } }));
    items.push({ label: "equal", apply: () => { st.x0 = new Array(n).fill(1); update(); } });
    presets(startBox, items);
  }

  // ---- computed state ----
  let an = analyse(st.A, st.inv ? st.sigma : null);
  let run = iterate(st.A, st.x0, st.inv ? st.sigma : null, KMAX);
  let ex: number[] = [], el: number[] = [];
  function update(): void {
    const sig = st.inv ? st.sigma : null;
    an = analyse(st.A, sig);
    run = iterate(st.A, st.x0, sig, KMAX);
    const t = an.target;
    ex = run.xs.map(x => t ? vecErr(x, t.v) : NaN);
    el = run.lam.map(l => t ? Math.abs(l - t.lam) : NaN);
    mBtns.forEach((b, i) => b.classList.toggle("active", (i === 1) === st.inv));
    sigButtons.btns.forEach((b, i) => b.classList.toggle("active", st.inv && Math.abs(sigButtons.vals[i] - st.sigma) < 1e-9));
    sS.input.closest("label")!.classList.toggle("pi-off", !st.inv);
    const two = st.A.length === 2;
    row4.style.display = two && editable ? "flex" : "none";
    if (two) entS.forEach((s, i) => s.set(st.A[i >> 1][i & 1]));
    draw();
  }

  // ---- panel (a): the plane (n = 2) ----
  const PY0 = 36, PY1 = H - 40;                           // plane: a square of side PY1 - PY0
  const cx = 60 + (PY1 - PY0) / 2, R = 1.3;
  const wx = scale(-R, R, cx - (PY1 - PY0) / 2, cx + (PY1 - PY0) / 2), wy = scale(-R, R, PY1, PY0);
  const toWorld = (ev: PointerEvent): [number, number] => {
    const [px, py] = svgPoint(svgA, ev);
    return [wx.invert(px), wy.invert(py)];
  };
  function arrow(g: SVGElement, x: number, y: number, stroke: string, w: number, op = 1): void {
    const L = Math.hypot(x, y);
    if (L < 1e-9) return;
    const s = Math.min(1, R * 0.97 / Math.max(Math.abs(x), Math.abs(y)));
    const X = wx(x * s), Y = wy(y * s), ux = (X - wx(0)) / Math.hypot(X - wx(0), Y - wy(0)), uy = (Y - wy(0)) / Math.hypot(X - wx(0), Y - wy(0));
    svgEl("line", { x1: wx(0), y1: wy(0), x2: X - ux * 8, y2: Y - uy * 8, stroke, "stroke-width": w, "stroke-opacity": op }, g);
    const hx = X, hy = Y, bx = X - ux * 16, by = Y - uy * 16;
    svgEl("polygon", { points: `${hx},${hy} ${bx - uy * 7},${by + ux * 7} ${bx + uy * 7},${by - ux * 7}`, fill: stroke, "fill-opacity": op }, g);
  }

  function drawPlane(): void {
    const k = st.k;
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent = "(a) direction of xₖ";
    const clipId = `pi-clip-${uid}`;
    const cp = svgEl("clipPath", { id: clipId }, svgA);
    svgEl("rect", { x: wx(-R), y: wy(R), width: wx(R) - wx(-R), height: wy(-R) - wy(R) }, cp);
    const g = svgEl("g", { "clip-path": `url(#${clipId})` }, svgA);
    svgEl("rect", { x: wx(-R), y: wy(R), width: wx(R) - wx(-R), height: wy(-R) - wy(R), fill: "none", class: "grid" }, svgA);
    svgEl("line", { x1: wx(-R), x2: wx(R), y1: wy(0), y2: wy(0), class: "axis" }, g);
    svgEl("line", { x1: wx(0), x2: wx(0), y1: wy(-R), y2: wy(R), class: "axis" }, g);
    svgEl("circle", { cx: wx(0), cy: wy(0), r: wx(1) - wx(0), fill: "none", stroke: "#bbb", "stroke-width": 1.5 }, g);
    svgEl("text", { x: wx(0), y: PY1 + 24, "text-anchor": "middle" }, svgA).textContent = st.names[0] + " →";
    const n1 = svgEl("text", { x: wx(-R) - 12, y: wy(0), "text-anchor": "middle", transform: `rotate(-90 ${wx(-R) - 12} ${wy(0)})` }, svgA);
    n1.textContent = st.names[1] + " →";
    // eigen-directions of A (real eigenvalues only)
    const real = an.eig.filter(e => e.im === 0);
    if (real.length === 2 && Math.abs(real[0].re - real[1].re) > 1e-12) {
      real.forEach((e, i) => {
        const B = st.A.map((r, a) => r.map((v, b) => v - (a === b ? e.re : 0)));
        // null vector of the 2 x 2 matrix B: from the larger row
        const r = Math.hypot(B[0][0], B[0][1]) >= Math.hypot(B[1][0], B[1][1]) ? B[0] : B[1];
        let v = [-r[1], r[0]];
        if (Math.hypot(v[0], v[1]) === 0) v = [1, 0];
        const s = Math.hypot(v[0], v[1]); v = [v[0] / s, v[1] / s];
        if (v[0] < 0 || (v[0] === 0 && v[1] < 0)) v = [-v[0], -v[1]];
        const c = i === 0 ? COL.v1 : COL.v2;
        svgEl("line", { x1: wx(-2 * v[0]), y1: wy(-2 * v[1]), x2: wx(2 * v[0]), y2: wy(2 * v[1]), stroke: c, "stroke-width": 2.5, "stroke-dasharray": "9 6" }, g);
        const lx = -1.12 * v[0], ly = -1.12 * v[1];
        const tx = svgEl("text", { x: wx(lx) + 10, y: wy(ly) + (v[1] >= 0 ? 18 : -4), "font-weight": 600 }, g);
        tx.style.fill = c;
        tx.textContent = `v${i === 0 ? "₁" : "₂"}`;
      });
    }
    // trail of iterates (unit directions)
    const unit = (x: number[]) => { const s = norm2(x); return [x[0] / s, x[1] / s]; };
    for (let j = 0; j < k; j++) {
      const u = unit(run.xs[j]);
      if (!u.every(Number.isFinite)) continue;
      const op = 0.15 + 0.5 * (j + 1) / (k + 1);
      svgEl("line", { x1: wx(0), y1: wy(0), x2: wx(u[0]), y2: wy(u[1]), stroke: COL.x, "stroke-width": 1.5, "stroke-opacity": op }, g);
      svgEl("circle", { cx: wx(u[0]), cy: wy(u[1]), r: 4, fill: COL.x, "fill-opacity": op }, g);
    }
    const xk = run.xs[k];
    const uk = unit(xk);
    if (uk.every(Number.isFinite)) {
      // y = A x_k (for unit x_k), or the solve for inverse iteration
      if (!st.inv && k < KMAX) {
        const y = matvec(st.A, uk);
        arrow(g, y[0], y[1], COL.y, 2.5, 0.85);
        const s = Math.min(1, R * 0.9 / Math.max(Math.abs(y[0]), Math.abs(y[1]), 1e-9));
        const ty = svgEl("text", { x: wx(y[0] * s) + 8, y: wy(y[1] * s) - 8, "font-weight": 600 }, g);
        ty.style.fill = COL.y; ty.textContent = "A xₖ";
      }
      arrow(g, uk[0], uk[1], COL.ink, 4);
      const px = uk[1], py = uk[0];                       // screen normal to the arrow
      const tk = svgEl("text", { x: wx(uk[0] * 0.5) + 18 * px, y: wy(uk[1] * 0.5) + 18 * py + 6, "font-weight": 700, "text-anchor": "middle" }, g);
      tk.textContent = `x${k === 0 ? "₀" : "ₖ"}`;
    }
    // start handle
    const u0 = unit(st.x0);
    const h = svgEl("circle", { cx: wx(u0[0]), cy: wy(u0[1]), r: 13, class: "handle", stroke: COL.x0 }, svgA);
    draggable(h, toWorld, (x, y) => {
      if (Math.hypot(x, y) < 0.05) return;
      st.x0 = [x, y];
      pBtns.forEach(b => b.classList.remove("active"));
      update();
    });
    const hint = svgEl("text", { x: wx(R) + 12, y: PY0 + 18, class: "sm" }, svgA);
    for (const [i, line] of ["drag ○", "to move", "the start", "x₀"].entries())
      svgEl("tspan", { x: wx(R) + 12, dy: i ? 22 : 0 }, hint).textContent = line;
  }

  // ---- panel (a): bars (n > 2) ----
  function drawBars(): void {
    const k = st.k, t = an.target, n = st.A.length;
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent =
      t ? `(a) xₖ (bars), ${st.inv ? "target" : "v₁"} (dashed)` : "(a) xₖ (bars)";
    const xk = run.xs[k];
    const pos = t ? (t.v.every(c => c >= 0) || t.v.every(c => c <= 0)) : xk.every(c => c >= 0);
    let shown: number[], vs: number[] | null = null;
    if (pos && Math.abs(sum(xk)) > 1e-300) {
      shown = xk.map(c => c / sum(xk));
      if (t) vs = t.v.map(c => c / sum(t.v));
    } else {
      const s = norm2(xk); shown = xk.map(c => c / s);
      if (t) { const sg = shown.reduce((a, c, i) => a + c * t.v[i], 0) >= 0 ? 1 : -1; vs = t.v.map(c => sg * c); }
    }
    const lo = pos ? 0 : -1, hi = 1;
    const y = scale(lo, hi, PY1 - 30, PY0 + 34);
    const x0 = 70, x1 = WA - 14, bw = (x1 - x0) / n;
    for (let v = lo; v <= hi + 1e-9; v += pos ? 0.25 : 0.5) {
      svgEl("line", { x1: x0, x2: x1, y1: y(v), y2: y(v), class: Math.abs(v) < 1e-9 ? "axis" : "grid" }, svgA);
      svgEl("text", { x: x0 - 6, y: y(v) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = minus(v.toFixed(2).replace(/0$/, ""));
    }
    const yl = svgEl("text", { x: 16, y: (PY0 + PY1) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(PY0 + PY1) / 2})` }, svgA);
    yl.textContent = pos ? "proportion" : "component (unit length)";
    for (let i = 0; i < n; i++) {
      const xc = x0 + (i + 0.5) * bw, w = bw * 0.6;
      const c = st.cols?.[i] ?? COL.x;
      const v = Number.isFinite(shown[i]) ? Math.max(lo, Math.min(hi, shown[i])) : 0;
      svgEl("rect", { x: xc - w / 2, y: Math.min(y(0), y(v)), width: w, height: Math.abs(y(v) - y(0)), fill: c, "fill-opacity": 0.85 }, svgA);
      if (vs) {
        const tv = vs[i];
        svgEl("rect", { x: xc - w / 2 - 5, y: Math.min(y(0), y(tv)), width: w + 10, height: Math.max(1, Math.abs(y(tv) - y(0))), fill: "none", stroke: COL.ink, "stroke-width": 2.5, "stroke-dasharray": "7 5" }, svgA);
      }
      svgEl("text", { x: xc, y: PY1 - 6, "text-anchor": "middle" }, svgA).textContent = st.names[i];
      svgEl("text", { x: xc, y: y(Math.max(0, v, vs ? vs[i] : 0)) - 8, "text-anchor": "middle", class: "sm" }, svgA).textContent =
        Number.isFinite(shown[i]) ? minus(shown[i].toFixed(3)) : "–";
    }
  }

  // ---- panel (b): error against k ----
  const B0 = 78, B1 = WB - 16, BT = 40, BB = H - 52;
  const EY0 = -16, EY1 = 1;
  const bx = scale(0, KMAX, B0, B1), by = scale(EY0, EY1, BB, BT);
  const ly = (e: number) => by(Math.max(EY0, Math.min(EY1, Math.log10(Math.max(e, 1e-300)))));
  function drawErr(): void {
    const k = st.k, t = an.target;
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgB).textContent = "(b) error against step k (log scale)";
    for (let e = EY0; e <= EY1; e += 2) {
      svgEl("line", { x1: B0, x2: B1, y1: by(e), y2: by(e), class: "grid" }, svgB);
      svgEl("text", { x: B0 - 6, y: by(e) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = pow10(e);
    }
    svgEl("line", { x1: B0, x2: B1, y1: BB, y2: BB, class: "axis" }, svgB);
    const tick = KMAX <= 60 ? 10 : KMAX <= 120 ? 20 : 50;
    for (let j = 0; j <= KMAX; j += tick) {
      svgEl("line", { x1: bx(j), x2: bx(j), y1: BB, y2: BB + 5, class: "axis" }, svgB);
      svgEl("text", { x: bx(j), y: BB + 22, "text-anchor": "middle", class: "sm" }, svgB).textContent = String(j);
    }
    svgEl("text", { x: (B0 + B1) / 2, y: H - 8, "text-anchor": "middle" }, svgB).textContent = "step k";
    const yl = svgEl("text", { x: 16, y: (BT + BB) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(BT + BB) / 2})` }, svgB);
    yl.textContent = t ? "error" : "residual";
    const clipId = `pi-clipb-${uid}`;
    const cp = svgEl("clipPath", { id: clipId }, svgB);
    svgEl("rect", { x: B0, y: BT - 6, width: B1 - B0 + 6, height: BB - BT + 12 }, cp);
    const g = svgEl("g", { "clip-path": `url(#${clipId})` }, svgB);
    // predicted line through the vector error at an anchor step
    if (t && an.rate > 0) {
      let k0 = 0;
      for (let j = Math.min(8, KMAX); j >= 0; j--) if (ex[j] > 1e-12 && Number.isFinite(ex[j])) { k0 = j; break; }
      const C0 = ex[k0] / an.rate ** k0;
      if (C0 > 0 && Number.isFinite(C0)) {
        const pts: string[] = [];
        for (let j = 0; j <= KMAX; j++) {
          const v = Math.log10(C0) + j * Math.log10(an.rate);
          if (v < EY0 - 1) break;
          pts.push(`${bx(j)},${by(Math.min(EY1 + 1, v))}`);
        }
        svgEl("polyline", { points: pts.join(" "), fill: "none", stroke: COL.pred, "stroke-width": 2.5, "stroke-dasharray": "9 6" }, g);
      }
    }
    const series = (vals: number[], c: string, shape: "circle" | "tri", from: number) => {
      const pts: string[] = [];
      for (let j = from; j <= k; j++) if (Number.isFinite(vals[j])) pts.push(`${bx(j)},${ly(vals[j])}`);
      svgEl("polyline", { points: pts.join(" "), fill: "none", stroke: c, "stroke-width": 2 }, g);
      for (let j = from; j <= k; j++) {
        if (!Number.isFinite(vals[j])) continue;
        const X = bx(j), Y = ly(vals[j]), r = j === k ? 7 : 4;
        if (shape === "circle") svgEl("circle", { cx: X, cy: Y, r, fill: c, stroke: j === k ? COL.ink : "none", "stroke-width": 2 }, g);
        else svgEl("polygon", { points: `${X},${Y - r * 1.2} ${X - r},${Y + r * 0.8} ${X + r},${Y + r * 0.8}`, fill: c, stroke: j === k ? COL.ink : "none", "stroke-width": 2 }, g);
      }
    };
    if (t) { series(ex, COL.x, "circle", 0); series(el, COL.lam, "tri", 1); }
    else series(run.res, COL.res, "circle", 1);
    svgEl("line", { x1: bx(k), x2: bx(k), y1: BT, y2: BB, stroke: COL.ink, "stroke-width": 1, "stroke-dasharray": "2 3" }, svgB);
    // legend
    const lg = (txt: string, row: number, c: string, dash = "") => {
      const X = B0 + 16, Y = BB - 14 - (2 - row) * 24;
      svgEl("rect", { x: X - 6, y: Y - 18, width: 306, height: 24, fill: "#fff", "fill-opacity": 0.8 }, svgB);
      svgEl("line", { x1: X, x2: X + 26, y1: Y - 6, y2: Y - 6, stroke: c, "stroke-width": 3, "stroke-dasharray": dash }, svgB);
      const tx = svgEl("text", { x: X + 34, y: Y }, svgB); tx.style.fill = c === COL.pred ? "#444" : c; tx.textContent = txt;
    };
    if (t) {
      lg(st.inv ? "error in x" : "error in x (vs v₁)", 0, COL.x);
      lg("error in λ estimate", 1, COL.lam);
      lg(st.inv ? `predicted rate ${sci(an.rate, 3)}ᵏ` : `predicted |λ₂/λ₁|ᵏ = ${sci(an.rate, 3)}ᵏ`, 2, COL.pred, "9 6");
    } else lg("residual ‖Ax − λx‖", 0, COL.res);
  }

  function draw(): void {
    svgA.replaceChildren(); svgB.replaceChildren();
    if (st.A.length === 2) drawPlane(); else drawBars();
    drawErr();
    // readout
    const k = st.k, t = an.target;
    ro.replaceChildren();
    const item = (key: string, v: string, cls = "") => {
      const s = htmlEl("span", cls, "", ro);
      if (key) htmlEl("b", "", key + " ", s);
      s.append(v);
    };
    const eigs = an.eig;
    item("eigenvalues:", eigList(eigs));
    if (st.inv) {
      const d = eigs.map(e => cabs({ re: e.re - st.sigma, im: e.im })).sort((a, b) => a - b);
      item("rate |λ−σ|/|λ′−σ| =", `${sci(d[0], 3)}/${sci(d[1], 3)} = ${sci(an.rate, 3)}`);
    } else if (eigs.length > 1) item("|λ₂/λ₁| =", `${sci(cabs(eigs[1]), 4)}/${sci(cabs(eigs[0]), 4)} = ${sci(an.rate, 3)}`);
    if (an.rate > 0 && an.rate < 1) item("steps per digit:", (1 / Math.log10(1 / an.rate)).toFixed(1));
    item("k =", String(k));
    item("λ estimate", k > 0 && Number.isFinite(run.lam[k]) ? minus(run.lam[k].toPrecision(8)) : "–");
    if (t) {
      item("target λ =", minus(t.lam.toPrecision(8)));
      item("error in x", sci(ex[k], 2));
      // measured rate over the last up to 5 steps with error above the rounding floor
      let j1 = k; while (j1 > 0 && !(ex[j1] > 1e-13)) j1--;
      const j0 = Math.max(0, j1 - 5);
      if (an.rate < 1 && j1 - j0 >= 2 && ex[j0] > 0) item("measured rate", sci((ex[j1] / ex[j0]) ** (1 / (j1 - j0)), 3));
    }
    item("residual", k > 0 ? sci(run.res[k], 2) : "–");
    item("normalised by", run.useSum ? "the sum (proportions)" : "the 2-norm");
    if (an.why) item("", an.why, "pi-warn");
    if (st.inv && run.tiny) item("", "A − σI is singular to working precision: the solve is huge, but its direction is right", "pi-warn");
  }

  loadPreset(start.preset, true);
  if (Number.isFinite(start.sigma)) { st.sigma = start.sigma; sS.set(st.sigma); update(); }
}

autoMount(NAME, mount);

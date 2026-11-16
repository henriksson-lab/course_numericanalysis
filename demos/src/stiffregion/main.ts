// Stability regions and stiffness, for the stiffness unit.
//  (a) the complex z = h·lambda plane with the stability regions |R(z)| <= 1 of Euler, RK4,
//      backward Euler and the trapezoidal rule (the selected method shaded), and the points
//      h·lambda for the two eigenvalues of the current problem. The step slider (or dragging
//      a point) moves both points along their rays from 0: the fast eigenvalue leaves an
//      explicit region long before the slow one gets anywhere near the edge.
//  (b) the solution of the 2x2 linear system y' = A y with the same h, exact (grey) and by
//      the selected method; computed in closed form, mode by mode, as R(h lambda)^n.
//  Readout table: per method, the stability limit on h, |R| for both eigenvalues, the
//  error at t = T with the current h, and the number of equal steps needed for error <= tol.
// Problems: eigenvalues -1 and -k (slider for the stiffness ratio k), a linear system with
// the enzyme model's eigenvalues at t = 0 (-0.0833 and -1201 per s, 30 s), and the 24-hour
// oscillator (+-i omega, omega = 2 pi / 24 per hour, 10 days).
// Auto-mounts on <div data-demo="stiffregion">. Optional attributes:
//   data-problem="stiff|enzyme|clock"  data-method="euler|rk4|beuler|trap"
//   data-logh="-3"  (log10 of the starting h)  data-ratio="1000"  data-tol="1e-3"
//   data-logt="0|1" (log time axis in (b))
import { svgEl, htmlEl, slider, presets, button, checkbox, scale, svgPoint, draggable, autoMount } from "../lib/ui";
import {
  C, cabs, Method, METHODS, LABEL, RFORM, R, rayLimit, regionPolygon, Problem, eigs, exact, numeric,
  errAtT, stepsNeeded, NMAX,
} from "./math";

const NAME = "stiffregion";
// Okabe-Ito, as in make-stiffness-figs.R (fig03: Euler blue, RK4 green, backward Euler red)
const COL: Record<Method, string> = { euler: "#0072B2", rk4: "#009E73", beuler: "#D55E00", trap: "#CC79A7" };
const GREY = "#8c8c8c", INK = "#222", WARN = "#b00020";
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .sr-row { display: flex; flex-wrap: wrap; gap: .35em 1.2em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .sr-group { display: flex; gap: .4em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .sr-row .qb-slider { flex: 1 1 300px; max-width: 440px; grid-template-columns: auto 1fr 5.4em; }
[data-demo="${NAME}"] .sr-row .qb-slider.off { opacity: .4; }
[data-demo="${NAME}"] .sr-plots { display: flex; flex-wrap: wrap; gap: 8px; width: 100%; }
[data-demo="${NAME}"] .sr-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] table.sr-tab { border-collapse: collapse; font-size: .9em; width: 100%; }
[data-demo="${NAME}"] .sr-tab th, [data-demo="${NAME}"] .sr-tab td { padding: .12em .55em; text-align: right; white-space: nowrap; }
[data-demo="${NAME}"] .sr-tab th { font-weight: 600; border-bottom: 1px solid #bbb; }
[data-demo="${NAME}"] .sr-tab th:first-child, [data-demo="${NAME}"] .sr-tab td:first-child { text-align: left; }
[data-demo="${NAME}"] .sr-tab tr.sel td { background: var(--qb-panel); font-weight: 600; }
[data-demo="${NAME}"] .sr-tab tbody tr { cursor: pointer; }
[data-demo="${NAME}"] .sr-tab td.bad { color: ${WARN}; }
[data-demo="${NAME}"] .sr-tab td .sw { display: inline-block; width: .8em; height: .8em; border-radius: 2px; margin-right: .35em; vertical-align: -.05em; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 18px; gap: .25em; }
.reveal [data-demo="${NAME}"] .sr-tab { font-size: 17px; }
.reveal [data-demo="${NAME}"] .sr-row .qb-slider { flex: 1 1 260px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- problems ----
type PKey = "stiff" | "enzyme" | "clock";
const OMEGA = 2 * Math.PI / 24;
const BASE: Record<PKey, Problem & { lh: number; label: string }> = {
  stiff: { key: "stiff", complex: false, ls: -1, lf: -1000, T: 4, unit: "", lh: -3, label: "λ = −1, −k" },
  enzyme: { key: "enzyme", complex: false, ls: -0.08327, lf: -1200.9, T: 30, unit: " s", lh: -3, label: "enzyme eigenvalues" },
  clock: { key: "clock", complex: true, ls: 0, lf: OMEGA, T: 240, unit: " h", lh: 0, label: "24-h clock, ±iω" },
};

// ---- formatting ----
const sup: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const pow10 = (k: number): string => k === 0 ? "1" : "10" + String(k).split("").map(c => sup[c]).join("");
function sci(v: number, d = 2): string {
  if (!Number.isFinite(v)) return v > 0 ? "∞" : "–";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e-3 && a < 1e4) return String(+v.toPrecision(d));
  const [m, e] = v.toExponential(d - 1).split("e");
  return `${m}e${Number(e)}`;
}
const fmtInt = (n: number): string => Math.round(n).toString().replace(/\B(?=(\d{3})+(?!\d))/g, " ");
const fmtN = (n: number): string => Number.isFinite(n) ? fmtInt(n) : "> " + pow10(Math.log10(NMAX));
function fmtErr(e: number): string {
  if (!Number.isFinite(e) || e > 1e100) return "overflow";
  return sci(e, 2);
}

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (s: string | undefined, def: number) => (s !== undefined && s !== "" && Number.isFinite(Number(s)) ? Number(s) : def);
  const pk0: PKey = ds.problem === "enzyme" || ds.problem === "clock" ? ds.problem : "stiff";
  const m0: Method = (METHODS as string[]).includes(ds.method ?? "") ? ds.method as Method : "euler";
  const init = {
    pk: pk0, m: m0,
    lh: num(ds.logh, BASE[pk0].lh),
    lr: Math.log10(Math.min(1e5, Math.max(1, num(ds.ratio, 1000)))),
    ltol: Math.log10(Math.min(0.1, Math.max(1e-8, num(ds.tol, 1e-3)))),
    logt: ds.logt === "1",
  };
  const st = { ...init };
  const prob = (): Problem => {
    const b = BASE[st.pk];
    if (st.pk === "stiff") return { ...b, lf: -(10 ** st.lr) };
    return b;
  };
  const lhRange = (): [number, number] => { const l = Math.log10(BASE[st.pk].T); return [Math.round((l - 6) * 100) / 100, Math.round(l * 100) / 100]; };
  const N = (): number => Math.max(1, Math.round(BASE[st.pk].T / 10 ** st.lh));
  const hEff = (): number => BASE[st.pk].T / N();

  // ---- controls ----
  const r1 = htmlEl("div", "sr-row", "", host);
  const g1 = htmlEl("div", "sr-group", "", r1);
  htmlEl("span", "", "problem:", g1);
  const pKeys: PKey[] = ["stiff", "enzyme", "clock"];
  const pBtns = presets(g1, pKeys.map(k => ({ label: BASE[k].label, apply: () => setProblem(k) })));
  const g2 = htmlEl("div", "sr-group", "", r1);
  htmlEl("span", "", "method:", g2);
  const mBtns = presets(g2, METHODS.map(m => ({ label: LABEL[m], apply: () => { st.m = m; draw(); } })));

  const r2 = htmlEl("div", "sr-row", "", host);
  const hS = slider(r2, "step h", lhRange()[0], lhRange()[1], 0.01, st.lh, v => { st.lh = v; draw(); },
    v => sci(10 ** v, 2));
  const limBtn = button(r2, "h → stability limit", () => {
    const lim = limitH(prob(), st.m);
    if (Number.isFinite(lim) && lim > 0) setLh(Math.log10(0.99 * lim));
  });
  const kS = slider(r2, "ratio k", 0, 5, 0.05, st.lr, v => { st.lr = v; draw(); }, v => fmtInt(10 ** v));
  const kLab = kS.input.parentElement as HTMLElement;

  const r3 = htmlEl("div", "sr-row", "", host);
  const tS = slider(r3, "tolerance", -8, -1, 0.5, st.ltol, v => { st.ltol = v; draw(); }, v => sci(10 ** v, 2));
  const ltBox = checkbox(r3, "log time axis", st.logt, v => { st.logt = v; draw(); });
  const g3 = htmlEl("div", "sr-group", "", r3);
  button(g3, "reset", () => {
    Object.assign(st, init);
    setRanges(); hS.set(st.lh); kS.set(st.lr); tS.set(st.ltol); ltBox.checked = st.logt; draw();
  });

  function setRanges(): void {
    const [a, b] = lhRange();
    hS.input.min = String(a); hS.input.max = String(b);
  }
  function setLh(v: number): void {
    const [a, b] = lhRange();
    st.lh = Math.min(b, Math.max(a, v)); hS.set(st.lh); draw();
  }
  function setProblem(k: PKey): void {
    st.pk = k; setRanges();
    st.lh = BASE[k].lh; hS.set(st.lh);
    st.logt = false; ltBox.checked = false;
    draw();
  }

  // ---- plots ----
  const plots = htmlEl("div", "sr-plots", "", host);
  const HH = 360, WA = 440, WB = 700, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${HH}`, role: "img", "aria-label": "Stability regions in the complex h lambda plane" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 4px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${HH}`, role: "img", "aria-label": "Solution of y' = Ay, exact and numerical" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 4px)`;
  const tab = htmlEl("table", "sr-tab", "", host);

  // (a) the z plane, equal aspect: 46 px per unit
  const U = 46, aL = 50, aR = 430, aT = 36, aB = 316;
  const X0 = -6.2, X1 = X0 + (aR - aL) / U, Y1 = (aB - aT) / 2 / U, Y0 = -Y1;
  const zx = scale(X0, X1, aL, aR), zy = scale(Y0, Y1, aB, aT);
  const rk4Poly = regionPolygon("rk4", -1.2);
  const uid = Math.random().toString(36).slice(2, 8);

  // (b) solution panel
  const bL = 70, bR = WB - 16, bT = 40, bB = HH - 56;

  // handles stay alive across redraws (they are being dragged)
  const handleLayer = svgEl("g", {}, svgA);
  const handles: SVGCircleElement[] = [0, 1].map(i => {
    const c = svgEl("circle", { r: 12, fill: "#fff", stroke: INK, "stroke-width": 3 }, handleLayer);
    draggable(c, ev => { const [px, py] = svgPoint(svgA, ev); return [zx.invert(px), zy.invert(py)]; },
      (x, y) => dragTo(i, x, y));
    return c;
  });
  function dragTo(i: number, x: number, y: number): void {
    const lam = eigs(prob())[i];
    const h = (x * lam[0] + y * lam[1]) / (lam[0] ** 2 + lam[1] ** 2);
    if (h > 0) setLh(Math.log10(h));
  }

  function limitH(p: Problem, m: Method): number {
    let lim = Infinity;
    for (const l of eigs(p)) {
      const a = cabs(l);
      lim = Math.min(lim, rayLimit(m, [l[0] / a, l[1] / a]) / a);
    }
    return lim;
  }

  const txt = (svg: Element, x: number, y: number, s: string, attrs: Record<string, string | number> = {}, fill?: string) => {
    const t = svgEl("text", { x, y, ...attrs }, svg);
    t.textContent = s;
    if (fill) t.style.fill = fill;
    return t;
  };

  function draw(): void {
    const p = prob(), m = st.m, n = N(), h = hEff();
    pBtns.forEach((b, i) => b.classList.toggle("active", pKeys[i] === st.pk));
    mBtns.forEach((b, i) => b.classList.toggle("active", METHODS[i] === m));
    kS.input.disabled = st.pk !== "stiff";
    kLab.classList.toggle("off", st.pk !== "stiff");
    if (st.pk === "enzyme") kS.set(Math.log10(p.lf / p.ls));
    const limSel = limitH(p, m);
    limBtn.disabled = !(Number.isFinite(limSel) && limSel > 0);
    ltBox.disabled = p.complex;
    const logt = st.logt && !p.complex;
    const L = eigs(p);
    const Z: C[] = L.map(l => [h * l[0], h * l[1]]);

    // ================= (a) =================
    for (const c of Array.from(svgA.childNodes)) if (c !== handleLayer) svgA.removeChild(c);
    const back = svgEl("g", {}, svgA);
    svgA.insertBefore(back, handleLayer);
    const defs = svgEl("defs", {}, back);
    const cp = svgEl("clipPath", { id: "sr-clipA-" + uid }, defs);
    svgEl("rect", { x: aL, y: aT, width: aR - aL, height: aB - aT }, cp);
    txt(back, 6, 24, "(a) stability regions, z = hλ", { class: "ttl" });
    const g = svgEl("g", { "clip-path": `url(#sr-clipA-${uid})` }, back);
    svgEl("rect", { x: aL, y: aT, width: aR - aL, height: aB - aT, fill: "#fff" }, g);
    // selected region shaded
    const fillOp = 0.28;
    const poly = (pts: C[]) => pts.map(q => `${zx(q[0]).toFixed(1)},${zy(q[1]).toFixed(1)}`).join(" ");
    if (m === "euler") svgEl("circle", { cx: zx(-1), cy: zy(0), r: U, fill: COL.euler, "fill-opacity": fillOp }, g);
    if (m === "rk4") svgEl("polygon", { points: poly(rk4Poly), fill: COL.rk4, "fill-opacity": fillOp }, g);
    if (m === "beuler") {
      const mask = svgEl("mask", { id: "sr-mask-" + uid }, defs);
      svgEl("rect", { x: aL, y: aT, width: aR - aL, height: aB - aT, fill: "#fff" }, mask);
      svgEl("circle", { cx: zx(1), cy: zy(0), r: U, fill: "#000" }, mask);
      svgEl("rect", { x: aL, y: aT, width: aR - aL, height: aB - aT, fill: COL.beuler, "fill-opacity": fillOp, mask: `url(#sr-mask-${uid})` }, g);
    }
    if (m === "trap") svgEl("rect", { x: aL, y: aT, width: zx(0) - aL, height: aB - aT, fill: COL.trap, "fill-opacity": fillOp }, g);
    // grid and axes
    for (let v = Math.ceil(X0); v <= X1; v++) {
      svgEl("line", { x1: zx(v), x2: zx(v), y1: aT, y2: aB, class: v === 0 ? "axis" : "grid" }, g);
      txt(back, zx(v), aB + 20, v === 0 ? "0" : String(v).replace("-", "−"), { "text-anchor": "middle", class: "sm" });
    }
    for (let v = Math.ceil(Y0); v <= Y1; v++) {
      svgEl("line", { x1: aL, x2: aR, y1: zy(v), y2: zy(v), class: v === 0 ? "axis" : "grid" }, g);
      txt(back, aL - 6, zy(v) + 6, v === 0 ? "0" : String(v).replace("-", "−") + "i", { "text-anchor": "end", class: "sm" });
    }
    txt(back, (aL + aR) / 2, HH - 10, "Re z", { "text-anchor": "middle" });
    // outlines of all four regions
    const w = (mm: Method) => mm === m ? 3.5 : 1.8;
    svgEl("circle", { cx: zx(-1), cy: zy(0), r: U, fill: "none", stroke: COL.euler, "stroke-width": w("euler") }, g);
    svgEl("polygon", { points: poly(rk4Poly), fill: "none", stroke: COL.rk4, "stroke-width": w("rk4") }, g);
    svgEl("circle", { cx: zx(1), cy: zy(0), r: U, fill: "none", stroke: COL.beuler, "stroke-width": w("beuler"), "stroke-dasharray": "7 4" }, g);
    svgEl("line", { x1: zx(0), x2: zx(0), y1: aT, y2: aB, stroke: COL.trap, "stroke-width": w("trap"), "stroke-dasharray": "7 4" }, g);
    svgEl("rect", { x: aL, y: aT, width: aR - aL, height: aB - aT, fill: "none", stroke: "#999" }, back);
    // legend
    const legend: [Method, string][] = [
      ["euler", "Euler: disc"], ["rk4", "RK4"],
      ["beuler", "backward Euler: outside circle"], ["trap", "trapezoid: Re z ≤ 0"],
    ];
    legend.forEach(([mm, s], i) => {
      const t = txt(back, aL + 8, aT + 22 + i * 22, s, { "font-weight": mm === m ? 700 : 400 }, COL[mm]);
      t.style.paintOrder = "stroke"; t.style.stroke = "#fff"; t.style.strokeWidth = "4px";
    });
    // eigenvalue points h*lambda
    const names = p.complex ? ["hλ₁ = +iωh", "hλ₂ = −iωh"] : ["hλ_slow", "hλ_fast"];
    Z.forEach((z, i) => {
      let px = zx(z[0]), py = zy(z[1]), off = false;
      const xmin = aL + 14, ymin = aT + 14, ymax = aB - 14;
      if (px < xmin) { px = xmin; off = true; }
      if (py < ymin) { py = ymin; off = true; }
      if (py > ymax) { py = ymax; off = true; }
      handles[i].setAttribute("cx", px.toFixed(1));
      handles[i].setAttribute("cy", py.toFixed(1));
      const stab = cabs(R(m, z)) <= 1 + 1e-9;
      handles[i].setAttribute("stroke", stab ? INK : WARN);
      handles[i].setAttribute("fill", off ? "#ffe9a8" : "#fff");
      // ray from 0 (dotted), so the radial motion is visible
      svgEl("line", { x1: zx(0), y1: zy(0), x2: px, y2: py, stroke: INK, "stroke-width": 1.2, "stroke-dasharray": "2 4" }, g);
      const zs = p.complex ? (i ? "−" : "") + sci(Math.abs(z[1]), 3) + "i" : sci(z[0], 3).replace("-", "−");
      // short tag at the handle, value in the lower-left corner
      const tag = p.complex ? (i ? "−iωh" : "+iωh") : (i ? "fast" : "slow");
      let lx = px, ly: number, anchor = "middle";
      if (p.complex) {
        lx = px + 18; ly = py + 6; anchor = "start";
        if (lx > aR - 60) { lx = px - 18; anchor = "end"; }
      } else {
        ly = i === 0 ? py - 20 : py + 34;
        if (px < aL + 30) { lx = px - 10; anchor = "start"; }
      }
      const t = txt(back, lx, ly, (off && !p.complex ? "◀ " : "") + tag, { "text-anchor": anchor, "font-weight": 700 }, stab ? INK : WARN);
      t.style.paintOrder = "stroke"; t.style.stroke = "#fff"; t.style.strokeWidth = "4px";
      const v = txt(back, aL + 8, aB - 34 + i * 24, `${names[i]} = ${zs}${off ? " (off scale)" : ""}`, { "font-weight": 600 }, stab ? INK : WARN);
      v.style.paintOrder = "stroke"; v.style.stroke = "#fff"; v.style.strokeWidth = "4px";
    });

    // ================= (b) =================
    svgB.replaceChildren();
    const unit = p.unit.trim();
    txt(svgB, 6, 24, `(b) y′ = Ay: exact and ${LABEL[m]}, h = ${sci(h, 2)}${p.unit}, N = ${fmtInt(n)} steps`, { class: "ttl" });
    const [ya, yb] = p.complex ? [-1.6, 1.6] : [-0.5, 2.3];
    const by = scale(ya, yb, bB, bT);
    const t0 = logt ? p.T * 1e-5 : 0;
    const bxLin = scale(0, p.T, bL, bR), bxLog = scale(Math.log10(t0 || 1), Math.log10(p.T), bL, bR);
    const bx = (t: number) => logt ? bxLog(Math.log10(Math.max(t, t0))) : bxLin(t);
    const defsB = svgEl("defs", {}, svgB);
    const cpB = svgEl("clipPath", { id: "sr-clipB-" + uid }, defsB);
    svgEl("rect", { x: bL, y: bT, width: bR - bL, height: bB - bT }, cpB);
    for (let v = Math.ceil(ya * 2) / 2; v <= yb + 1e-9; v += 0.5) {
      svgEl("line", { x1: bL, x2: bR, y1: by(v), y2: by(v), class: Math.abs(v) < 1e-9 ? "axis" : "grid" }, svgB);
      txt(svgB, bL - 6, by(v) + 6, String(v).replace("-", "−"), { "text-anchor": "end", class: "sm" });
    }
    if (logt) {
      for (let k = Math.ceil(Math.log10(t0)); k <= Math.log10(p.T); k++) {
        svgEl("line", { x1: bx(10 ** k), x2: bx(10 ** k), y1: bT, y2: bB, class: "grid" }, svgB);
        txt(svgB, bx(10 ** k), bB + 22, pow10(k), { "text-anchor": "middle", class: "sm" });
      }
    } else {
      const step = p.T / (p.T === 240 ? 10 : p.T === 30 ? 6 : 4);
      for (let v = 0; v <= p.T + 1e-9; v += step) {
        svgEl("line", { x1: bx(v), x2: bx(v), y1: bT, y2: bB, class: v === 0 ? "axis" : "grid" }, svgB);
        txt(svgB, bx(v), bB + 22, String(+v.toFixed(2)), { "text-anchor": "middle", class: "sm" });
      }
    }
    txt(svgB, (bL + bR) / 2, HH - 10, `time t${unit ? ` (${unit})` : ""}${logt ? ", log axis" : ""}`, { "text-anchor": "middle" });
    const pg = svgEl("g", { "clip-path": `url(#sr-clipB-${uid})` }, svgB);
    const cl = (v: number) => Math.max(ya - 1, Math.min(yb + 1, Number.isFinite(v) ? v : (v > 0 ? yb + 1 : ya - 1)));
    // exact: dense in t, plus extra points through the fast transient
    let ts: number[] = [];
    if (logt) for (let i = 0; i <= 500; i++) ts.push(t0 * (p.T / t0) ** (i / 500));
    else {
      for (let i = 0; i <= 500; i++) ts.push(p.T * i / 500);
      if (!p.complex) for (let j = 1; j <= 40; j++) ts.push(j * 0.25 / Math.abs(p.lf));
    }
    ts = ts.filter(t => t <= p.T).sort((a, b) => a - b);
    const ex = ts.map(t => exact(p, t));
    for (const c of [0, 1]) {
      svgEl("polyline", { points: ts.map((t, i) => `${bx(t).toFixed(1)},${by(cl(ex[i][c])).toFixed(1)}`).join(" "),
        fill: "none", stroke: GREY, "stroke-width": 6, "stroke-opacity": 0.75 }, pg);
    }
    // numerical: closed form at selected step indices
    let idx: number[] = [];
    if (n <= 3000) idx = Array.from({ length: n + 1 }, (_, i) => i);
    else if (logt) {
      const s = new Set<number>([0]);
      for (let i = 0; i <= 2000; i++) s.add(Math.round((n) ** (i / 2000)));
      for (let i = 0; i <= 1000; i++) s.add(Math.round(n * i / 1000));
      idx = Array.from(s).sort((a, b) => a - b);
    } else idx = Array.from({ length: 2001 }, (_, i) => Math.round(n * i / 2000));
    if (logt) idx = idx.filter(i => i > 0);
    const ys = idx.map(i => numeric(p, m, h, i));
    for (const c of [0, 1]) {
      svgEl("polyline", { points: idx.map((k, i) => `${bx(k * h).toFixed(1)},${by(cl(ys[i][c])).toFixed(1)}`).join(" "),
        fill: "none", stroke: COL[m], "stroke-width": 2.5, "stroke-dasharray": c ? "8 5" : "" }, pg);
      if (n <= 60) idx.forEach((k, i) => svgEl("circle", { cx: bx(k * h), cy: by(cl(ys[i][c])), r: 4.5, fill: COL[m] }, pg));
    }
    svgEl("rect", { x: bL, y: bT, width: bR - bL, height: bB - bT, fill: "none", stroke: "#999" }, svgB);
    // legend
    const lgx = bR - 250, lgy = bT + 24;
    const lgBox = svgEl("rect", { x: lgx - 10, y: lgy - 18, width: 254, height: 72, fill: "#fff", "fill-opacity": 0.85 }, svgB);
    void lgBox;
    svgEl("line", { x1: lgx, x2: lgx + 26, y1: lgy - 6, y2: lgy - 6, stroke: GREY, "stroke-width": 6, "stroke-opacity": 0.75 }, svgB);
    txt(svgB, lgx + 34, lgy, "exact", {});
    svgEl("line", { x1: lgx, x2: lgx + 26, y1: lgy + 16, y2: lgy + 16, stroke: COL[m], "stroke-width": 2.5 }, svgB);
    txt(svgB, lgx + 34, lgy + 22, p.complex ? `${LABEL[m]}: y₁ (cos)` : `${LABEL[m]}: y₁`, {}, COL[m]);
    svgEl("line", { x1: lgx, x2: lgx + 26, y1: lgy + 38, y2: lgy + 38, stroke: COL[m], "stroke-width": 2.5, "stroke-dasharray": "8 5" }, svgB);
    txt(svgB, lgx + 34, lgy + 44, p.complex ? "y₂ (sin)" : "y₂", {}, COL[m]);
    // blow-up / damping note
    const gmax = Math.max(...Z.map(z => cabs(R(m, z))));
    let note = "", noteCol = INK;
    if (gmax > 1 + 1e-9) {
      const nb = Math.ceil(6 * Math.LN10 / Math.log(gmax));
      noteCol = WARN;
      note = nb <= n ? `|R| = ${sci(gmax, 3)} > 1: |y| passes 10⁶ after ${fmtInt(nb)} steps (t = ${sci(nb * h, 2)})`
        : `|R| = ${sci(gmax, 4)} > 1: grows ×${sci(gmax ** n, 2)} over the N steps`;
    } else if (p.complex) {
      note = `amplitude after N steps: ${sci(gmax ** n, 3)} (exact: 1)`;
    } else {
      const rf = R(m, Z[1])[0];
      note = `fast mode: R = ${sci(rf, 3)} per step${rf < 0 ? " (sign flips every step)" : ""}`;
    }
    const nt = txt(svgB, bL + 10, bB - 12, note, { "font-weight": 600 }, noteCol);
    nt.style.paintOrder = "stroke"; nt.style.stroke = "#fff"; nt.style.strokeWidth = "5px";

    // ================= table =================
    const tol = 10 ** st.ltol;
    tab.replaceChildren();
    const thead = htmlEl("thead", "", "", tab), hr = htmlEl("tr", "", "", thead);
    const heads = ["method", "stable for h ≤", p.complex ? "|R(hλ)|" : "|R(hλ_slow)|", p.complex ? "amplitude at T" : "|R(hλ_fast)|",
      `error at T = ${p.T}${p.unit}`, `steps for error ≤ ${sci(tol, 1)}`];
    heads.forEach(s => htmlEl("th", "", s, hr));
    const tb = htmlEl("tbody", "", "", tab);
    for (const mm of METHODS) {
      const tr = htmlEl("tr", mm === m ? "sel" : "", "", tb);
      tr.addEventListener("click", () => { st.m = mm; draw(); });
      const td0 = htmlEl("td", "", "", tr);
      const sw = htmlEl("span", "sw", "", td0); sw.style.background = COL[mm];
      td0.append(LABEL[mm]);
      const lim = limitH(p, mm);
      htmlEl("td", "", lim === 0 ? "none" : Number.isFinite(lim) ? sci(lim, 3) + p.unit : "any h", tr);
      const r0 = cabs(R(mm, Z[0])), r1v = cabs(R(mm, Z[1]));
      const second = p.complex ? r0 ** n : r1v;
      htmlEl("td", r0 > 1 + 1e-9 ? "bad" : "", sci(r0, 4), tr);
      htmlEl("td", (p.complex ? second > 1 + 1e-9 : r1v > 1 + 1e-9) ? "bad" : "", sci(second, 4), tr);
      const e = errAtT(p, mm, n);
      htmlEl("td", e > tol ? "bad" : "", fmtErr(e), tr);
      htmlEl("td", "", fmtN(stepsNeeded(p, mm, tol)), tr);
    }
    const cap = htmlEl("caption", "qb-note", "", tab);
    cap.style.captionSide = "bottom"; cap.style.textAlign = "left";
    cap.textContent = `Same h for every row (N = ${fmtInt(n)} steps). ${RFORM[m]}. Red: |R| > 1 (unstable) or error above the tolerance. Click a row or drag a point in (a) to change method or h.`;
  }

  setRanges();
  draw();
}

autoMount(NAME, mount);

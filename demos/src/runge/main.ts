// Runge explorer for the polynomial-interpolation unit. The lecture's logistic OD curve,
// N(t) = K / (1 + (K/N0 - 1) e^{-rt}), K = 1.2, N0 = 0.05, r = 0.8, read at n + 1 times on
// [0, T] (equally spaced or Chebyshev), and the degree-n interpolating polynomial,
// evaluated with the barycentric formula.
//  (a) the true curve, the readings (draggable: in time, or in value) and the interpolant;
//      optionally reading noise (sd 0.02 OD, fixed seed) and a least-squares fit of low degree;
//  (b) the largest error on [0, T] against n = 1 ... 30 on a log axis, for both node sets,
//      with a marker at the current design.
// With T = 12, n = 12, equally spaced and noise on, the readings are exactly the lecture's
// y12 (R's set.seed(1)); other designs draw their noise from a fixed seed.
// Auto-mounts on <div data-demo="runge">. Optional attributes:
//   data-n="12"  data-nodes="equi|cheb"  data-tmax="48|12|8"  data-noise="0|1"
//   data-fit="0|1"  data-deg="5"
import { svgEl, htmlEl, slider, presets, checkbox, button, scale, svgPoint, autoMount } from "../lib/ui";
import { NodeSet, od, nodes, noiseFor, baryWeights, baryEval, lebesgue, lsFit, grid, maxErr, SD } from "./math";

const NAME = "runge";
// Okabe-Ito, as in make-polynomial-interpolation-figs.R: equispaced vermillion, Chebyshev blue.
const COL = { equi: "#D55E00", cheb: "#0072B2", fit: "#009E73", cust: "#CC79A7", truth: "#333333", ref: "#737373" };
const NMAX = 30;
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .rg-row { display: flex; flex-wrap: wrap; gap: .35em 1.2em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .rg-group { display: flex; gap: .35em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .rg-row .qb-slider { min-width: 300px; flex: 1 1 300px; max-width: 440px; grid-template-columns: auto 1fr 2.6em; }
[data-demo="${NAME}"] .rg-row .qb-slider.rg-small { min-width: 220px; flex: 0 1 260px; }
[data-demo="${NAME}"] .rg-plots { display: flex; flex-wrap: wrap; gap: 10px; width: 100%; }
[data-demo="${NAME}"] .rg-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .rg-read { display: flex; flex-wrap: wrap; gap: .15em 1.4em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .rg-read b { font-weight: 600; }
[data-demo="${NAME}"] .rg-warn { color: #b00020; font-weight: 600; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; gap: .3em; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 18px; padding: .3em .6em; }
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
function sci(v: number, d = 3): string {
  if (!Number.isFinite(v)) return String(v);
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e-3 && a < 1e5) return v.toPrecision(d);
  const [m, e] = v.toExponential(d - 1).split("e");
  return `${m}e${Number(e)}`;
}
const fmtT = (t: number): string => t.toFixed(t < 10 ? 2 : 1);

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (s: string | undefined, def: number) => (s !== undefined && s !== "" && Number.isFinite(Number(s)) ? Number(s) : def);
  const tOpt = [8, 12, 48];
  const init = {
    n: Math.round(Math.min(NMAX, Math.max(1, num(ds.n, 12)))),
    set: (ds.nodes === "cheb" ? "cheb" : "equi") as NodeSet,
    T: tOpt.includes(num(ds.tmax, 48)) ? num(ds.tmax, 48) : 48,
    noise: ds.noise === "1",
    fit: ds.fit === "1" || ds.noise === "1",
    deg: Math.round(Math.min(8, Math.max(1, num(ds.deg, 5)))),
  };
  const st = { ...init, drag: "time" as "time" | "value" };
  const uid = ++mountCount;

  // current readings: times, reading errors (noise) and hand-made offsets (value drag)
  let tt: number[] = [], nz: number[] = [], off: number[] = [];
  let custom = false;
  const rebuild = () => {
    tt = nodes(st.set, st.n, st.T);
    nz = st.noise ? noiseFor(st.set, st.n, st.T) : tt.map(() => 0);
    off = tt.map(() => 0);
    custom = false;
  };
  const yy = () => tt.map((t, i) => od(t) + nz[i] + off[i]);

  // error against n for both node sets (for the current T and noise setting)
  let curves: Record<NodeSet, number[]> = { equi: [], cheb: [] };
  let curvesKey = "";
  const ensureCurves = () => {
    const key = `${st.T}|${st.noise}`;
    if (key === curvesKey) return;
    curvesKey = key;
    const g = grid(st.T);
    for (const set of ["equi", "cheb"] as NodeSet[]) {
      curves[set] = [NaN];
      for (let n = 1; n <= NMAX; n++) {
        const t = nodes(set, n, st.T), e = st.noise ? noiseFor(set, n, st.T) : t.map(() => 0);
        const y = t.map((v, i) => od(v) + e[i]), w = baryWeights(t, st.T);
        curves[set].push(maxErr(g, x => baryEval(t, y, w, x)).e);
      }
    }
  };

  // ---- controls ----
  const row1 = htmlEl("div", "rg-row", "", host);
  const nS = slider(row1, "degree n", 1, NMAX, 1, st.n, v => { st.n = v; rebuild(); draw(); });
  const g1 = htmlEl("div", "rg-group", "", row1);
  htmlEl("span", "", "times:", g1);
  const setBtns = presets(g1, [
    { label: "equally spaced", apply: () => { st.set = "equi"; rebuild(); draw(); } },
    { label: "Chebyshev", apply: () => { st.set = "cheb"; rebuild(); draw(); } },
  ]);
  const g2 = htmlEl("div", "rg-group", "", row1);
  htmlEl("span", "", "run:", g2);
  const tBtns = presets(g2, tOpt.map(T => ({ label: `0–${T} h`, apply: () => { st.T = T; rebuild(); draw(); } })));

  const row2 = htmlEl("div", "rg-row", "", host);
  const noiseCb = checkbox(row2, `noise (sd ${SD} OD)`, st.noise, v => {
    st.noise = v;
    if (v && !st.fit) { st.fit = true; fitCb.checked = true; }
    rebuild(); draw();
  });
  const fitCb = checkbox(row2, "least-squares fit", st.fit, v => { st.fit = v; draw(); });
  const dS = slider(row2, "fit degree", 1, 8, 1, st.deg, v => { st.deg = v; draw(); });
  dS.input.parentElement!.classList.add("rg-small");
  const g3 = htmlEl("div", "rg-group", "", row2);
  htmlEl("span", "", "drag moves:", g3);
  const dragBtns = presets(g3, [
    { label: "time", apply: () => { st.drag = "time"; } },
    { label: "value", apply: () => { st.drag = "value"; } },
  ]);
  dragBtns[0].classList.add("active");

  const row3 = htmlEl("div", "rg-row", "", host);
  const g4 = htmlEl("div", "rg-group", "", row3);
  htmlEl("span", "", "presets:", g4);
  const apply = (p: Partial<typeof init>) => {
    Object.assign(st, init, p);
    nS.set(st.n); dS.set(st.deg); noiseCb.checked = st.noise; fitCb.checked = st.fit;
    rebuild(); draw();
  };
  const preBtns = presets(g4, [
    { label: "48 h, every 4 h", apply: () => apply({ n: 12, set: "equi", T: 48, noise: false, fit: false }) },
    { label: "48 h, Chebyshev", apply: () => apply({ n: 12, set: "cheb", T: 48, noise: false, fit: false }) },
    { label: "0–8 h, every 2 h", apply: () => apply({ n: 4, set: "equi", T: 8, noise: false, fit: false }) },
    { label: "noisy hourly, 12 h", apply: () => apply({ n: 12, set: "equi", T: 12, noise: true, fit: true, deg: 5 }) },
  ]);
  button(row3, "reset", () => { preBtns.forEach(b => b.classList.remove("active")); apply({}); });

  // ---- plots ----
  const plots = htmlEl("div", "rg-plots", "", host);
  const H = 340, WA = 660, WB = 470, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "OD curve, readings and interpolating polynomial" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  svgA.style.touchAction = "none";
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Largest error against degree, log scale" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro = htmlEl("div", "qb-readout rg-read", "", host);

  const top = 34, bot = H - 50;
  const aL = 62, aR = WA - 12;
  const Y0 = -0.5, Y1 = 1.7;
  const ay = scale(Y0, Y1, bot, top);
  const ax = () => scale(0, st.T, aL, aR);
  const bL = 70, bR = WB - 12;
  const EY0 = -9, EY1 = 1;
  const bx = scale(0.5, NMAX + 0.5, bL, bR), by = scale(EY0, EY1, bot, top);
  const ly = (e: number) => by(Math.max(EY0 - 0.3, Math.min(EY1 + 0.2, Math.log10(Math.max(e, 1e-30)))));

  // ---- dragging a reading (pointer handled on the svg, so redraws do not break it) ----
  let dragI = -1;
  const toWorld = (ev: PointerEvent): [number, number] => {
    const [px, py] = svgPoint(svgA, ev);
    return [ax().invert(px), ay.invert(py)];
  };
  svgA.addEventListener("pointerdown", ev => {
    const [px, py] = svgPoint(svgA, ev);
    const y = yy(), x = ax();
    let best = -1, bd = 22 * 22;
    tt.forEach((t, i) => {
      const d = (x(t) - px) ** 2 + (ay(Math.max(Y0, Math.min(Y1, y[i]))) - py) ** 2;
      if (d < bd) { bd = d; best = i; }
    });
    if (best < 0) return;
    ev.preventDefault();
    dragI = best;
    svgA.setPointerCapture(ev.pointerId);
    const move = (e: PointerEvent) => {
      const [wt, wy] = toWorld(e);
      if (st.drag === "time") {
        let t = Math.max(0, Math.min(st.T, wt));
        if (tt.some((s, j) => j !== dragI && Math.abs(s - t) < st.T * 1e-4)) return;   // keep times distinct
        tt[dragI] = t;
      } else {
        const t = tt[dragI];
        off[dragI] = Math.max(Y0, Math.min(Y1, wy)) - od(t) - nz[dragI];
      }
      custom = true;
      preBtns.forEach(b => b.classList.remove("active"));
      draw();
    };
    const up = () => {
      dragI = -1;
      svgA.removeEventListener("pointermove", move);
      svgA.removeEventListener("pointerup", up);
      svgA.removeEventListener("pointercancel", up);
      draw();
    };
    svgA.addEventListener("pointermove", move);
    svgA.addEventListener("pointerup", up);
    svgA.addEventListener("pointercancel", up);
  });

  function draw(): void {
    ensureCurves();
    nS.set(st.n);
    setBtns.forEach((b, i) => b.classList.toggle("active", !custom && (i === 0) === (st.set === "equi")));
    tBtns.forEach((b, i) => b.classList.toggle("active", tOpt[i] === st.T));
    const T = st.T, x = ax(), y = yy(), w = baryWeights(tt, T), g = grid(T);
    const p = (s: number) => baryEval(tt, y, w, s);
    const col = custom ? COL.cust : COL[st.set];
    const fit = st.fit ? lsFit(tt, y, st.deg, T) : null;

    // ---------- (a) ----------
    svgA.replaceChildren();
    const defs = svgEl("defs", {}, svgA);
    const cp = svgEl("clipPath", { id: `rg-clip-${uid}` }, defs);
    svgEl("rect", { x: aL, y: top, width: aR - aL, height: bot - top }, cp);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgA).textContent =
      `(a) ${st.n + 1} readings, degree-${st.n} interpolant`;
    for (let v = -0.4; v <= Y1 + 1e-9; v += 0.4) {
      const vv = Math.round(v * 10) / 10;
      svgEl("line", { x1: aL, x2: aR, y1: ay(vv), y2: ay(vv), class: vv === 0 ? "axis" : "grid" }, svgA);
      svgEl("text", { x: aL - 6, y: ay(vv) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = vv.toFixed(1);
    }
    const tstep = T === 48 ? 8 : T === 12 ? 2 : 1;
    for (let t = 0; t <= T + 1e-9; t += tstep) {
      svgEl("line", { x1: x(t), x2: x(t), y1: bot, y2: bot + 5, class: "axis" }, svgA);
      svgEl("text", { x: x(t), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgA).textContent = String(t);
    }
    svgEl("text", { x: (aL + aR) / 2, y: H - 6, "text-anchor": "middle" }, svgA).textContent = "time t (h)";
    svgEl("text", { x: 16, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(top + bot) / 2})` }, svgA).textContent = "OD600";
    const pg = svgEl("g", { "clip-path": `url(#rg-clip-${uid})` }, svgA);
    const NP = 900;
    const path = (f: (s: number) => number) => {
      let d = "";
      for (let i = 0; i <= NP; i++) {
        const s = T * i / NP, v = Math.max(Y0 - 1, Math.min(Y1 + 1, f(s)));
        d += `${i ? "L" : "M"}${x(s).toFixed(1)},${ay(v).toFixed(1)}`;
      }
      return d;
    };
    svgEl("path", { d: path(od), fill: "none", stroke: COL.truth, "stroke-width": 3 }, pg);
    svgEl("path", { d: path(p), fill: "none", stroke: col, "stroke-width": 3 }, pg);
    if (fit) svgEl("path", { d: path(fit.f), fill: "none", stroke: COL.fit, "stroke-width": 3, "stroke-dasharray": "9 6" }, pg);
    // where the interpolant leaves the plot
    let lo = Infinity, loAt = 0, hi = -Infinity, hiAt = 0;
    for (const s of g) { const v = p(s); if (v < lo) { lo = v; loAt = s; } if (v > hi) { hi = v; hiAt = s; } }
    const edge = (v: number, at: number, down: boolean) => {
      const tx = svgEl("text", { x: Math.max(aL + 70, Math.min(aR - 70, x(at))), y: down ? bot - 8 : top + 18, "text-anchor": "middle", "font-weight": 600 }, svgA);
      tx.style.fill = col; tx.textContent = `${down ? "↓" : "↑"} to ${sci(v, 3)}`;
    };
    if (lo < Y0) edge(lo, loAt, true);
    if (hi > Y1) edge(hi, hiAt, false);
    // readings: Chebyshev/equispaced ticks on the time axis, then the dots
    tt.forEach((t, i) => {
      svgEl("line", { x1: x(t), x2: x(t), y1: bot - 10, y2: bot, stroke: col, "stroke-width": 2.5 }, svgA);
      const yv = Math.max(Y0, Math.min(Y1, y[i]));
      svgEl("circle", { cx: x(t), cy: ay(yv), r: i === dragI ? 10 : 8, fill: col, stroke: "#fff", "stroke-width": 2, style: "cursor: grab" }, svgA);
    });
    // legend
    const leg = (txt: string, yy0: number, c: string, dash = "") => {
      svgEl("line", { x1: aR - 250, x2: aR - 222, y1: yy0 - 6, y2: yy0 - 6, stroke: c, "stroke-width": 3, "stroke-dasharray": dash }, svgA);
      const t = svgEl("text", { x: aR - 214, y: yy0 }, svgA);
      t.style.fill = c; t.textContent = txt;
    };
    const lgY = T === 48 ? ay(0.62) : ay(0.42);
    leg("true OD curve", lgY, COL.truth);
    leg(custom ? "interpolant (your times)" : st.set === "equi" ? "interpolant, equispaced" : "interpolant, Chebyshev", lgY + 24, col);
    if (fit) leg(`least squares, degree ${Math.min(st.deg, tt.length - 1)}`, lgY + 48, COL.fit, "9 6");

    // ---------- (b) ----------
    svgB.replaceChildren();
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgB).textContent = `(b) largest error on 0–${T} h`;
    for (let k = EY0; k <= EY1; k++) {
      if ((k - EY1) % 2) continue;
      svgEl("line", { x1: bL, x2: bR, y1: by(k), y2: by(k), class: "grid" }, svgB);
      svgEl("text", { x: bL - 6, y: by(k) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = pow10(k);
    }
    svgEl("line", { x1: bL, x2: bR, y1: bot, y2: bot, class: "axis" }, svgB);
    for (let n = 5; n <= NMAX; n += 5) {
      svgEl("line", { x1: bx(n), x2: bx(n), y1: bot, y2: bot + 5, class: "axis" }, svgB);
      svgEl("text", { x: bx(n), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgB).textContent = String(n);
    }
    svgEl("text", { x: (bL + bR) / 2, y: H - 6, "text-anchor": "middle" }, svgB).textContent = "degree n  (n + 1 readings)";
    svgEl("text", { x: 16, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(top + bot) / 2})` }, svgB).textContent = "max |p − N| (OD)";
    if (st.noise) {
      svgEl("line", { x1: bL, x2: bR, y1: by(Math.log10(SD)), y2: by(Math.log10(SD)), stroke: COL.ref, "stroke-width": 2, "stroke-dasharray": "3 5" }, svgB);
      svgEl("text", { x: bR - 4, y: by(Math.log10(SD)) + 22, "text-anchor": "end", class: "sm" }, svgB).textContent = "noise sd";
    }
    svgEl("line", { x1: bx(st.n), x2: bx(st.n), y1: top, y2: bot, stroke: "#999", "stroke-width": 1, "stroke-dasharray": "2 3" }, svgB);
    for (const set of ["equi", "cheb"] as NodeSet[]) {
      const c = COL[set], cv = curves[set];
      const pts = cv.map((e, n) => n ? `${bx(n)},${ly(e)}` : "").filter(Boolean).join(" ");
      svgEl("polyline", { points: pts, fill: "none", stroke: c, "stroke-width": 2.2 }, svgB);
      cv.forEach((e, n) => {
        if (!n) return;
        if (set === "equi") svgEl("circle", { cx: bx(n), cy: ly(e), r: 3.5, fill: c }, svgB);
        else svgEl("rect", { x: bx(n) - 3.5, y: ly(e) - 3.5, width: 7, height: 7, fill: c }, svgB);
      });
    }
    const lgB = (txt: string, yy0: number, c: string) => {
      svgEl("line", { x1: bL + 12, x2: bL + 38, y1: yy0 - 6, y2: yy0 - 6, stroke: c, "stroke-width": 3 }, svgB);
      const t = svgEl("text", { x: bL + 44, y: yy0 }, svgB);
      t.style.fill = c; t.textContent = txt;
    };
    lgB("equally spaced", by(-6.6), COL.equi);
    lgB("Chebyshev", by(-7.4), COL.cheb);
    const me = maxErr(g, p);
    svgEl("circle", { cx: bx(st.n), cy: ly(me.e), r: 9, fill: col, stroke: "#222", "stroke-width": 2.5, "fill-opacity": 0.85 }, svgB);

    // ---------- readout ----------
    let L = 0;
    for (let i = 0; i < g.length; i += 2) L = Math.max(L, lebesgue(tt, w, g[i]));
    ro.replaceChildren();
    const item = (k: string, v: string, cls = "") => {
      const s = htmlEl("span", cls, "", ro);
      htmlEl("b", "", k + " ", s); s.append(v);
    };
    item("interpolant: max error", `${sci(me.e)} OD at t = ${fmtT(me.at)} h`);
    item("lowest OD", `${sci(lo)} at ${fmtT(loAt)} h`, lo < 0 ? "rg-warn" : "");
    item("Lebesgue constant", sci(L, 3));
    if (fit) {
      const fe = maxErr(g, fit.f);
      const d = Math.min(st.deg, tt.length - 1);
      item(`fit (degree ${d}): max error`, d >= tt.length - 1 ? `${sci(fe.e)} (= interpolant)` : `${sci(fe.e)} OD, residual sd ${sci(fit.sd)}`);
    }
  }
  rebuild();
  draw();
}

autoMount(NAME, mount);

// Taylor explorer: a function f (eˣ, log(1+x) or √(1+x)) and its Taylor polynomial P_n
// around 0. (a) f and P_n on [-1, 2], with the error f(x) − P_n(x) at the chosen x drawn as
// a vertical bar. (b) |f − P_n| against |x| on log–log axes (computed in double precision),
// the band allowed by the remainder term, and the first omitted term (slope n + 1).
// Click or drag in either panel to move x. The readout gives f(x), P_n(x), the error, the
// remainder bound and the measured slope between x/10 and x. At small x the computed error
// hits the rounding floor (about 10⁻¹⁶) and the slope breaks: a preview of floating point.
// Auto-mounts on <div data-demo="taylor">. Optional: data-fn="exp|log|sqrt", data-n="2",
// data-x="0.1".
import { svgEl, htmlEl, slider, button, presets, scale, svgPoint, autoMount } from "../lib/ui";
import { FNS, FnKey, TaylorFn, taylor, formula, err, bound, nextTerm, slope, sup } from "./math";

const NAME = "taylor";
const COL = { f: "#000000", p: "#0072B2", e: "#D55E00", band: "#009E73", ref: "#555555", grey: "#8c8c8c" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { gap: 1em; }
[data-demo="${NAME}"].qb-demo .qb-plot { flex: 1 1 100%; max-width: 820px; }
[data-demo="${NAME}"].qb-demo .qb-panel { flex: 1 1 100%; max-width: 560px; gap: .5em; }
[data-demo="${NAME}"].qb-demo svg { width: 100%; height: auto; touch-action: none; cursor: crosshair; }
[data-demo="${NAME}"].qb-demo svg text { font-size: 19px; }
[data-demo="${NAME}"].qb-demo svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"].qb-demo svg .sm { font-size: 18px; fill: #444; }
[data-demo="${NAME}"].qb-demo .ty-row { display: flex; flex-wrap: wrap; gap: .35em; align-items: center; }
[data-demo="${NAME}"].qb-demo .ty-formula { font-size: .95em; }
[data-demo="${NAME}"].qb-demo .qb-readout { line-height: 1.45; }
.reveal [data-demo="${NAME}"].qb-demo { flex-wrap: nowrap; font-size: 20px; }
.reveal [data-demo="${NAME}"].qb-demo .qb-plot { flex: 0 0 770px; max-width: 770px; }
.reveal [data-demo="${NAME}"].qb-demo .qb-panel { flex: 0 0 360px; max-width: none; }
.reveal [data-demo="${NAME}"].qb-demo .qb-readout { font-size: 18px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- number formatting ----
/** x for display: 0.1, 1.5, 0.003, 1e-5, 3.2e-4. */
function fmtX(x: number): string {
  if (Math.abs(x) >= 0.001) return String(+x.toFixed(3));
  const [m, e] = x.toExponential(1).split("e");
  return `${+m}e${Number(e)}`;
}
/** Full precision, so the agreeing digits of f(x) and P_n(x) can be compared. */
function fmtFull(v: number): string {
  if (!Number.isFinite(v)) return String(v);
  if (v !== 0 && (Math.abs(v) < 1e-3 || Math.abs(v) >= 1e6)) return v.toExponential(15).replace("e+", "e");
  return v.toPrecision(16);
}
const fmtE = (v: number): string => (v === 0 ? "0" : v.toExponential(3).replace("e+", "e"));
/** Superscript power of ten for axis ticks: 10⁻¹⁶. */
const pow10 = (k: number): string => (k === 0 ? "1" : "10" + (k < 0 ? "⁻" : "") + sup(Math.abs(k)));

const XMAX = 2;          // right end of panel (a) and of the x slider
const UMIN = -6;         // panel (b): |x| from 10⁻⁶ ...
const UMAX = Math.log10(XMAX);
const VMIN = -20, VMAX = 2;   // ... and |error| from 10⁻²⁰ to 10²

function mount(host: HTMLElement): void {
  injectStyle();
  const d = host.dataset;
  const fnKey = (k: string | undefined): FnKey => (k === "log" || k === "sqrt" ? k : "exp");
  const num = (s: string | undefined, def: number) => (s !== undefined && s !== "" && Number.isFinite(Number(s)) ? Number(s) : def);
  const init = { fn: fnKey(d.fn), n: Math.round(Math.min(6, Math.max(0, num(d.n, 2)))), x: num(d.x, 0.1) };
  const st = { ...init };

  // ---- layout: plots left, controls and readout right ----
  const plot = htmlEl("div", "qb-plot", "", host);
  const panel = htmlEl("div", "qb-panel", "", host);
  const W = 790, H = 450, PW = 375;    // two panels side by side
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
    "aria-label": "Left: f and its Taylor polynomial. Right: the error against x on log–log axes." }, plot);
  const top = 40, bottom = H - 58;
  // panel (a)
  const aL = 46, aR = PW - 8;
  const ax = scale(-1, XMAX, aL, aR);
  // panel (b)
  const bL = PW + 84, bR = W - 8;
  const bx = scale(UMIN, UMAX, bL, bR);
  const by = scale(VMIN, VMAX, bottom, top);

  const defs = svgEl("defs", {}, svg);
  const clip = svgEl("clipPath", { id: "ty-clip-a" }, defs);
  svgEl("rect", { x: aL, y: top, width: aR - aL, height: bottom - top }, clip);
  const clipB = svgEl("clipPath", { id: "ty-clip-b" }, defs);
  svgEl("rect", { x: bL, y: top, width: bR - bL, height: bottom - top }, clipB);
  const layer = svgEl("g", {}, svg);

  // ---- controls ----
  const fnRow = htmlEl("div", "ty-row", "", panel);
  htmlEl("span", "", "f(x) =", fnRow);
  const keys: FnKey[] = ["exp", "log", "sqrt"];
  const fnBtns = presets(fnRow, keys.map(k => ({ label: FNS[k].label, apply: () => { st.fn = k; clampX(); sync(); draw(); } })));
  const nS = slider(panel, "degree n", 0, 6, 1, st.n, v => { st.n = v; draw(); });
  const xS = slider(panel, "point x", -0.9, XMAX, 0.01, st.x, v => {
    st.x = Math.round(v * 100) / 100; clampX(); xS.set(st.x); draw();
  }, fmtX);
  const stepRow = htmlEl("div", "ty-row", "", panel);
  const mul = (k: number) => { st.x = Number((st.x * k).toPrecision(12)); clampX(); xS.set(st.x); draw(); };
  button(stepRow, "x ÷ 10", () => mul(0.1));
  button(stepRow, "x ÷ 2", () => mul(0.5));
  button(stepRow, "x × 2", () => mul(2));
  button(stepRow, "x × 10", () => mul(10));
  button(stepRow, "x → −x", () => mul(-1));
  const formulaEl = htmlEl("div", "ty-formula", "", panel);
  const ro = htmlEl("div", "qb-readout", "", panel);
  htmlEl("div", "qb-note", "Click or drag in either plot to move x (in (b) on a log scale). Slope: measured on (b) between x/10 and x.", panel);
  const setAll = (fn: FnKey, n: number, x: number) => { st.fn = fn; st.n = n; st.x = x; clampX(); sync(); draw(); };
  const pre = presets(panel, [
    { label: "reset", apply: () => setAll(init.fn, init.n, init.x) },
    { label: "e^0.2, n = 3", apply: () => setAll("exp", 3, 0.2) },
    { label: "log(1.05)", apply: () => setAll("log", 2, 0.05) },
    { label: "√1.02", apply: () => setAll("sqrt", 2, 0.02) },
    { label: "log(1+x), x = 1.5", apply: () => setAll("log", 3, 1.5) },
    { label: "x = 10⁻⁵", apply: () => setAll(st.fn, st.n, 1e-5) },
  ]);
  pre.forEach(b => b.addEventListener("click", () => window.setTimeout(() => b.classList.remove("active"), 600)));

  function clampX(): void {
    const fn = FNS[st.fn];
    st.x = Math.min(XMAX, Math.max(Math.max(fn.xmin, -0.9), st.x));
    if (Math.abs(st.x) < 1e-6) st.x = st.x < 0 ? -1e-6 : 1e-6;
  }
  function sync(): void {
    nS.set(st.n);
    xS.set(st.x);
    fnBtns.forEach((b, i) => b.classList.toggle("active", keys[i] === st.fn));
  }

  // ---- pointer: click or drag in a panel moves x ----
  function pointTo(ev: PointerEvent): void {
    const [px, py] = svgPoint(svg, ev);
    if (py < top - 20 || py > bottom + 30) return;
    if (px < PW) {
      st.x = Math.round(ax.invert(px) * 1000) / 1000;
    } else {
      const u = Math.min(UMAX, Math.max(UMIN, bx.invert(px)));
      const s = st.x < 0 ? -1 : 1;
      st.x = s * Number((10 ** u).toPrecision(2));
    }
    clampX(); xS.set(st.x); draw();
  }
  svg.addEventListener("pointerdown", ev => {
    ev.preventDefault();
    svg.setPointerCapture(ev.pointerId);
    pointTo(ev);
    const move = (e: PointerEvent) => pointTo(e);
    const up = () => {
      svg.removeEventListener("pointermove", move);
      svg.removeEventListener("pointerup", up);
      svg.removeEventListener("pointercancel", up);
    };
    svg.addEventListener("pointermove", move);
    svg.addEventListener("pointerup", up);
    svg.addEventListener("pointercancel", up);
  });

  const txt = (x: number, y: number, s: string, attrs: Record<string, string | number> = {}, fill?: string) => {
    const t = svgEl("text", { x, y, ...attrs }, layer);
    t.textContent = s;
    if (fill) t.style.fill = fill;
    return t;
  };
  const path = (pts: [number, number][]): string => pts.map((p, i) => (i ? "L" : "M") + p[0].toFixed(1) + "," + p[1].toFixed(1)).join("");

  function draw(): void {
    const fn: TaylorFn = FNS[st.fn];
    const n = st.n, x = st.x;
    layer.replaceChildren();

    // ================= (a) f and P_n =================
    const [y0, y1] = fn.ylim;
    const ay = scale(y0, y1, bottom, top);
    txt(4, 24, `(a) ${fn.label} and P${subN(n)}`, { class: "ttl" });
    const yStep = y1 - y0 > 5 ? 2 : y1 - y0 > 3 ? 1 : 0.5;
    for (let v = Math.ceil(y0 / yStep) * yStep; v <= y1 + 1e-9; v += yStep) {
      svgEl("line", { x1: aL, x2: aR, y1: ay(v), y2: ay(v), class: Math.abs(v) < 1e-9 ? "axis" : "grid" }, layer);
      txt(aL - 6, ay(v) + 6, String(+v.toFixed(1)), { "text-anchor": "end", class: "sm" });
    }
    for (let v = -1; v <= XMAX; v += 0.5) {
      svgEl("line", { x1: ax(v), x2: ax(v), y1: top, y2: bottom, class: v === 0 ? "axis" : "grid" }, layer);
      if (Number.isInteger(v)) txt(ax(v), bottom + 22, String(v), { "text-anchor": "middle", class: "sm" });
    }
    txt((aL + aR) / 2, bottom + 46, "x", { "text-anchor": "middle" });
    const g = svgEl("g", { "clip-path": "url(#ty-clip-a)" }, layer);
    const N = 300;
    const xsA = Array.from({ length: N + 1 }, (_, i) => -1 + (XMAX + 1) * i / N);
    const curve = (h: (v: number) => number) => {
      const segs: [number, number][][] = [[]];
      for (const v of xsA) {
        const y = h(v);
        if (Number.isFinite(y)) segs[segs.length - 1].push([ax(v), ay(Math.max(y0 - 5, Math.min(y1 + 5, y)))]);
        else if (segs[segs.length - 1].length) segs.push([]);
      }
      return segs.filter(s => s.length > 1).map(path).join("");
    };
    svgEl("path", { d: curve(v => (v > -1 || fn.key === "exp" ? fn.f(v) : NaN)), fill: "none", stroke: COL.f, "stroke-width": 3.5 }, g);
    svgEl("path", { d: curve(v => taylor(fn, n, v)), fill: "none", stroke: COL.p, "stroke-width": 3, "stroke-dasharray": "9 5" }, g);
    const fx = fn.f(x), px = taylor(fn, n, x), e = fx - px;
    svgEl("line", { x1: ax(x), x2: ax(x), y1: top, y2: bottom, stroke: COL.grey, "stroke-width": 1.5, "stroke-dasharray": "3 4" }, g);
    const cy = (v: number) => ay(Math.max(y0 - 0.2, Math.min(y1 + 0.2, v)));
    svgEl("line", { x1: ax(x), x2: ax(x), y1: cy(fx), y2: cy(px), stroke: COL.e, "stroke-width": 6, "stroke-linecap": "round" }, g);
    svgEl("circle", { cx: ax(x), cy: cy(fx), r: 6.5, fill: COL.f, stroke: "#fff", "stroke-width": 1.5 }, g);
    svgEl("circle", { cx: ax(x), cy: cy(px), r: 6.5, fill: COL.p, stroke: "#fff", "stroke-width": 1.5 }, g);
    // handle on the x axis
    const hy = Math.min(bottom, Math.max(top, ay(0)));
    svgEl("rect", { x: ax(x) - 8, y: hy - 8, width: 16, height: 16, rx: 3, fill: "#fff", stroke: COL.e, "stroke-width": 3 }, layer);
    // legend
    const lgx = aL + 10;
    svgEl("line", { x1: lgx, x2: lgx + 30, y1: top + 14, y2: top + 14, stroke: COL.f, "stroke-width": 3.5 }, layer);
    txt(lgx + 36, top + 20, fn.label);
    svgEl("line", { x1: lgx, x2: lgx + 30, y1: top + 40, y2: top + 40, stroke: COL.p, "stroke-width": 3, "stroke-dasharray": "9 5" }, layer);
    txt(lgx + 36, top + 46, `P${subN(n)}`, {}, COL.p);
    svgEl("line", { x1: lgx + 15, x2: lgx + 15, y1: top + 58, y2: top + 76, stroke: COL.e, "stroke-width": 6, "stroke-linecap": "round" }, layer);
    txt(lgx + 36, top + 72, "error at x", {}, COL.e);

    // ================= (b) |error| against |x| =================
    const sgn = x < 0 ? -1 : 1;
    const uHi = sgn < 0 ? Math.log10(0.9) : UMAX;
    txt(PW + 8, 24, `(b) |error| against ${sgn < 0 ? "−x (x < 0)" : "x"}, log–log`, { class: "ttl" });
    for (let k = VMIN; k <= VMAX; k += 2) {
      const major = k % 4 === 0;
      svgEl("line", { x1: bL, x2: bR, y1: by(k), y2: by(k), class: k === 0 ? "axis" : "grid", "stroke-opacity": major ? 1 : 0.5 }, layer);
      if (major) txt(bL - 6, by(k) + 6, pow10(k), { "text-anchor": "end", class: "sm" });
    }
    for (let k = UMIN; k <= 0; k++) {
      svgEl("line", { x1: bx(k), x2: bx(k), y1: top, y2: bottom, class: "grid", "stroke-opacity": k % 2 === 0 ? 1 : 0.5 }, layer);
      if (k % 2 === 0) txt(bx(k), bottom + 22, pow10(k), { "text-anchor": "middle", class: "sm" });
    }
    svgEl("line", { x1: bL, x2: bL, y1: top, y2: bottom, class: "axis" }, layer);
    svgEl("line", { x1: bL, x2: bR, y1: bottom, y2: bottom, class: "axis" }, layer);
    txt((bL + bR) / 2, bottom + 46, sgn < 0 ? "|x|" : "x", { "text-anchor": "middle" });
    const yl = txt(PW + 14, (top + bottom) / 2, "", { "text-anchor": "middle", transform: `rotate(-90 ${PW + 14} ${(top + bottom) / 2})` });
    yl.textContent = `|f(x) − P${subN(n)}(x)|`;

    const gb = svgEl("g", { "clip-path": "url(#ty-clip-b)" }, layer);
    const lg = (v: number) => Math.log10(v);
    const clampV = (v: number) => Math.max(VMIN - 3, Math.min(VMAX + 3, v));
    const M = 360;
    const us = Array.from({ length: M + 1 }, (_, i) => UMIN + (uHi - UMIN) * i / M);
    // remainder band
    const lo: [number, number][] = [], hi: [number, number][] = [];
    for (const u of us) {
      const [b0, b1] = bound(fn, n, sgn * 10 ** u);
      lo.push([bx(u), by(clampV(lg(b0)))]);
      hi.push([bx(u), by(clampV(lg(b1)))]);
    }
    svgEl("path", { d: path([...hi, ...lo.reverse()]) + "Z", fill: COL.band, "fill-opacity": 0.22, stroke: "none" }, gb);
    // first omitted term: slope n + 1
    svgEl("path", { d: path(us.map(u => [bx(u), by(clampV(lg(nextTerm(fn, n, sgn * 10 ** u))))])),
      fill: "none", stroke: COL.ref, "stroke-width": 2, "stroke-dasharray": "7 5" }, gb);
    // rounding level
    const eps = Math.log10(2.220446049250313e-16);
    svgEl("line", { x1: bL, x2: bR, y1: by(eps), y2: by(eps), stroke: COL.grey, "stroke-width": 2, "stroke-dasharray": "2 4" }, gb);
    txt(bR - 4, by(eps) + 22, "rounding level 2.2·10⁻¹⁶", { "text-anchor": "end", class: "sm" });
    // computed error: break the line where it is exactly 0
    // exactly 0 (f and P_n round to the same double) cannot go on a log axis: tick marks
    // on the bottom edge instead
    const segs: [number, number][][] = [[]];
    let zeros = "";
    for (const u of us) {
      const a = Math.abs(err(fn, n, sgn * 10 ** u));
      if (a > 0 && Number.isFinite(a)) segs[segs.length - 1].push([bx(u), by(clampV(lg(a)))]);
      else {
        if (segs[segs.length - 1].length) segs.push([]);
        if (a === 0) zeros += `M${bx(u).toFixed(1)},${bottom}v-12`;
      }
    }
    for (const s of segs) if (s.length === 1) s.push([s[0][0] + 1.5, s[0][1]]);   // isolated points
    svgEl("path", { d: segs.filter(s => s.length > 1).map(path).join(""), fill: "none", stroke: COL.e, "stroke-width": 2.5, "stroke-linejoin": "round", "stroke-linecap": "round" }, gb);
    if (zeros) {
      svgEl("path", { d: zeros, stroke: COL.e, "stroke-width": 2, "stroke-opacity": 0.7 }, gb);
      txt(bL + 6, bottom - 18, "▮ error exactly 0", { class: "sm" }, COL.e);
    }
    // marker at the current x
    const ux = Math.log10(Math.abs(x));
    svgEl("line", { x1: bx(ux), x2: bx(ux), y1: top, y2: bottom, stroke: COL.grey, "stroke-width": 1.5, "stroke-dasharray": "3 4" }, gb);
    const ae = Math.abs(e);
    if (ae > 0) svgEl("circle", { cx: bx(ux), cy: by(clampV(lg(ae))), r: 7, fill: COL.e, stroke: "#fff", "stroke-width": 1.5 }, gb);
    svgEl("rect", { x: bx(ux) - 8, y: bottom - 8, width: 16, height: 16, rx: 3, fill: "#fff", stroke: COL.e, "stroke-width": 3 }, layer);
    // legend, top left (errors are small there)
    const lx = bL + 10;
    svgEl("line", { x1: lx, x2: lx + 30, y1: top + 14, y2: top + 14, stroke: COL.e, "stroke-width": 2.5 }, layer);
    txt(lx + 36, top + 20, "computed error", {}, COL.e);
    svgEl("line", { x1: lx, x2: lx + 30, y1: top + 40, y2: top + 40, stroke: COL.ref, "stroke-width": 2, "stroke-dasharray": "7 5" }, layer);
    txt(lx + 36, top + 46, `next term, slope ${n + 1}`);
    svgEl("rect", { x: lx, y: top + 56, width: 30, height: 16, fill: COL.band, "fill-opacity": 0.3 }, layer);
    txt(lx + 36, top + 70, "remainder bound");

    // ================= readout =================
    formulaEl.textContent = `P${subN(n)}(x) = ${formula(fn, n)}`;
    const [b0, b1] = bound(fn, n, x);
    const sl = slope(fn, n, x);
    const pn = `P${subN(n)}(x)`;
    ro.textContent =
      `x        = ${fmtX(x)}\n` +
      `f(x)     = ${fmtFull(fx)}\n` +
      `${pn.padEnd(8)} = ${fmtFull(px)}\n` +
      `error    = ${e === 0 ? "0 (same double)" : fmtE(e)}\n` +
      `bound    = ${fmtE(b0)} … ${fmtE(b1)}\n` +
      `slope    = ${Number.isFinite(sl) ? sl.toFixed(2) : "— (error 0 at x/10)"}   (n + 1 = ${n + 1})`;
  }

  sync();
  draw();
}

const SUB = "₀₁₂₃₄₅₆₇₈₉";
const subN = (n: number): string => String(n).split("").map(c => SUB[Number(c)]).join("");

autoMount(NAME, mount);

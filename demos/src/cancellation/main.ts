// Cancellation explorer for "Error, conditioning and stability".
// View 1 (variance): a fixed sample of n values offset + sd·z (z seeded, standardised to
//   mean 0, sd 1). Population variance in plain doubles by the one-pass formula
//   sum(x²)/n − mean², two-pass sum((x − mean)²)/n and Welford's update, compared with a
//   reference computed in double-double arithmetic from the same stored numbers. Plot:
//   correct digits = −log10(relative error) against log10(offset), with the rule
//   16 − 2·log10(offset/sd) and the conditioning limit 16 − log10(2·offset/sd).
// View 2 (kappa): condition number κ = |x f'(x)/f(x)| of x − 1, log x, √x and the dose
//   from a response D = K R/(1 − R), with digits lost = log10 κ.
// Auto-mounts on <div data-demo="cancellation">. Optional attributes:
//   data-view="variance|kappa"  data-offset="2.4e8"  data-sd="3"  data-n="1000"
//   data-seed="3"  data-func="sub|log|sqrt|dose"  data-x="1.001"  data-tabs="0" (hide tabs)
import { svgEl, htmlEl, slider, presets, scale, rng, autoMount, Slider } from "../lib/ui";
import { onePass, twoPass, welford, refVariance, digits } from "./math";

const NAME = "cancellation";
// Okabe–Ito, as in the unit's figures (make-error-conditioning-figs.R)
const COL = { one: "#D55E00", two: "#0072B2", wel: "#E69F00", cond: "#009E73", rule: "#666" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .5em; }
[data-demo="${NAME}"] .cc-body { display: flex; flex-wrap: wrap; gap: 1em; width: 100%; align-items: flex-start; }
[data-demo="${NAME}"] .cc-body .qb-panel { flex: 0 1 430px; }
[data-demo="${NAME}"] .cc-body .qb-plot { flex: 1 1 420px; max-width: 700px; }
[data-demo="${NAME}"] svg text { font-size: 21px; }
[data-demo="${NAME}"] svg .sm { font-size: 19px; fill: #444; }
[data-demo="${NAME}"] table.cc-tab { border-collapse: collapse; font-size: .9em; width: 100%; }
[data-demo="${NAME}"] table.cc-tab td, [data-demo="${NAME}"] table.cc-tab th { padding: .1em .35em; text-align: left; border: none; }
[data-demo="${NAME}"] table.cc-tab th { font-weight: 600; font-size: .85em; color: #555; }
[data-demo="${NAME}"] table.cc-tab td.num { font-family: ui-monospace, monospace; }
[data-demo="${NAME}"] table.cc-tab td.num b { color: #0b6b3a; }
[data-demo="${NAME}"] table.cc-tab td.num span { color: #999; }
[data-demo="${NAME}"] table.cc-tab td.dg { font-family: ui-monospace, monospace; text-align: right; font-weight: 600; }
[data-demo="${NAME}"] .cc-sw { display: inline-block; width: .8em; height: .8em; border-radius: 2px; margin-right: .3em; vertical-align: -.05em; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 20px; }
.reveal [data-demo="${NAME}"] .cc-body { flex-wrap: nowrap; }
.reveal [data-demo="${NAME}"] .cc-body .qb-panel { flex: 0 0 410px; }
.reveal [data-demo="${NAME}"] .cc-body .qb-plot { flex: 1 1 620px; max-width: none; }
.reveal [data-demo="${NAME}"] table.cc-tab { font-size: .9em; }
.reveal [data-demo="${NAME}"] .qb-note { font-size: .9em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

/** 2.4e8-style short format for offsets and sd. */
function fmtSci(v: number): string {
  if (v === 0) return "0";
  const e = Math.floor(Math.log10(Math.abs(v)) + 1e-9);
  if (e >= -2 && e < 3) return String(+v.toPrecision(2));
  const m = +(v / 10 ** e).toFixed(1);
  return m === 10 ? `1e${e + 1}` : `${m}e${e}`;
}
const num = (s: string | undefined, d: number) => (s !== undefined && s !== "" && Number.isFinite(Number(s)) ? Number(s) : d);

/** 16 significant digits, the first `d` significant digits (the correct ones) in bold. */
function markDigits(v: number, d: number): string {
  const a = v === 0 ? "0" : v.toPrecision(16);
  let k = 0, seen = 0, started = false;
  const want = Math.floor(d + 1e-9);
  while (k < a.length && seen < want) {
    const c = a[k];
    if (c >= "1" && c <= "9") started = true;
    if (started && c >= "0" && c <= "9") seen++;
    k++;
  }
  if (seen === 0) k = 0;
  return `<b>${a.slice(0, k)}</b><span>${a.slice(k)}</span>`;
}

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  if (ds.tabs !== "0") {
    const tabs = presets(host, [
      { label: "variance: one-pass vs two-pass", apply: () => show("variance") },
      { label: "condition number κ", apply: () => show("kappa") },
    ]);
    tabs[ds.view === "kappa" ? 1 : 0].classList.add("active");
  }
  const body = htmlEl("div", "cc-body", "", host);
  const show = (v: string) => { body.replaceChildren(); (v === "kappa" ? kappaView : varianceView)(body, ds); };
  show(ds.view === "kappa" ? "kappa" : "variance");
}

// ======================================================================================
// View 1: the variance
// ======================================================================================
const GRID = Array.from({ length: 241 }, (_, i) => i * 0.05);   // log10(offset) 0..12

function varianceView(body: HTMLElement, ds: DOMStringMap): void {
  const n = Math.max(10, Math.round(num(ds.n, 1000)));
  const r = rng(Math.round(num(ds.seed, 3)));
  const z = new Float64Array(n);
  for (let i = 0; i < n; i++) z[i] = r.norm();
  let m = 0; z.forEach(v => (m += v)); m /= n;
  let q = 0; z.forEach(v => (q += (v - m) ** 2)); const s = Math.sqrt(q / n);
  for (let i = 0; i < n; i++) z[i] = (z[i] - m) / s;              // mean 0, sd 1

  const start = { offset: num(ds.offset, 2.4e8), sd: num(ds.sd, 3) };
  const st = { ...start };
  const sample = (off: number, sd: number) => z.map(v => off + sd * v);
  const evalAt = (off: number, sd: number) => {
    const x = sample(off, sd), ref = refVariance(x);
    const one = onePass(x), two = twoPass(x), wel = welford(x);
    return { ref, one, two, wel, d1: digits(one, ref), d2: digits(two, ref), dw: digits(wel, ref) };
  };
  let cacheSd = NaN;
  let curves: { d1: number; d2: number; dw: number }[] = [];
  const getCurves = () => {
    if (st.sd !== cacheSd) { curves = GRID.map(g => evalAt(10 ** g, st.sd)); cacheSd = st.sd; }
    return curves;
  };

  // ---- controls ----
  const panel = htmlEl("div", "qb-panel", "", body);
  const oS = slider(panel, "offset", 0, 12, 0.05, Math.log10(st.offset), v => { st.offset = 10 ** v; clearPresets(); draw(); }, v => fmtSci(10 ** v));
  const sS = slider(panel, "sd", -2, 4, 0.05, Math.log10(st.sd), v => { st.sd = +(10 ** v).toPrecision(3); clearPresets(); draw(); }, v => fmtSci(10 ** v));
  const items = [
    { label: "positions: 2.4e8, sd 3", off: 2.4e8, sd: 3 },
    { label: "counts: 1e8, sd 1e4", off: 1e8, sd: 1e4 },
    { label: "small: 10, sd 3", off: 10, sd: 3 },
    { label: "reset", off: start.offset, sd: start.sd },
  ];
  const pb = presets(panel, items.map(it => ({ label: it.label, apply: () => {
    st.offset = it.off; st.sd = it.sd; oS.set(Math.log10(it.off)); sS.set(Math.log10(it.sd)); draw();
  } })));
  const clearPresets = () => pb.forEach(b => b.classList.remove("active"));
  const tab = htmlEl("table", "cc-tab", "", panel);
  const note = htmlEl("div", "qb-note", "", panel);

  // ---- plot ----
  const W = 680, H = 440, L = 62, R = W - 30, T = 34, B = H - 62;
  const plot = htmlEl("div", "qb-plot", "", body);
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Correct digits of the variance against the offset" }, plot);
  const sx = scale(0, 12, L, R), sy = scale(0, 16.5, B, T);

  function draw(): void {
    const cur = evalAt(st.offset, st.sd);
    const cv = getCurves();
    const lr = Math.log10(st.offset / st.sd);

    // table
    const row = (name: string, col: string, v: number, d: number | null) =>
      `<tr><td><span class="cc-sw" style="background:${col}"></span>${name}</td>` +
      `<td class="num">${d === null ? v.toPrecision(16) : markDigits(v, d)}</td>` +
      `<td class="dg" style="color:${col}">${d === null ? "" : d >= 15.95 ? "16" : d.toFixed(1)}</td></tr>`;
    tab.innerHTML = `<tr><th>method</th><th>variance</th><th>digits</th></tr>` +
      row("reference", "#222", cur.ref, null) +
      row("one-pass", COL.one, cur.one, cur.d1) +
      row("two-pass", COL.two, cur.two, cur.d2) +
      row("Welford", COL.wel, cur.wel, cur.dw);
    note.innerHTML =
      `n = ${n}, offset/sd = ${fmtSci(st.offset / st.sd)}. Rule for one-pass: 16 − 2·log₁₀(offset/sd) = ` +
      `<b>${Math.max(0, 16 - 2 * lr).toFixed(1)}</b>; problem allows 16 − log₁₀(2·offset/sd) = ` +
      `<b>${Math.min(16, 16 - Math.log10(2) - lr).toFixed(1)}</b>. ` +
      `Reference: variance of the stored numbers, in double-double (≈32 digits). Bold: the correct digits.`;

    // plot
    svg.replaceChildren();
    svgEl("text", { x: L - 50, y: 22, "font-weight": 600 }, svg).textContent =
      `correct digits of the variance (sd = ${fmtSci(st.sd)})`;
    for (let d = 0; d <= 16; d += 4) {
      svgEl("line", { x1: L, x2: R, y1: sy(d), y2: sy(d), class: d === 0 ? "axis" : "grid" }, svg);
      svgEl("text", { x: L - 8, y: sy(d) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = String(d);
    }
    for (let e = 0; e <= 12; e += 2) {
      svgEl("line", { x1: sx(e), x2: sx(e), y1: B, y2: B + 6, class: "axis" }, svg);
      const t = svgEl("text", { x: sx(e), y: B + 26, "text-anchor": "middle", class: "sm" }, svg);
      t.textContent = e === 0 ? "1" : `1e${e}`;
    }
    svgEl("text", { x: (L + R) / 2, y: H - 10, "text-anchor": "middle" }, svg).textContent = "offset (log axis)";
    // reference lines
    const line = (f: (g: number) => number, col: string, dash: string, w: number) => {
      const pts = GRID.map(g => `${sx(g)},${sy(Math.max(0, Math.min(16, f(g))))}`).join(" ");
      svgEl("polyline", { points: pts, fill: "none", stroke: col, "stroke-width": w, "stroke-dasharray": dash }, svg);
    };
    const lsd = Math.log10(st.sd);
    line(g => 16 - Math.log10(2) - (g - lsd), COL.cond, "2 5", 3.5);
    line(g => 16 - 2 * (g - lsd), COL.rule, "9 6", 2.5);
    const curve = (k: "d1" | "d2" | "dw", col: string, w: number, dash = "") =>
      svgEl("polyline", { points: GRID.map((g, i) => `${sx(g)},${sy(cv[i][k])}`).join(" "), fill: "none",
        stroke: col, "stroke-width": w, "stroke-dasharray": dash, "stroke-linejoin": "round" }, svg);
    curve("dw", COL.wel, 2.5);
    curve("d2", COL.two, 3);
    curve("d1", COL.one, 3);
    // current position
    const g0 = Math.log10(st.offset);
    if (g0 >= 0 && g0 <= 12) {
      svgEl("line", { x1: sx(g0), x2: sx(g0), y1: T, y2: B, stroke: "#222", "stroke-width": 1.5, "stroke-dasharray": "4 4" }, svg);
      const dot = (d: number, col: string, rr: number) =>
        svgEl("circle", { cx: sx(g0), cy: sy(d), r: rr, fill: col, stroke: "#fff", "stroke-width": 2 }, svg);
      dot(cur.dw, COL.wel, 7); dot(cur.d2, COL.two, 8); dot(cur.d1, COL.one, 8);
    }
    // legend in the lower left, which the curves leave empty
    const leg: [string, string, string, number][] = [
      ["two-pass", COL.two, "", 3], ["Welford", COL.wel, "", 2.5], ["one-pass", COL.one, "", 3],
      ["rule for one-pass", COL.rule, "9 6", 2.5], ["limit set by κ", COL.cond, "2 5", 3.5],
    ];
    const lx = L + 12, ly = B - 12 - 23 * (leg.length - 1);
    leg.forEach(([txt, col, dash, w], k) => {
      const y = ly + 23 * k;
      svgEl("line", { x1: lx, x2: lx + 30, y1: y - 6, y2: y - 6, stroke: col, "stroke-width": w, "stroke-dasharray": dash }, svg);
      svgEl("text", { x: lx + 38, y, class: "sm" }, svg).textContent = txt;
    });
  }
  draw();
}

// ======================================================================================
// View 2: condition numbers
// ======================================================================================
interface Fn { key: string; label: string; f: (x: number) => number; kappa: (x: number) => number;
  lo: number; hi: number; step: number; x0: number; xname: string; formula: string }
const K_DOSE = 10;
const FNS: Fn[] = [
  { key: "sub", label: "x − 1", f: x => x - 1, kappa: x => Math.abs(x / (x - 1)), lo: 0.5, hi: 2, step: 0.001, x0: 1.001, xname: "x", formula: "f = x − 1,  κ = |x / (x − 1)|" },
  { key: "log", label: "log x", f: Math.log, kappa: x => 1 / Math.abs(Math.log(x)), lo: 0.5, hi: 2, step: 0.001, x0: 1.001, xname: "x", formula: "f = log x,  κ = 1 / |log x|" },
  { key: "sqrt", label: "√x", f: Math.sqrt, kappa: () => 0.5, lo: 0.5, hi: 2, step: 0.001, x0: 1.001, xname: "x", formula: "f = √x,  κ = 1/2" },
  { key: "dose", label: "dose from response", f: R => K_DOSE * R / (1 - R), kappa: R => 1 / (1 - R), lo: 0.01, hi: 0.999, step: 0.001, x0: 0.99, xname: "R", formula: "D = K R/(1 − R),  κ = 1/(1 − R)" },
];

function kappaView(body: HTMLElement, ds: DOMStringMap): void {
  let fn = FNS.find(f => f.key === ds.func) ?? FNS[0];
  let x = num(ds.x, fn.x0);
  const panel = htmlEl("div", "qb-panel", "", body);
  const fb = presets(panel, FNS.map(f => ({ label: f.label, apply: () => { fn = f; x = f.x0; buildSlider(); draw(); } })));
  fb[FNS.indexOf(fn)].classList.add("active");
  const sHost = htmlEl("div", "", "", panel);
  let xS: Slider;
  const buildSlider = () => {
    sHost.replaceChildren();
    xS = slider(sHost, fn.xname, fn.lo, fn.hi, fn.step, x, v => { x = v; draw(); }, v => v.toFixed(3));
  };
  buildSlider();
  const ro = htmlEl("div", "qb-readout", "", panel);
  htmlEl("div", "qb-note", "Digits lost ≈ log₁₀ κ, whatever the algorithm. Doubles start with 16 digits, " +
    "a measurement with 1 % error with 2.", panel);

  const W = 680, H = 440, L = 62, R = W - 30, T = 34, B = H - 62;
  const plot = htmlEl("div", "qb-plot", "", body);
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "log10 of the condition number" }, plot);
  const YLO = -1, YHI = 4;
  const sy = scale(YLO, YHI, B, T);

  function draw(): void {
    const sx = scale(fn.lo, fn.hi, L, R);
    const k = fn.kappa(x), lk = Math.log10(k);
    const h = 1e-7, fx = fn.f(x);
    const kn = Math.abs(fn.f(x * (1 + h)) - fx) / Math.abs(fx) / h;
    const fmtK = (v: number) => (!Number.isFinite(v) ? "∞" : v >= 1e4 ? v.toExponential(2) : String(+v.toPrecision(5)));
    ro.textContent =
      `${fn.formula}\n` +
      `at ${fn.xname} = ${x.toFixed(3)}:  κ = ${fmtK(k)}   (kappa_num ${fmtK(kn)})\n` +
      `digits lost  log₁₀ κ = ${Number.isFinite(lk) ? lk.toFixed(1) : "∞"}\n` +
      `left of 16 (double):  ${Number.isFinite(lk) ? Math.max(0, 16 - Math.max(0, lk)).toFixed(1) : "0"}\n` +
      `left of 2 (1 % data): ${Number.isFinite(lk) ? Math.max(0, 2 - Math.max(0, lk)).toFixed(1) : "0"}`;

    svg.replaceChildren();
    svgEl("text", { x: L - 50, y: 22, "font-weight": 600 }, svg).textContent = `digits lost, log₁₀ κ, for f = ${fn.label}`;
    for (let d = YLO; d <= YHI; d++) {
      svgEl("line", { x1: L, x2: R, y1: sy(d), y2: sy(d), class: d === 0 ? "axis" : "grid" }, svg);
      svgEl("text", { x: L - 8, y: sy(d) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = String(d);
    }
    const ticks = fn.key === "dose" ? [0.2, 0.4, 0.6, 0.8, 1] : [0.5, 1, 1.5, 2];
    for (const t of ticks) {
      svgEl("line", { x1: sx(t), x2: sx(t), y1: B, y2: B + 6, class: "axis" }, svg);
      svgEl("text", { x: sx(t), y: B + 26, "text-anchor": "middle", class: "sm" }, svg).textContent = String(t);
    }
    svgEl("text", { x: (L + R) / 2, y: H - 10, "text-anchor": "middle" }, svg).textContent =
      fn.key === "dose" ? "response R (fraction of maximum)" : "x";
    // a 1 %-measurement has 2 digits: above this line none are left
    svgEl("line", { x1: L, x2: R, y1: sy(2), y2: sy(2), stroke: COL.one, "stroke-width": 2, "stroke-dasharray": "8 6" }, svg);
    const left = fn.key === "dose";
    const tl = svgEl("text", { x: left ? L + 8 : R - 8, y: sy(2) - 8, class: "sm", "text-anchor": left ? "start" : "end" }, svg);
    tl.style.fill = COL.one; tl.textContent = "1 % data: no digits left above this line";
    // dense grid, plus points crowding in on the singularity at 1 and the current x
    const N = 800;
    const xs = Array.from({ length: N + 1 }, (_, i) => fn.lo + (fn.hi - fn.lo) * i / N);
    for (let e = 1; e <= 6; e += 0.125) xs.push(1 - 10 ** -e, 1 + 10 ** -e);
    xs.push(x);
    xs.sort((a, b) => a - b);
    const seg: string[][] = [[]];
    for (const xv of xs) {
      if (xv < fn.lo || xv > fn.hi) continue;
      const v = Math.log10(fn.kappa(xv));
      if (!Number.isFinite(v) || v > YHI + 0.5) { if (seg[seg.length - 1].length) seg.push([]); continue; }
      seg[seg.length - 1].push(`${sx(xv)},${sy(Math.max(YLO, Math.min(YHI, v)))}`);
    }
    for (const s of seg) if (s.length > 1)
      svgEl("polyline", { points: s.join(" "), fill: "none", stroke: COL.two, "stroke-width": 3.5 }, svg);
    svgEl("line", { x1: sx(x), x2: sx(x), y1: T, y2: B, stroke: "#222", "stroke-width": 1.5, "stroke-dasharray": "4 4" }, svg);
    const yv = Number.isFinite(lk) ? Math.max(YLO, Math.min(YHI, lk)) : YHI;
    svgEl("circle", { cx: sx(x), cy: sy(yv), r: 9, fill: COL.two, stroke: "#fff", "stroke-width": 2 }, svg);
    if (!Number.isFinite(lk) || lk > YHI) {
      svgEl("text", { x: sx(x) + 12, y: sy(YHI) + 20, class: "sm" }, svg).textContent = "off the scale";
    }
  }
  draw();
}

autoMount(NAME, mount);

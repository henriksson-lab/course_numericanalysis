// Deconvolution explorer for the svd unit. The lecture's simulated signature S (40 genes x
// 5 cell types; CD8T = CD4T except genes 1 and 2, which differ by +-d, d = 2 % by default),
// true proportions p = (0.30, 0.15, 0.15, 0.10, 0.30), bulk b0 = S p, and a noisy sample
// b = b0 * (1 + s * noise): either the lecture's fixed pattern sin(37 g) or N(0, 1) noise.
// SVD by one-sided Jacobi on S itself (math.ts), recomputed when d changes.
//  (a) proportions: estimate for the sample (bars), the truth (dashed), and the estimates
//      from 100 fixed noise replicates (dots), so bias and spread are both visible;
//  (b) filter factors f_i = phi_i sigma_i against sigma (log axis): least squares (all 1),
//      truncated SVD (a step after sigma_k), ridge sigma^2 / (sigma^2 + lambda^2);
//  (c) expected error over the noise, split into bias ||E x - p|| and noise sqrt(E||x - E x||^2),
//      against lambda (ridge) or k (TSVD), computed exactly from the filtered pseudoinverse.
// Fixed seeds, so every run shows the same replicates.
// Auto-mounts on <div data-demo="deconv">. Optional attributes:
//   data-method="ls|tsvd|ridge"  data-k="4"  data-lambda="10"  data-noise="0.05"
//   data-diff="0.02"  data-sample="pattern|random"  data-reps="100"  data-seed="1"
import { svgEl, htmlEl, slider, presets, button, scale, rng, autoMount } from "../lib/ui";
import { signature, svd, phis, filterMatrix, matVec, norm, expectedError, Method, SVD, M, N, TYPES, P_TRUE, EPS } from "./math";

const NAME = "deconv";
// Okabe-Ito, as in make-svd-figs.R (fig07: least squares grey, TSVD orange, ridge sky blue; truth vermillion dashed)
const COL = { ls: "#999999", tsvd: "#E69F00", ridge: "#56B4E9", truth: "#D55E00", bias: "#CC79A7", noise: "#0072B2", ink: "#222", dot: "#555" };
const MCOL: Record<Method, string> = { ls: COL.ls, tsvd: COL.tsvd, ridge: COL.ridge };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .dc-row { display: flex; flex-wrap: wrap; gap: .35em 1.2em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .dc-group { display: flex; gap: .4em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .dc-row .qb-slider { min-width: 250px; flex: 1 1 250px; max-width: 400px; grid-template-columns: auto 1fr 3.8em; margin-right: .6em; }
[data-demo="${NAME}"] .dc-row .qb-slider.off { opacity: .45; }
[data-demo="${NAME}"] .dc-plots { display: flex; flex-wrap: wrap; gap: 6px; width: 100%; }
[data-demo="${NAME}"] .dc-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .dc-read { display: flex; flex-wrap: wrap; gap: .1em 1.3em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .dc-read b { font-weight: 600; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; gap: .25em; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; padding: .25em .5em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- formatting ----
const minus = (s: string): string => s.replace(/^-/, "−");
const f3 = (v: number): string => minus((Math.abs(v) < 5e-4 ? 0 : v).toFixed(3));
function g3(v: number): string {
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e5 || a < 1e-3) { const [m, e] = v.toExponential(1).split("e"); return `${minus(m)}e${Number(e)}`; }
  return minus(a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : a >= 1 ? v.toFixed(2) : v.toPrecision(2));
}
const lamFmt = (l: number): string => g3(10 ** l);
const pct = (v: number): string => `${(100 * v).toFixed(Math.abs(v * 100 - Math.round(v * 100)) < 1e-9 ? 0 : 1)}%`;
const SUB = ["₁", "₂", "₃", "₄", "₅"];

// ---- fixed noise ----
const PATTERN = Array.from({ length: M }, (_, g) => Math.sin(37 * (g + 1)));   // the lecture's "noise"
function normals(seed: number, count: number): number[][] {
  const r = rng(seed);
  return Array.from({ length: count }, () => Array.from({ length: M }, () => r.norm()));
}

const LAM0 = -1, LAM1 = 3;                     // log10 lambda range
const LGRID = Array.from({ length: 81 }, (_, i) => LAM0 + i * 0.05);

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (v: string | undefined, def: number, lo: number, hi: number): number => {
    const x = Number(v);
    return v !== undefined && Number.isFinite(x) ? Math.min(hi, Math.max(lo, x)) : def;
  };
  const start = {
    method: (["ls", "tsvd", "ridge"].includes(ds.method ?? "") ? ds.method : "ls") as Method,
    k: Math.round(num(ds.k, 4, 1, 5)),
    ll: Math.round(Math.log10(num(ds.lambda, 10, 0.1, 1000)) / 0.05) * 0.05,
    s: num(ds.noise, 0.05, 0, 0.2),
    d: num(ds.diff, 0.02, 0, 0.3),
    random: ds.sample === "random",
    seed: Math.round(num(ds.seed, 1, 1, 1e6)),
  };
  const nrep = Math.round(num(ds.reps, 100, 0, 500));
  const st = { ...start };
  const reps = normals(1000 + start.seed, nrep);
  let sampleZ = normals(start.seed, 1)[0];

  // ---- controls ----
  const row1 = htmlEl("div", "dc-row", "", host);
  const mg = htmlEl("div", "dc-group", "", row1);
  htmlEl("span", "", "method:", mg);
  const METHODS: Method[] = ["ls", "tsvd", "ridge"];
  const mBtns = presets(mg, [
    { label: "least squares", apply: () => { st.method = "ls"; draw(); } },
    { label: "truncated SVD", apply: () => { st.method = "tsvd"; draw(); } },
    { label: "ridge", apply: () => { st.method = "ridge"; draw(); } },
  ]);
  const kS = slider(row1, "keep k", 1, 5, 1, st.k, v => { st.k = v; st.method = "tsvd"; draw(); });
  const lS = slider(row1, "ridge λ", LAM0, LAM1, 0.05, st.ll, v => { st.ll = v; st.method = "ridge"; draw(); }, lamFmt);
  const row2 = htmlEl("div", "dc-row", "", host);
  const sS = slider(row2, "noise", 0, 0.2, 0.005, st.s, v => { st.s = v; recompute(); draw(); }, pct);
  const dS = slider(row2, "CD4−CD8 difference", 0, 0.3, 0.01, st.d, v => { st.d = v; recompute(); draw(); }, pct);
  const sg = htmlEl("div", "dc-group", "", row2);
  htmlEl("span", "", "sample:", sg);
  const sBtns = presets(sg, [
    { label: "lecture pattern", apply: () => { st.random = false; recompute(); draw(); } },
    { label: "random", apply: () => { st.random = true; recompute(); draw(); } },
  ]);
  button(sg, "new random sample", () => { st.seed += 1; sampleZ = normals(st.seed, 1)[0]; st.random = true; recompute(); draw(); });
  button(sg, "reset", () => {
    Object.assign(st, start);
    sampleZ = normals(st.seed, 1)[0];
    kS.set(st.k); lS.set(st.ll); sS.set(st.s); dS.set(st.d);
    recompute(); draw();
  });

  // ---- plots ----
  const plots = htmlEl("div", "dc-plots", "", host);
  const H = 345, WA = 430, WB = 330, WC = 390, WT = WA + WB + WC;
  const mk = (w: number, label: string) => {
    const s = svgEl("svg", { viewBox: `0 0 ${w} ${H}`, role: "img", "aria-label": label }, plots);
    s.style.width = `calc(${(100 * w / WT).toFixed(2)}% - 4px)`;
    return s;
  };
  const svgA = mk(WA, "Estimated cell-type proportions against the truth");
  const svgB = mk(WB, "Filter factors against singular value");
  const svgC = mk(WC, "Expected error: bias and noise against the amount of regularisation");
  const ro = htmlEl("div", "qb-readout dc-read", "", host);

  // ---- state derived from s and d ----
  let S: number[][] = [], sv: SVD = { d: [], U: [], V: [] }, b0: number[] = [], b: number[] = [];
  let curve: { bias: number; noise: number; total: number }[] = [];
  let kCurve: { bias: number; noise: number; total: number }[] = [];
  let lsErr = { bias: 0, noise: 0, total: 0 };
  function recompute(): void {
    S = signature(st.d);
    sv = svd(S);
    b0 = matVec(S, P_TRUE);
    b = b0.map((v, g) => v * (1 + st.s * (st.random ? sampleZ[g] : PATTERN[g])));
    curve = LGRID.map(l => expectedError(filterMatrix(sv, phis(sv, "ridge", 5, 10 ** l)), b0, P_TRUE, st.s));
    kCurve = [1, 2, 3, 4, 5].map(k => expectedError(filterMatrix(sv, phis(sv, "tsvd", k, 0)), b0, P_TRUE, st.s));
    lsErr = expectedError(filterMatrix(sv, phis(sv, "ls", 5, 0)), b0, P_TRUE, st.s);
  }

  // panel geometry
  const top = 40, bot = H - 62;
  // (a)
  const GROUPS = [...TYPES, "CD4+CD8"];
  const gx = (i: number) => 62 + (i < 5 ? i : 5.3) * 55 + 27;
  const Y0 = -0.6, Y1 = 1.2;
  const ay = scale(Y0, Y1, bot, top);
  // (b)
  const SX0 = -1, SX1 = 3.2;
  const bx = scale(SX0, SX1, 58, WB - 14), byF = scale(0, 1, bot, top + 26);
  const sxOf = (s: number) => bx(Math.max(SX0, Math.min(SX1, Math.log10(Math.max(s, 1e-30)))));
  // (c)
  const EY0 = -3, EY1 = 1;
  const cx = scale(LAM0, LAM1, 66, WC - 28), cy = scale(EY0, EY1, bot, top + 26);
  const kx = scale(5, 1, 66, WC - 28);
  const cyOf = (e: number) => cy(Math.max(EY0, Math.min(EY1, Math.log10(Math.max(e, 1e-30)))));

  function axisTitle(svg: SVGSVGElement, w: number, xlab: string, ylab: string, title: string): void {
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svg).textContent = title;
    svgEl("text", { x: (58 + w) / 2, y: H - 8, "text-anchor": "middle" }, svg).textContent = xlab;
    const yl = svgEl("text", { x: 16, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(top + bot) / 2})` }, svg);
    yl.textContent = ylab;
  }

  function draw(): void {
    const m = st.method, lam = 10 ** st.ll;
    mBtns.forEach((bt, i) => bt.classList.toggle("active", METHODS[i] === m));
    sBtns.forEach((bt, i) => bt.classList.toggle("active", (i === 1) === st.random));
    kS.input.parentElement!.classList.toggle("off", m !== "tsvd");
    lS.input.parentElement!.classList.toggle("off", m !== "ridge");
    const phi = phis(sv, m, st.k, lam);
    const f = phi.map((p, i) => p * sv.d[i]);
    const Mx = filterMatrix(sv, phi);
    const x = matVec(Mx, b);
    const ee = expectedError(Mx, b0, P_TRUE, st.s);
    const mc = MCOL[m];

    // ---------- (a) proportions ----------
    svgA.replaceChildren();
    axisTitle(svgA, WA, "", "proportion", "(a) estimate vs truth");
    for (let v = -0.4; v <= 1.2 + 1e-9; v += 0.4) {
      svgEl("line", { x1: 62, x2: WA - 8, y1: ay(v), y2: ay(v), class: Math.abs(v) < 1e-9 ? "axis" : "grid" }, svgA);
      svgEl("text", { x: 56, y: ay(v) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = minus(v.toFixed(1));
    }
    const repX = reps.map(z => matVec(Mx, b0.map((v, g) => v * (1 + st.s * z[g]))));
    const vals = (xx: number[]) => [...xx, xx[0] + xx[1]];
    const truth = vals(P_TRUE), est = vals(x);
    for (let i = 0; i < GROUPS.length; i++) {
      const c = gx(i);
      svgEl("text", { x: c, y: bot + 24, "text-anchor": "middle", class: "sm", style: i === 5 ? "font-size:16px" : "" }, svgA).textContent = GROUPS[i];
      const y0 = ay(0), ye = ay(Math.max(Y0, Math.min(Y1, est[i])));
      svgEl("rect", { x: c - 17, width: 34, y: Math.min(y0, ye), height: Math.abs(ye - y0), fill: mc, stroke: COL.ink, "stroke-width": 1.2, opacity: 0.9 }, svgA);
      let up = 0, dn = 0;
      repX.forEach((xr, r) => {
        const v = vals(xr)[i];
        if (v > Y1) { up++; return; }
        if (v < Y0) { dn++; return; }
        const jit = ((r * 0.6180339887) % 1 - 0.5) * 40;
        svgEl("circle", { cx: c + jit, cy: ay(v), r: 2.6, fill: COL.dot, opacity: 0.45 }, svgA);
      });
      if (up) svgEl("text", { x: c, y: top - 2, "text-anchor": "middle", class: "sm", style: "font-size:16px" }, svgA).textContent = `↑${up}`;
      if (dn) svgEl("text", { x: c, y: bot - 4, "text-anchor": "middle", class: "sm", style: "font-size:16px" }, svgA).textContent = `↓${dn}`;
      svgEl("line", { x1: c - 25, x2: c + 25, y1: ay(truth[i]), y2: ay(truth[i]), stroke: COL.truth, "stroke-width": 3.5, "stroke-dasharray": "7 4" }, svgA);
    }
    svgEl("line", { x1: gx(5) - 40, x2: gx(5) - 40, y1: top + 4, y2: bot, stroke: "#bbb", "stroke-dasharray": "3 3" }, svgA);
    // legend
    const lg = (x0: number, y0: number, sw: () => void, t: string) => { sw(); svgEl("text", { x: x0 + 22, y: y0 + 6, class: "sm", style: "font-size:16px" }, svgA).textContent = t; };
    const ly = H - 12;
    lg(62, ly, () => svgEl("rect", { x: 62, y: ly - 6, width: 16, height: 13, fill: mc, stroke: COL.ink }, svgA), "sample");
    lg(150, ly, () => svgEl("line", { x1: 150, x2: 168, y1: ly, y2: ly, stroke: COL.truth, "stroke-width": 3.5, "stroke-dasharray": "6 3" }, svgA), "truth");
    if (nrep) lg(232, ly, () => svgEl("circle", { cx: 240, cy: ly, r: 4, fill: COL.dot, opacity: 0.6 }, svgA), `${nrep} replicates`);

    // ---------- (b) filter factors ----------
    svgB.replaceChildren();
    axisTitle(svgB, WB, "singular value σ", "filter factor f", "(b) filter factors");
    for (const v of [0, 0.5, 1]) {
      svgEl("line", { x1: 58, x2: WB - 14, y1: byF(v), y2: byF(v), class: v === 0 ? "axis" : "grid" }, svgB);
      svgEl("text", { x: 52, y: byF(v) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = String(v);
    }
    for (let e = -1; e <= 3; e++) {
      svgEl("line", { x1: bx(e), x2: bx(e), y1: bot, y2: bot + 5, class: "axis" }, svgB);
      svgEl("text", { x: bx(e), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgB).textContent = e === 0 ? "1" : e === 1 ? "10" : `1e${e}`;
    }
    let pts = "";
    if (m === "ridge") {
      for (let i = 0; i <= 120; i++) { const l = SX0 + i * (SX1 - SX0) / 120, s = 10 ** l; pts += `${bx(l)},${byF(s * s / (s * s + lam * lam))} `; }
    } else {
      const kk = m === "tsvd" ? st.k : 5;
      const tol = Math.max(M, N) * EPS * sv.d[0];
      let kept = 0; for (let i = 0; i < kk; i++) if (sv.d[i] > tol) kept = i + 1;
      const cut = kept >= 5 ? (sv.d[4] > 0 ? Math.log10(sv.d[4]) - 0.15 : SX0) : Math.log10(Math.sqrt(sv.d[kept - 1] * Math.max(sv.d[kept], 1e-30)));
      const c = Math.max(kept >= 5 ? SX0 : SX0 + 0.1, Math.min(SX1, cut));
      pts = `${bx(SX0)},${byF(c <= SX0 ? 1 : 0)} ${bx(c)},${byF(c <= SX0 ? 1 : 0)} ${bx(c)},${byF(1)} ${bx(SX1)},${byF(1)}`;
    }
    svgEl("polyline", { points: pts, fill: "none", stroke: m === "ls" ? "#666" : mc, "stroke-width": 3.5 }, svgB);
    if (m === "ridge") {
      svgEl("line", { x1: bx(st.ll), x2: bx(st.ll), y1: top + 26, y2: bot, stroke: mc, "stroke-dasharray": "4 4", "stroke-width": 1.5 }, svgB);
      svgEl("text", { x: bx(st.ll) + 4, y: top + 22, class: "sm", "text-anchor": st.ll > 2.2 ? "end" : "start", style: "font-size:16px" }, svgB).textContent = `σ = λ`;
    }
    sv.d.forEach((s, i) => {
      const xx = sxOf(s), yy = byF(f[i]);
      svgEl("circle", { cx: xx, cy: yy, r: 7, fill: "#fff", stroke: COL.ink, "stroke-width": 2.5 }, svgB);
      if (i === 4) svgEl("text", { x: Math.max(xx, 72), y: yy + (f[i] > 0.5 ? 26 : -12), "text-anchor": "middle", class: "sm" }, svgB).textContent = s < 10 ** SX0 ? "σ₅ ≈ 0" : `σ${SUB[i]}`;
    });
    svgEl("text", { x: sxOf(Math.sqrt(sv.d[0] * sv.d[3])), y: byF(f[3]) + (f[3] > 0.5 ? 26 : -12), "text-anchor": "middle", class: "sm" }, svgB).textContent = "σ₁..σ₄";

    // ---------- (c) bias vs noise ----------
    svgC.replaceChildren();
    const tsv = m === "tsvd";
    axisTitle(svgC, WC, tsv ? "components kept k" : "ridge λ", "expected error", "(c) bias vs noise");
    for (let e = EY0; e <= EY1; e++) {
      svgEl("line", { x1: 66, x2: WC - 14, y1: cy(e), y2: cy(e), class: "grid" }, svgC);
      svgEl("text", { x: 60, y: cy(e) + 6, "text-anchor": "end", class: "sm" }, svgC).textContent = e === 0 ? "1" : e === 1 ? "10" : `1e${e}`;
    }
    svgEl("line", { x1: 66, x2: WC - 14, y1: bot, y2: bot, class: "axis" }, svgC);
    const series: [keyof typeof ee, string, string][] = [["bias", COL.bias, ""], ["noise", COL.noise, "6 4"], ["total", COL.ink, ""]];
    if (tsv) {
      for (let k = 1; k <= 5; k++) {
        svgEl("line", { x1: kx(k), x2: kx(k), y1: bot, y2: bot + 5, class: "axis" }, svgC);
        svgEl("text", { x: kx(k), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgC).textContent = String(k);
      }
      for (const [key, col, dash] of series) {
        svgEl("polyline", { points: kCurve.map((c, i) => `${kx(i + 1)},${cyOf(c[key])}`).join(" "), fill: "none", stroke: col, "stroke-width": key === "total" ? 3.5 : 3, "stroke-dasharray": dash }, svgC);
        kCurve.forEach((c, i) => svgEl("circle", { cx: kx(i + 1), cy: cyOf(c[key]), r: 3.5, fill: col }, svgC));
      }
      svgEl("circle", { cx: kx(st.k), cy: cyOf(ee.total), r: 9, fill: "none", stroke: COL.tsvd, "stroke-width": 4 }, svgC);
    } else {
      for (let e = LAM0; e <= LAM1; e++) {
        svgEl("line", { x1: cx(e), x2: cx(e), y1: bot, y2: bot + 5, class: "axis" }, svgC);
        svgEl("text", { x: cx(e), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgC).textContent = e === 0 ? "1" : e === 1 ? "10" : e === 2 ? "100" : e === 3 ? "1000" : "0.1";
      }
      for (const [key, col, dash] of series) {
        svgEl("polyline", { points: curve.map((c, i) => `${cx(LGRID[i])},${cyOf(c[key])}`).join(" "), fill: "none", stroke: col, "stroke-width": key === "total" ? 3.5 : 3, "stroke-dasharray": dash, opacity: m === "ls" ? 0.35 : 1 }, svgC);
      }
      if (m === "ridge") svgEl("circle", { cx: cx(st.ll), cy: cyOf(ee.total), r: 9, fill: "none", stroke: COL.ridge, "stroke-width": 4 }, svgC);
      else {
        svgEl("circle", { cx: 66, cy: cyOf(lsErr.total), r: 9, fill: "none", stroke: "#666", "stroke-width": 4 }, svgC);
        svgEl("text", { x: 82, y: cyOf(lsErr.total) + (lsErr.total > 0.03 ? 26 : -14), class: "sm", style: "font-size:16px" }, svgC).textContent = "least squares (λ = 0)";
      }
    }
    // legend (c)
    const lx = WC - 120;
    series.forEach(([key, col, dash], i) => {
      const yy = top + 40 + i * 22;
      svgEl("line", { x1: lx, x2: lx + 24, y1: yy, y2: yy, stroke: col, "stroke-width": 3.5, "stroke-dasharray": dash }, svgC);
      svgEl("text", { x: lx + 30, y: yy + 6, class: "sm", style: "font-size:16px" }, svgC).textContent = key;
    });

    // ---------- readout ----------
    const best = curve.reduce((bi, c, i) => (c.total < curve[bi].total ? i : bi), 0);
    const bestK = kCurve.reduce((bi, c, i) => (c.total < kCurve[bi].total ? i : bi), 0) + 1;
    const res = norm(matVec(S, x).map((v, g) => v - b[g]));
    const err = norm(x.map((v, i) => v - P_TRUE[i]));
    ro.innerHTML = "";
    const span = (h: string) => { const e = htmlEl("span", "", "", ro); e.innerHTML = h; };
    span(`<b>x</b> = (${x.map(f3).join(", ")}), CD4+CD8 = ${f3(x[0] + x[1])}`);
    span(`‖x − p‖ = ${f3(err)}, ‖Sx − b‖ = ${g3(res)}`);
    span(`σ₅ = ${g3(sv.d[4])}, κ = ${sv.d[4] > 0 ? g3(sv.d[0] / sv.d[4]) : "∞"}, f₅ = ${f3(f[4])}, gain f₅/σ₅ = ${g3(phi[4])}`);
    span(`expected: bias ${f3(ee.bias)}, noise ${f3(ee.noise)}, total ${f3(ee.total)}`);
    span(`smallest total: ridge λ ≈ ${g3(10 ** LGRID[best])} (${f3(curve[best].total)}), TSVD k = ${bestK} (${f3(kCurve[bestK - 1].total)})`);
  }

  recompute();
  draw();
}

autoMount(NAME, mount);

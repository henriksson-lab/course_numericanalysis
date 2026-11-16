// QR algorithm explorer for the eigen-algorithms unit.
// Tab "QR algorithm": A_k - mu I = Q R, A_{k+1} = R Q + mu I with deflation (exactly qr_eig of
// exercise C2), for the unit's matrices: the lecture's 3 x 3, a 5 x 5 symmetric matrix, a 4 x 4
// with two close eigenvalues, B = [[2,1],[1,2]], the switch matrix P, the predator-prey
// Jacobian J and a 4 x 4 community Jacobian with a complex pair.
//  (a) A_k as a heat map (colour = sign, strength = log10 |a_ij|; diagonal outlined, active
//      block dashed, deflated eigenvalues in green);
//  (b) |subdiagonal entries| against k on a log axis, with the predicted rates
//      |lambda_{i+1}/lambda_i|^k (dashed) when no shift is used.
// Tab "Sensitivity": a chain of n transit compartments with a feedback eps (or eps times a
// random matrix) against a symmetric matrix with the same perturbation. Eigenvalues are
// computed live (Hessenberg + Francis double shift): the chain's scatter on a circle of
// radius eps^(1/n), the symmetric matrix moves by at most eps.
// Auto-mounts on <div data-demo="qralgo">. Optional attributes:
//   data-preset="lecture|sym5|close|b2|switch|predprey|community"
//   data-shift="none|rayleigh|wilkinson"  data-k="0"
//   data-tab="qr|sens"  data-sens="1|0" (show the sensitivity tab)
//   data-n="10"  data-logeps="-10"  data-pert="entry|random"
import { svgEl, htmlEl, slider, button, presets, scale, rng, autoMount } from "../lib/ui";
import { Mat, Cx, Shift, runQR, eigvals, sortByModulus, cabs, Run } from "./math";

const NAME = "qralgo";
// Okabe-Ito, as in make-eigen-algorithms-figs.R
const OI = { orange: "#E69F00", sky: "#56B4E9", green: "#009E73", blue: "#0072B2", red: "#D55E00", purple: "#CC79A7", grey: "#999999" };
const SUBCOL = [OI.blue, OI.red, OI.green, OI.orange];
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .qa-row { display: flex; flex-wrap: wrap; gap: .35em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .qa-group { display: flex; gap: .35em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .qa-row .qb-slider { min-width: 200px; flex: 1 1 200px; max-width: 320px; grid-template-columns: auto 1fr 4.2em; }
[data-demo="${NAME}"] .qa-plots { display: flex; flex-wrap: wrap; gap: 10px; width: 100%; }
[data-demo="${NAME}"] .qa-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .qa-read { display: flex; flex-wrap: wrap; gap: .15em 1.3em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .qa-read b { font-weight: 600; }
[data-demo="${NAME}"] .qa-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .qa-ok { color: #00704f; font-weight: 600; }
[data-demo="${NAME}"] .qa-tabs .qb-btn { font-size: .95em; }
[data-demo="${NAME}"] .qa-box { display: flex; flex-direction: column; gap: .4em; width: 100%; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 18px; gap: .25em; }
.reveal [data-demo="${NAME}"] .qa-box { gap: .25em; }
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

// ---- the unit's matrices ----
interface Preset { id: string; label: string; A: Mat; note: string }
const PRESETS: Preset[] = [
  { id: "lecture", label: "3 × 3 lecture", A: [[4, 1, 0], [1, 3, 1], [0, 1, 2]], note: "the lecture's symmetric 3 × 3" },
  { id: "sym5", label: "5 × 5 symmetric", A: [[6, 2, 1, 0, 1], [2, 5, 1, 1, 0], [1, 1, 4, 1, 0], [0, 1, 1, 3, 1], [1, 0, 0, 1, 2]], note: "a symmetric 5 × 5" },
  { id: "close", label: "close λ₁, λ₂", A: [[2, 0.5, 0, 0], [0.5, 1, 0, 0.5], [0, 0, 2, 0.5], [0, 0.5, 0.5, 1]], note: "symmetric 4 × 4 with two close eigenvalues" },
  { id: "b2", label: "B (2 × 2)", A: [[2, 1], [1, 2]], note: "B of exercise C2" },
  { id: "switch", label: "switch P", A: [[0, 1], [1, 0]], note: "the two-state switch P" },
  { id: "predprey", label: "predator–prey J", A: [[-0.1, -1], [0.5, 0]], note: "the predator–prey Jacobian" },
  { id: "community", label: "4 × 4 community", A: [[-2, 0.5, -1, -1.5], [-1, -1, 1, -1], [-2, -1, -1, -2], [-1, -0.5, 1, -1.5]], note: "a nonsymmetric 4-species Jacobian" },
];
const SHIFTS: { id: Shift; label: string }[] = [
  { id: "none", label: "no shift" },
  { id: "rayleigh", label: "Rayleigh μ = aₙₙ" },
  { id: "wilkinson", label: "Wilkinson" },
];
const KMAX = 600, TOL = 1e-12;

// ---- formatting ----
const minus = (s: string): string => s.replace(/-/g, "−");
const SUP: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const SUB = "₀₁₂₃₄₅₆₇₈₉";
const sub = (i: number): string => String(i).split("").map(c => SUB[Number(c)]).join("");
const pow10 = (k: number): string => k === 0 ? "1" : "10" + String(k).split("").map(c => SUP[c]).join("");
function sci(v: number, d = 2): string {
  if (v === 0) return "0";
  if (!Number.isFinite(v)) return "—";
  const a = Math.abs(v);
  if (a >= 1e-3 && a < 1e5) return minus(v.toPrecision(d));
  const [m, e] = v.toExponential(d - 1).split("e");
  return minus(`${m}e${Number(e)}`);
}
function cell(v: number): string {
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 99.95) return minus(v.toFixed(0));
  if (a >= 0.01) { let s = v.toPrecision(3); if (s.includes(".")) s = s.replace(/0+$/, "").replace(/\.$/, ""); return minus(s); }
  const [m, e] = v.toExponential(0).split("e");
  return minus(`${m}e${Number(e)}`);
}
const fx = (v: number, d = 6): string => {
  const s = Math.abs(v) < 0.5 * 10 ** -d ? (0).toFixed(d) : v.toFixed(d);
  return minus(s);
};
function fmtC(z: Cx, d = 4): string {
  if (Math.abs(z.im) < 1e-12) return fx(z.re, d).replace(/\.?0+$/, "") || "0";
  const re = fx(z.re, d).replace(/\.?0+$/, "");
  return `${re} ± ${fx(Math.abs(z.im), d).replace(/\.?0+$/, "")}i`;
}
/** List eigenvalues (sorted by modulus), writing a conjugate pair once. */
function listEv(ev: Cx[], d = 4): string {
  const out: string[] = [];
  for (let i = 0; i < ev.length; i++) {
    if (Math.abs(ev[i].im) > 1e-12) { out.push(fmtC(ev[i], d)); if (i + 1 < ev.length && Math.abs(ev[i + 1].im + ev[i].im) < 1e-9) i++; }
    else out.push(fmtC(ev[i], d));
  }
  return out.join(",  ");
}

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const qrBox = htmlEl("div", "qa-box", "");
  const sensBox = htmlEl("div", "qa-box", "");
  if (ds.sens !== "0") {
    const tabs = htmlEl("div", "qa-row qa-tabs", "", host);
    const tb = presets(tabs, [
      { label: "QR algorithm", apply: () => show("qr") },
      { label: "Sensitivity of eigenvalues", apply: () => show("sens") },
    ]);
    const show = (t: string) => {
      qrBox.style.display = t === "qr" ? "flex" : "none";
      sensBox.style.display = t === "sens" ? "flex" : "none";
      tb.forEach((b, i) => b.classList.toggle("active", (i === 0) === (t === "qr")));
    };
    host.append(qrBox, sensBox);
    buildSens(sensBox, ds);
    show(ds.tab === "sens" ? "sens" : "qr");
  } else host.append(qrBox);
  buildQR(qrBox, ds);
}

// =====================================================================================
// Tab 1: the QR algorithm
// =====================================================================================
function buildQR(box: HTMLElement, ds: DOMStringMap): void {
  const start = {
    preset: PRESETS.find(p => p.id === ds.preset) ?? PRESETS[0],
    shift: (SHIFTS.find(s => s.id === ds.shift)?.id ?? "none") as Shift,
    k: Math.max(0, Math.min(KMAX, Math.round(Number(ds.k) || 0))),
  };
  const st = { preset: start.preset, shift: start.shift, k: start.k };

  const row1 = htmlEl("div", "qa-row", "", box);
  const pg = htmlEl("div", "qa-group", "", row1);
  htmlEl("span", "", "matrix:", pg);
  const pBtns = presets(pg, PRESETS.map(p => ({ label: p.label, apply: () => { stop(); st.preset = p; st.k = 0; recompute(); } })));

  const row3 = htmlEl("div", "qa-row", "", box);
  const sg = htmlEl("div", "qa-group", "", row3);
  htmlEl("span", "", "shift:", sg);
  const sBtns = presets(sg, SHIFTS.map(s => ({ label: s.label, apply: () => { stop(); st.shift = s.id; st.k = 0; recompute(); } })));

  const kS = slider(row3, "step k", 0, KMAX, 1, st.k, v => { stop(); st.k = v; draw(); });
  const kg = htmlEl("div", "qa-group", "", row3);
  button(kg, "k = 0", () => { stop(); setK(0); });
  button(kg, "step +1", () => { stop(); setK(st.k + 1); });
  button(kg, "+10", () => { stop(); setK(st.k + 10); });
  button(kg, "+100", () => { stop(); setK(st.k + 100); });
  const playB = button(kg, "play", () => { if (timer) stop(); else play(); });
  button(kg, "reset", () => { stop(); st.preset = start.preset; st.shift = start.shift; st.k = start.k; recompute(); });

  const plots = htmlEl("div", "qa-plots", "", box);
  const H = 360, WA = 440, WB = 660, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Heat map of the matrix A_k" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Subdiagonal entries against QR step, log scale" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro = htmlEl("div", "qb-readout qa-read", "", box);

  let timer = 0;
  function stop(): void { if (timer) { clearInterval(timer); timer = 0; } playB.textContent = "play"; }
  function play(): void {
    if (st.k >= lastK()) setK(0);
    playB.textContent = "pause";
    timer = window.setInterval(() => {
      if (st.k >= lastK() || !box.isConnected) { stop(); return; }
      setK(st.k + 1);
    }, 300);
  }
  const lastK = () => run.mats.length - 1;
  function setK(k: number): void { st.k = Math.max(0, Math.min(lastK(), k)); kS.set(st.k); draw(); }

  // ---- computed state ----
  let run: Run = runQR(st.preset.A, st.shift, KMAX, TOL);
  let ev: Cx[] = sortByModulus(eigvals(st.preset.A));
  function recompute(): void {
    run = runQR(st.preset.A, st.shift, KMAX, TOL);
    ev = sortByModulus(eigvals(st.preset.A));
    pBtns.forEach((b, i) => b.classList.toggle("active", PRESETS[i] === st.preset));
    sBtns.forEach((b, i) => b.classList.toggle("active", SHIFTS[i].id === st.shift));
    setK(st.k);
  }

  function draw(): void {
    const A = run.mats[st.k], n = A.length, m = run.active[st.k];
    drawHeat(A, n, m);
    drawConv(n);
    // ---- readout ----
    ro.replaceChildren();
    const item = (html: string, cls = "") => { const s = htmlEl("span", cls, "", ro); s.innerHTML = html; return s; };
    const mu = st.k > 0 ? run.mus[st.k - 1] : NaN;
    item(`<b>k = ${st.k}</b>${st.shift !== "none" && st.k > 0 ? `   shift μ = ${fx(mu, 6)}` : ""}`);
    item(`diagonal of A<sub>k</sub>: ${A.map((r, i) => fx(r[i], 6)).join("  ")}`);
    item(`eigenvalues: ${listEv(ev, 6)}`);
    if (st.shift === "none" && n > 1) {
      const rates = ev.slice(1).map((z, i) => cabs(z) / cabs(ev[i]));
      item(`rates |λ<sub>i+1</sub>/λ<sub>i</sub>|: ${rates.map(r => r.toFixed(3)).join("  ")}`);
    }
    const done = run.doneK >= 0 && st.k >= run.doneK;
    const nDef = n - m;
    if (done) item(`all ${n} eigenvalues found after ${run.doneK} steps (deflate when |a<sub>m,m−1</sub>| &lt; 10⁻¹²)`, "qa-ok");
    else if (nDef > 0) item(`${nDef} eigenvalue${nDef > 1 ? "s" : ""} deflated (green); active block ${m} × ${m}`, "qa-ok");
    // warnings
    const pair = ev.findIndex(z => z.im > 1e-12);
    if (pair >= 0) {
      item(`complex pair ${fmtC(ev[pair], 4)}: equal modulus, so a ${sub(pair + 2)}${sub(pair + 1)} never shrinks; real arithmetic leaves a 2 × 2 block (LAPACK reads the pair off it)`, "qa-warn");
    } else if (st.shift === "none") {
      const tie = ev.slice(1).findIndex((z, i) => Math.abs(cabs(z) - cabs(ev[i])) < 1e-9 * Math.max(1, cabs(z)));
      if (tie >= 0) item(`|λ${sub(tie + 1)}| = |λ${sub(tie + 2)}| = ${fx(cabs(ev[tie]), 4).replace(/\.?0+$/, "")}: a${sub(tie + 2)}${sub(tie + 1)} never shrinks`, "qa-warn");
    }
    if (!done && pair < 0 && st.k >= 30 && m >= 2) {
      const now = Math.abs(A[m - 1][m - 2]), before = Math.abs(run.mats[st.k - 20][m - 1][m - 2]);
      if (now > 0.5 * before && run.active[st.k - 20] === m && st.shift !== "none")
        item(`not converging: |a${sub(m)}${sub(m - 1)}| has not shrunk in 20 steps (the shift sits halfway between two eigenvalues)`, "qa-warn");
    }
    if (st.k >= lastK() && !done && lastK() === KMAX) item(`stopped at the step limit ${KMAX}`, "qa-warn");
  }

  // ---- (a) heat map ----
  function drawHeat(A: Mat, n: number, m: number): void {
    svgA.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent = `(a) A_k after k = ${st.k} step${st.k === 1 ? "" : "s"}`;
    const cw = Math.min(100, 380 / n), ch = Math.min(100, (H - 110) / n);
    const gx = 44 + (380 - cw * n) / 2, gy = 56;
    for (let i = 0; i < n; i++) {
      svgEl("text", { x: gx - 10, y: gy + (i + 0.5) * ch + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = String(i + 1);
      svgEl("text", { x: gx + (i + 0.5) * cw, y: gy - 8, "text-anchor": "middle", class: "sm" }, svgA).textContent = String(i + 1);
    }
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
      const v = A[i][j], a = Math.abs(v);
      const t = v === 0 ? 0 : Math.max(0, Math.min(1, (Math.log10(a) + 12) / 12.5));
      const x = gx + j * cw, y = gy + i * ch;
      svgEl("rect", { x, y, width: cw, height: ch, fill: "#fff", stroke: "#ccc", "stroke-width": 1 }, svgA);
      if (t > 0) svgEl("rect", { x: x + 1, y: y + 1, width: cw - 2, height: ch - 2, fill: v > 0 ? OI.blue : OI.red, "fill-opacity": (0.08 + 0.85 * t).toFixed(3) }, svgA);
      const deflated = i === j && i >= m;
      if (i === j) svgEl("rect", { x: x + 2, y: y + 2, width: cw - 4, height: ch - 4, fill: "none", stroke: deflated ? OI.green : "#222", "stroke-width": deflated ? 5 : 2.5 }, svgA);
      const txt = svgEl("text", { x: x + cw / 2, y: y + ch / 2 + 6, "text-anchor": "middle", fill: t > 0.6 ? "#fff" : "#222" }, svgA);
      txt.style.fill = t > 0.6 ? "#fff" : "#222";
      if (i === j) txt.style.fontWeight = "600";
      txt.textContent = cell(v);
    }
    if (m < n && m >= 1) svgEl("rect", { x: gx - 3, y: gy - 3, width: m * cw + 6, height: m * ch + 6, fill: "none", stroke: OI.purple, "stroke-width": 3, "stroke-dasharray": "8 5" }, svgA);
    const ly = gy + n * ch + 30;
    svgEl("rect", { x: 12, y: ly - 15, width: 18, height: 18, fill: OI.blue, "fill-opacity": 0.8 }, svgA);
    svgEl("text", { x: 36, y: ly, class: "sm" }, svgA).textContent = "+";
    svgEl("rect", { x: 56, y: ly - 15, width: 18, height: 18, fill: OI.red, "fill-opacity": 0.8 }, svgA);
    svgEl("text", { x: 80, y: ly, class: "sm" }, svgA).textContent = "−  pale: small";
    if (m < n) {
      svgEl("rect", { x: 240, y: ly + 9 - 24, width: 18, height: 18, fill: "none", stroke: OI.green, "stroke-width": 4 }, svgA);
      svgEl("text", { x: 264, y: ly, class: "sm" }, svgA).textContent = "deflated";
    }
  }

  // ---- (b) |subdiagonal| against k ----
  const Y0 = -22, Y1 = 2;
  const NICE = [10, 20, 30, 40, 60, 80, 100, 150, 200, 300, 400, 600];
  function drawConv(n: number): void {
    svgB.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgB).textContent = "(b) |subdiagonal entries| against step k";
    if (n < 2) return;
    const base = run.doneK > 0 ? run.doneK : 60;
    const xmax = NICE.find(v => v >= Math.max(base, st.k)) ?? KMAX;
    const L = 78, R = 470, top = 42, bot = H - 50;
    const bx = scale(0, xmax, L, R), by = scale(Y0, Y1, bot, top);
    const ly = (v: number) => by(Math.max(Y0, Math.min(Y1, Math.log10(Math.max(v, 1e-300)))));
    for (let e = Y0 + 2; e <= 0; e += 4) {
      svgEl("line", { x1: L, x2: R, y1: by(e), y2: by(e), class: "grid" }, svgB);
      svgEl("text", { x: L - 6, y: by(e) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = pow10(e);
    }
    svgEl("line", { x1: L, x2: R, y1: bot, y2: bot, class: "axis" }, svgB);
    svgEl("line", { x1: L, x2: L, y1: top, y2: bot, class: "axis" }, svgB);
    const step = xmax <= 20 ? 5 : xmax <= 60 ? 10 : xmax <= 100 ? 20 : xmax <= 200 ? 50 : 100;
    for (let k = 0; k <= xmax; k += step) {
      svgEl("line", { x1: bx(k), x2: bx(k), y1: bot, y2: bot + 5, class: "axis" }, svgB);
      svgEl("text", { x: bx(k), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgB).textContent = String(k);
    }
    svgEl("text", { x: (L + R) / 2, y: H - 6, "text-anchor": "middle" }, svgB).textContent = "QR step k";
    // deflation tolerance
    svgEl("line", { x1: L, x2: R, y1: by(-12), y2: by(-12), stroke: "#555", "stroke-width": 1.5, "stroke-dasharray": "2 4" }, svgB);
    svgEl("text", { x: R - 4, y: by(-12) - 6, "text-anchor": "end", class: "sm" }, svgB).textContent = "deflate below 10⁻¹²";
    // current k
    svgEl("line", { x1: bx(Math.min(st.k, xmax)), x2: bx(Math.min(st.k, xmax)), y1: top, y2: bot, stroke: "#bbb", "stroke-width": 1.5 }, svgB);

    const unshifted = st.shift === "none";
    const nsub = n - 1;
    for (let i = 0; i < nsub; i++) {
      const col = SUBCOL[i % SUBCOL.length];
      if (unshifted) {   // predicted |lambda_{i+1}/lambda_i|^k
        const r = cabs(ev[i + 1]) / cabs(ev[i]);
        const kEnd = r < 1 ? Math.min(xmax, (Y0 / Math.log10(r))) : xmax;
        svgEl("line", { x1: bx(0), y1: by(0), x2: bx(kEnd), y2: by(r < 1 ? kEnd * Math.log10(r) : 0), stroke: col, "stroke-width": 2, "stroke-dasharray": "7 5", opacity: 0.8 }, svgB);
      }
      const kEnd = Math.min(st.k, run.defK[i] >= 0 ? run.defK[i] : st.k);
      const pts: string[] = [];
      for (let k = 0; k <= kEnd; k++) pts.push(`${bx(k).toFixed(1)},${ly(Math.abs(run.mats[k][i + 1][i])).toFixed(1)}`);
      svgEl("polyline", { points: pts.join(" "), fill: "none", stroke: col, "stroke-width": 2.5 }, svgB);
      if (xmax <= 100) for (let k = 0; k <= kEnd; k++)
        svgEl("circle", { cx: bx(k), cy: ly(Math.abs(run.mats[k][i + 1][i])), r: 3.5, fill: col }, svgB);
      if (run.defK[i] >= 0 && run.defK[i] <= st.k)
        svgEl("circle", { cx: bx(run.defK[i]), cy: ly(Math.abs(run.mats[run.defK[i]][i + 1][i])), r: 9, fill: "none", stroke: OI.green, "stroke-width": 3 }, svgB);
    }
    // legend
    const lx = R + 18;
    let yy = top + 14;
    for (let i = 0; i < nsub; i++) {
      const col = SUBCOL[i % SUBCOL.length];
      svgEl("line", { x1: lx, x2: lx + 26, y1: yy - 6, y2: yy - 6, stroke: col, "stroke-width": 3 }, svgB);
      svgEl("text", { x: lx + 34, y: yy }, svgB).textContent = `|a${sub(i + 2)}${sub(i + 1)}|`;
      if (unshifted) {
        yy += 24;
        svgEl("line", { x1: lx, x2: lx + 26, y1: yy - 6, y2: yy - 6, stroke: col, "stroke-width": 2, "stroke-dasharray": "7 5" }, svgB);
        const r = cabs(ev[i + 1]) / cabs(ev[i]);
        svgEl("text", { x: lx + 34, y: yy, class: "sm" }, svgB).textContent = `${r.toFixed(3)}ᵏ`;
      }
      yy += 30;
    }
    svgEl("circle", { cx: lx + 13, cy: yy - 6, r: 9, fill: "none", stroke: OI.green, "stroke-width": 3 }, svgB);
    svgEl("text", { x: lx + 34, y: yy, class: "sm" }, svgB).textContent = "deflated";
  }

  recompute();
}

// =====================================================================================
// Tab 2: sensitivity of eigenvalues
// =====================================================================================
type Pert = "entry" | "random";
function chainMat(n: number): Mat {
  return Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? -1 : i === j + 1 ? 1 : 0)));
}
function symMat(n: number): Mat {
  return Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? -1 : Math.abs(i - j) === 1 ? 0.5 : 0)));
}
function perturbed(n: number, eps: number, pert: Pert): { K: Mat; S: Mat } {
  const K = chainMat(n), S = symMat(n);
  if (pert === "entry") {
    K[0][n - 1] += eps;
    S[0][n - 1] += eps; S[n - 1][0] += eps;
    if (n === 1) S[0][0] -= eps;
  } else {
    const g = rng(1);
    const G = Array.from({ length: n }, () => Array.from({ length: n }, () => g.norm()));
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
      K[i][j] += eps * G[i][j];
      S[i][j] += eps * (G[i][j] + G[j][i]) / 2;
    }
  }
  return { K, S };
}
const sortedRe = (ev: Cx[]): number[] => ev.map(z => z.re).sort((a, b) => a - b);
function moves(n: number, eps: number, pert: Pert, s0: number[]): { evK: Cx[]; evS: Cx[]; dK: number; dS: number } {
  const { K, S } = perturbed(n, eps, pert);
  const evK = eigvals(K), evS = eigvals(S);
  const dK = Math.max(...evK.map(z => Math.hypot(z.re + 1, z.im)));
  const sS = sortedRe(evS);
  const dS = Math.max(...sS.map((v, i) => Math.abs(v - s0[i])));
  return { evK, evS, dK, dS };
}

function buildSens(box: HTMLElement, ds: DOMStringMap): void {
  const n0 = Math.max(2, Math.min(20, Math.round(Number(ds.n) || 10)));
  const le0 = Number.isFinite(Number(ds.logeps)) && ds.logeps !== undefined ? Math.max(-16, Math.min(-1, Number(ds.logeps))) : -10;
  const start = { n: n0, le: le0, pert: (ds.pert === "random" ? "random" : "entry") as Pert };
  const st = { ...start };

  const row1 = htmlEl("div", "qa-row", "", box);
  const pg = htmlEl("div", "qa-group", "", row1);
  htmlEl("span", "", "perturbation:", pg);
  const pBtns = presets(pg, [
    { label: "feedback ε in entry (1, n)", apply: () => { st.pert = "entry"; draw(); } },
    { label: "ε × random matrix", apply: () => { st.pert = "random"; draw(); } },
  ]);
  const nS = slider(row1, "compartments n", 2, 20, 1, st.n, v => { st.n = v; draw(); });
  const row2 = htmlEl("div", "qa-row", "", box);
  const eS = slider(row2, "ε", -16, -1, 0.5, st.le, v => { st.le = v; draw(); }, v => Number.isInteger(v) ? `1e${v}` : sci(10 ** v, 2));
  const eg = htmlEl("div", "qa-group", "", row2);
  const eVals = [-16, -12, -10, -4];
  const eBtns = presets(eg, eVals.map(l => ({ label: "ε = " + pow10(l), apply: () => { st.le = l; eS.set(l); draw(); } })));
  button(eg, "reset", () => { Object.assign(st, start); nS.set(st.n); eS.set(st.le); draw(); });

  const plots = htmlEl("div", "qa-plots", "", box);
  const H = 370, WA = 420, WB = 680, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Eigenvalues in the complex plane" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Eigenvalue change against perturbation size, log-log" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro = htmlEl("div", "qb-readout qa-read", "", box);

  const CK = OI.blue, CS = OI.red;

  function draw(): void {
    pBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.pert === "entry")));
    eBtns.forEach((b, i) => b.classList.toggle("active", Math.abs(eVals[i] - st.le) < 1e-9));
    const n = st.n, eps = 10 ** st.le;
    const s0 = sortedRe(eigvals(symMat(n)));
    const cur = moves(n, eps, st.pert, s0);
    const rad = eps ** (1 / n);

    // ---------- (a) complex plane ----------
    svgA.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent = "(a) eigenvalues, complex plane";
    const L = 62, top = 40, side = 284;
    const ax = scale(-2.2, 0.2, L, L + side), ay = scale(-1.2, 1.2, top + side, top);
    for (const v of [-2, -1.5, -1, -0.5, 0]) {
      svgEl("line", { x1: ax(v), x2: ax(v), y1: top, y2: top + side, class: v === 0 ? "axis" : "grid" }, svgA);
      if (Number.isInteger(v)) svgEl("text", { x: ax(v), y: top + side + 22, "text-anchor": "middle", class: "sm" }, svgA).textContent = minus(String(v));
    }
    for (const v of [-1, -0.5, 0, 0.5, 1]) {
      svgEl("line", { x1: L, x2: L + side, y1: ay(v), y2: ay(v), class: v === 0 ? "axis" : "grid" }, svgA);
      if (Number.isInteger(v)) svgEl("text", { x: L - 6, y: ay(v) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = v === 0 ? "0" : minus(String(v)) + "i";
    }
    svgEl("text", { x: L + side / 2, y: H - 4, "text-anchor": "middle" }, svgA).textContent = "real part";
    const sc = side / 2.4;
    svgEl("circle", { cx: ax(-1), cy: ay(0), r: Math.min(rad * sc, 3 * side), fill: "none", stroke: CK, "stroke-width": 2, "stroke-dasharray": "6 5" }, svgA);
    for (const v of s0) svgEl("circle", { cx: ax(v), cy: ay(0), r: 9, fill: "#fff", stroke: CS, "stroke-width": 2.5 }, svgA);
    for (const z of cur.evS) svgEl("circle", { cx: ax(z.re), cy: ay(z.im), r: 4, fill: CS }, svgA);
    const xs = 8;
    svgEl("path", { d: `M${ax(-1) - xs},${ay(0) - xs}L${ax(-1) + xs},${ay(0) + xs}M${ax(-1) - xs},${ay(0) + xs}L${ax(-1) + xs},${ay(0) - xs}`, stroke: "#222", "stroke-width": 3 }, svgA);
    for (const z of cur.evK) {
      const x = Math.max(L - 4, Math.min(L + side + 4, ax(z.re))), y = Math.max(top - 4, Math.min(top + side + 4, ay(z.im)));
      svgEl("circle", { cx: x, cy: y, r: 6, fill: CK, stroke: "#fff", "stroke-width": 1 }, svgA);
    }

    // ---------- (b) change against eps ----------
    svgB.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgB).textContent = "(b) how far the eigenvalues move";
    const BL = 78, BR = 470, bt = 44, bb = H - 50, Y0 = -17, Y1 = 0;
    const bx = scale(-16, -1, BL, BR), by = scale(Y0, Y1, bb, bt);
    const ly = (v: number) => by(Math.max(Y0, Math.min(Y1, Math.log10(Math.max(v, 1e-300)))));
    for (let e = Y0 + 1; e <= Y1; e += 4) {
      svgEl("line", { x1: BL, x2: BR, y1: by(e), y2: by(e), class: "grid" }, svgB);
      svgEl("text", { x: BL - 6, y: by(e) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = pow10(e);
    }
    for (let e = -16; e <= -1; e += 3) {
      svgEl("line", { x1: bx(e), x2: bx(e), y1: bt, y2: bb, class: "grid" }, svgB);
      svgEl("text", { x: bx(e), y: bb + 22, "text-anchor": "middle", class: "sm" }, svgB).textContent = pow10(e);
    }
    svgEl("line", { x1: BL, x2: BR, y1: bb, y2: bb, class: "axis" }, svgB);
    svgEl("line", { x1: BL, x2: BL, y1: bt, y2: bb, class: "axis" }, svgB);
    svgEl("text", { x: (BL + BR) / 2, y: H - 8, "text-anchor": "middle" }, svgB).textContent = "perturbation ε";
    const yl = svgEl("text", { x: 16, y: (bt + bb) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(bt + bb) / 2})` }, svgB);
    yl.textContent = "largest change";
    // predicted lines
    svgEl("line", { x1: bx(-16), y1: by(-16 / n), x2: bx(-1), y2: by(-1 / n), stroke: CK, "stroke-width": 2, "stroke-dasharray": "7 5" }, svgB);
    svgEl("line", { x1: bx(-16), y1: by(-16), x2: bx(-1), y2: by(-1), stroke: CS, "stroke-width": 2, "stroke-dasharray": "7 5" }, svgB);
    // computed points
    svgEl("line", { x1: bx(st.le), x2: bx(st.le), y1: bt, y2: bb, stroke: "#bbb", "stroke-width": 1.5 }, svgB);
    for (let e = -16; e <= -1; e++) {
      const r = moves(n, 10 ** e, st.pert, s0);
      for (const [d, c] of [[r.dK, CK], [r.dS, CS]] as [number, string][]) {
        const zero = d < 10 ** Y0;
        svgEl("circle", { cx: bx(e), cy: ly(d), r: 5, fill: zero ? "#fff" : c, stroke: c, "stroke-width": 2 }, svgB);
      }
    }
    for (const [d, c] of [[cur.dK, CK], [cur.dS, CS]] as [number, string][])
      svgEl("circle", { cx: bx(st.le), cy: ly(d), r: 10, fill: "none", stroke: c, "stroke-width": 3 }, svgB);
    // legend
    const lx = BR + 16;
    let yy = bt + 14;
    const leg = (draw: (y: number) => void, text: string, cls = "") => {
      draw(yy - 6);
      const t = svgEl("text", { x: lx + 34, y: yy, class: cls }, svgB); t.textContent = text; yy += 26;
    };
    leg(y => svgEl("circle", { cx: lx + 13, cy: y, r: 6, fill: CK }, svgB), "chain");
    leg(y => svgEl("line", { x1: lx, x2: lx + 26, y1: y, y2: y, stroke: CK, "stroke-width": 2, "stroke-dasharray": "7 5" }, svgB), `ε^(1/${n})`, "sm");
    yy += 8;
    leg(y => svgEl("circle", { cx: lx + 13, cy: y, r: 6, fill: CS }, svgB), "symmetric");
    leg(y => svgEl("line", { x1: lx, x2: lx + 26, y1: y, y2: y, stroke: CS, "stroke-width": 2, "stroke-dasharray": "7 5" }, svgB), "ε", "sm");
    yy += 8;
    leg(y => svgEl("circle", { cx: lx + 13, cy: y, r: 5, fill: "#fff", stroke: "#555", "stroke-width": 2 }, svgB), "moved < 10⁻¹⁷", "sm");
    yy += 14;
    for (const [t, c] of [["chain: all −1;", CK], ["drug moves 1 → 2", CK], ["→ … → n", CK], ["symmetric: −1 on", CS], ["the diagonal, 0.5", CS], ["on both sides", CS]] as [string, string][]) {
      const tt = svgEl("text", { x: lx, y: yy, class: "sm" }, svgB); tt.textContent = t; tt.style.fill = c; yy += 22;
    }

    // ---------- readout ----------
    ro.replaceChildren();
    const item = (html: string, cls = "") => { const s = htmlEl("span", cls, "", ro); s.innerHTML = html; };
    const epsTxt = Number.isInteger(st.le) ? pow10(st.le) : sci(eps, 2);
    item(`<b>n = ${n}, ε = ${epsTxt}</b>`);
    item(`<span style="color:${CK}">chain</span>: eigenvalues moved by up to <b>${sci(cur.dK, 3)}</b>  (ε<sup>1/n</sup> = ${sci(rad, 3)})`);
    item(`<span style="color:${CS}">symmetric</span>: moved by up to <b>${sci(cur.dS, 2)}</b>`);
    if (cur.dK > 0 && cur.dS > 0) item(`chain / symmetric: ${sci(cur.dK / cur.dS, 2)} ×`);
    if (st.pert === "entry" && cur.dK < 1e-6 * rad)
      item(`ε is at the level of rounding: the solver returned −1 exactly, but any point inside the dashed circle is an equally good answer`, "qa-warn");
  }
  draw();
}

autoMount(NAME, mount);

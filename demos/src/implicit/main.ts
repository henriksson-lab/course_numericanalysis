// Implicit vs explicit for the diffusion-implicit unit. Two tabs.
//
// "Race": the chamber of the previous unit (L = 500 µm, D = 1000 µm²/s, N = 50, dx = 10 µm)
// solved by FTCS, backward Euler and Crank–Nicolson with ONE shared time step, chosen by
// the mesh ratio r = D dt / dx² from 0.1 to 100 (dt = 0.01 ... 10 s).
//  initial data "sharp start": source c = 1 at x = 0, sink c = 0 at x = L, empty channel
//                (exact: Fourier series); error measured at t = 60 s, as on the slides;
//  initial data "smooth sine": both ends 0, c = sin(pi x / L); exact e^{-D k² t} sin(kx),
//                k = pi / L; error measured at t = 30 s, as on the slides.
//  (a) the profile at time t = n dt (n steps): exact line, the three schemes;
//  (b) log–log max error at the reference time against dt for the three schemes, with
//      slope-1 and slope-2 guides, the FTCS limit r = 1/2 and a marker at the current dt.
//  Readout: zigzag growth factors, steps = Thomas solves, operation counts, and the steps
//  FTCS would need at r = 1/2.
//
// "Embryo": the 1D embryo between a 20 °C and a 27 °C stream (N = 200, dx = 2.5 µm,
// D = 1.4e5 µm²/s, tau = R² / (2.405² D), no flux at the tips, T_out = 23.5 at the middle
// node). "Streams on" starts from 20 °C everywhere; "flow stopped" starts from the sharp
// 20|27 step without exchange (exact: cosine series). Backward Euler or Crank–Nicolson.
//  (c) the embryo between the streams, coloured by its temperature, and the profile T(x)
//      with the stream temperature, the steady formula and the 10–90 % band;
//  (d) T(t) at x = 125 and 375 µm, with tau and 4 tau marked.
//
// Auto-mounts on <div data-demo="implicit">. Optional attributes:
//   data-tab="race|embryo"  data-ic="step|sine"  data-r="40"  data-embryo="1|0"
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, autoMount } from "../lib/ui";

const NAME = "implicit";
// Okabe-Ito, as in make-diffusion-implicit-figs.R: FTCS orange, BE blue, CN vermillion
const COL = {
  ftcs: "#E69F00", be: "#0072B2", cn: "#D55E00", exact: "#222222", grey: "#999999",
  green: "#009E73", yellow: "#F0E442", purple: "#CC79A7", sky: "#56B4E9", warn: "#b00020", ink: "#222",
};
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .35em; }
[data-demo="${NAME}"] .im-row { display: flex; flex-wrap: wrap; gap: .3em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .im-group { display: flex; gap: .35em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .im-row .qb-slider { flex: 1 1 380px; min-width: 280px; max-width: 540px; grid-template-columns: auto 1fr 7.5em; }
[data-demo="${NAME}"] .im-plots { display: flex; flex-wrap: wrap; gap: 8px; width: 100%; }
[data-demo="${NAME}"] .im-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .im-read { width: 100%; box-sizing: border-box; display: flex; flex-wrap: wrap; gap: .1em 1.3em; }
[data-demo="${NAME}"] .im-read b { font-weight: 600; }
[data-demo="${NAME}"] .im-warn { color: ${COL.warn}; font-weight: 600; }
[data-demo="${NAME}"] .im-tabs .qb-btn { font-size: .95em; }
[data-demo="${NAME}"] .qb-btn.im-play { min-width: 5.5em; font-weight: 600; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; padding: .2em .5em; }
.reveal [data-demo="${NAME}"].qb-demo .qb-btn { font-size: .9em; padding: .12em .5em; }
.reveal [data-demo="${NAME}"].qb-demo .qb-slider output { font-size: 18px; font-weight: 600; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---------------------------------------------------------------- numerics

/** Thomas algorithm: a sub-diagonal (a[0] unused), b diagonal, c super-diagonal (c[m-1] unused). */
function thomas(a: Float64Array, b0: Float64Array, c: Float64Array, d0: Float64Array): Float64Array {
  const n = d0.length, b = Float64Array.from(b0), d = Float64Array.from(d0);
  for (let i = 1; i < n; i++) {
    const w = a[i] / b[i - 1];
    b[i] -= w * c[i - 1];
    d[i] -= w * d[i - 1];
  }
  const x = new Float64Array(n);
  x[n - 1] = d[n - 1] / b[n - 1];
  for (let i = n - 2; i >= 0; i--) x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
  return x;
}

type Scheme = "ftcs" | "be" | "cn";
const SCHEMES: Scheme[] = ["ftcs", "be", "cn"];
const LABEL: Record<Scheme, string> = { ftcs: "FTCS", be: "backward Euler", cn: "Crank–Nicolson" };
const OPS: Record<Scheme, number> = { ftcs: 5, be: 8, cn: 13 };   // operations per grid point per step

/** One step of the chosen scheme on c (including the two fixed boundary values), in place. */
function step(s: Scheme, c: Float64Array, r: number): void {
  const n = c.length - 1, m = n - 1;
  if (s === "ftcs") {
    let prev = c[0];
    for (let i = 1; i < n; i++) {
      const ci = c[i];
      c[i] = ci + r * (prev - 2 * ci + c[i + 1]);
      prev = ci;
    }
    return;
  }
  const th = s === "be" ? 1 : 0.5, rt = th * r, re = (1 - th) * r;
  const a = new Float64Array(m).fill(-rt), b = new Float64Array(m).fill(1 + 2 * rt), cc = new Float64Array(m).fill(-rt);
  const d = new Float64Array(m);
  for (let i = 1; i < n; i++) d[i - 1] = c[i] + re * (c[i - 1] - 2 * c[i] + c[i + 1]);
  d[0] += rt * c[0];
  d[m - 1] += rt * c[n];
  const x = thomas(a, b, cc, d);
  for (let i = 1; i < n; i++) c[i] = x[i - 1];
}

// The chamber
const D = 1000, L = 500, N = 50, DX = L / N;
type IC = "step" | "sine";
const TREF: Record<IC, number> = { step: 60, sine: 30 };
const XS = Array.from({ length: N + 1 }, (_, i) => i * DX);

function initial(ic: IC): Float64Array {
  const c = new Float64Array(N + 1);
  if (ic === "step") c[0] = 1;
  else for (let i = 1; i < N; i++) c[i] = Math.sin(Math.PI * XS[i] / L);
  return c;
}
/** Exact solution at x, t. */
function exact(ic: IC, x: number, t: number): number {
  if (ic === "sine") return Math.exp(-D * (Math.PI / L) ** 2 * t) * Math.sin(Math.PI * x / L);
  if (x <= 0) return 1;
  if (x >= L) return 0;
  if (t <= 0) return 0;
  let s = 0;
  for (let n = 1; n <= 20000; n++) {
    const e = Math.exp(-D * (n * Math.PI / L) ** 2 * t);
    if (e < 1e-18) break;
    s += 2 / (n * Math.PI) * Math.sin(n * Math.PI * x / L) * e;
  }
  return 1 - x / L - s;
}
function maxErr(ic: IC, c: Float64Array, t: number): number {
  let m = 0;
  for (let i = 0; i <= N; i++) {
    const e = Math.abs(c[i] - exact(ic, XS[i], t));
    if (!Number.isFinite(e)) return Infinity;
    if (e > m) m = e;
  }
  return m;
}

// dt grid: log10 r from -1 to 2
const LR0 = -1, LR1 = 2, LRS = 0.05;
const LRGRID = Array.from({ length: Math.round((LR1 - LR0) / LRS) + 1 }, (_, i) => LR0 + i * LRS);
const rOf = (lr: number): number => 10 ** lr;
const dtOf = (r: number): number => r * DX * DX / D;

interface ErrCurve { dt: number[]; err: Record<Scheme, number[]> }
const errCache: Partial<Record<IC, ErrCurve>> = {};
/** Effective dt and the number of steps for the reference time. */
function refSteps(ic: IC, dt: number): [number, number] {
  const n = Math.max(1, Math.round(TREF[ic] / dt));
  return [n, TREF[ic] / n];
}
function errCurve(ic: IC): ErrCurve {
  const hit = errCache[ic];
  if (hit) return hit;
  const out: ErrCurve = { dt: [], err: { ftcs: [], be: [], cn: [] } };
  for (const lr of LRGRID) {
    const [n, dte] = refSteps(ic, dtOf(rOf(lr)));
    out.dt.push(dte);
    for (const s of SCHEMES) out.err[s].push(errAt(ic, s, n, dte));
  }
  return (errCache[ic] = out);
}
function errAt(ic: IC, s: Scheme, n: number, dte: number): number {
  const c = initial(ic), r = D * dte / (DX * DX);
  for (let k = 0; k < n; k++) {
    step(s, c, r);
    if (s === "ftcs" && !Number.isFinite(c[1])) return Infinity;
  }
  return maxErr(ic, c, TREF[ic]);
}

// ---------------------------------------------------------------- formatting
function sci(v: number, d = 2): string {
  if (!Number.isFinite(v)) return v > 0 ? "∞" : "NaN";
  if (v === 0) return "0";
  if (Math.abs(v) >= 1e4 || Math.abs(v) < 1e-3) return v.toExponential(d - 1).replace("e+", "e");
  return String(+v.toPrecision(d));
}
const fmtInt = (n: number): string => Math.round(n).toString().replace(/\B(?=(\d{3})+(?!\d))/g, " ");
const sup: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
const pow10 = (k: number): string => k === 0 ? "1" : k === 1 ? "10" : "10" + String(k).split("").map(c => sup[c]).join("");
const fmtR = (r: number): string => r >= 10 ? r.toFixed(0) : r >= 1 ? r.toFixed(1) : r.toFixed(2);
const fmtDt = (dt: number): string => dt >= 1 ? `${+dt.toPrecision(3)} s` : `${+(dt * 1000).toPrecision(3)} ms`;

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const uid = ++mountCount;
  const withEmbryo = ds.embryo !== "0";

  const raceBox = htmlEl("div", "", "");
  const embBox = htmlEl("div", "", "");
  for (const b of [raceBox, embBox]) Object.assign(b.style, { display: "flex", flexDirection: "column", gap: ".35em", width: "100%" });
  let stopAll: () => void = () => {};
  if (withEmbryo) {
    const tabs = htmlEl("div", "im-row im-tabs", "", host);
    const tb = presets(tabs, [
      { label: "Race: FTCS, backward Euler, Crank–Nicolson", apply: () => show("race") },
      { label: "Embryo between two streams", apply: () => show("embryo") },
    ]);
    const show = (t: string) => {
      stopAll();
      raceBox.style.display = t === "race" ? "flex" : "none";
      embBox.style.display = t === "embryo" ? "flex" : "none";
      tb.forEach((b, i) => b.classList.toggle("active", (i === 0) === (t === "race")));
    };
    host.append(raceBox, embBox);
    const stopR = buildRace(raceBox, ds, uid);
    const stopE = buildEmbryo(embBox, uid);
    stopAll = () => { stopR(); stopE(); };
    show(ds.tab === "embryo" ? "embryo" : "race");
  } else {
    host.append(raceBox);
    buildRace(raceBox, ds, uid);
  }
}

// ---------------------------------------------------------------- race tab
function buildRace(box: HTMLElement, ds: DOMStringMap, uid: number): () => void {
  const r0 = Number(ds.r);
  const START = {
    ic: (ds.ic === "sine" ? "sine" : "step") as IC,
    lr: Number.isFinite(r0) && r0 > 0 ? Math.min(LR1, Math.max(LR0, Math.log10(r0))) : Math.log10(40),
  };
  const st = { ic: START.ic, lr: START.lr, t: TREF[START.ic], show: { ftcs: true, be: true, cn: true } as Record<Scheme, boolean> };

  // incremental simulation state
  let sim: { key: string; n: number; c: Record<Scheme, Float64Array>; blownAt: Record<Scheme, number> } | null = null;
  const dtNow = () => dtOf(rOf(st.lr));
  const advance = () => {
    const dt = dtNow(), r = rOf(st.lr), key = `${st.ic}|${st.lr.toFixed(3)}`;
    const target = Math.floor(st.t / dt + 1e-9);
    if (!sim || sim.key !== key || sim.n > target) {
      sim = { key, n: 0, c: { ftcs: initial(st.ic), be: initial(st.ic), cn: initial(st.ic) }, blownAt: { ftcs: -1, be: -1, cn: -1 } };
    }
    while (sim.n < target) {
      sim.n++;
      for (const s of SCHEMES) {
        if (sim.blownAt[s] >= 0) continue;
        step(s, sim.c[s], r);
        if (!Number.isFinite(sim.c[s][N >> 1]) || Math.abs(sim.c[s][N >> 1]) > 1e150) sim.blownAt[s] = sim.n;
      }
    }
    return sim;
  };

  // ---- controls ----
  const row1 = htmlEl("div", "im-row", "", box);
  const g1 = htmlEl("div", "im-group", "", row1);
  htmlEl("span", "", "start:", g1);
  const icBtns = presets(g1, [
    { label: "sharp start (chamber)", apply: () => { setIC("step"); } },
    { label: "smooth sine", apply: () => { setIC("sine"); } },
  ]);
  const rS = slider(row1, "time step", LR0, LR1, LRS, st.lr, v => { st.lr = v; draw(); },
    v => `r=${fmtR(rOf(v))}, ${fmtDt(dtOf(rOf(v)))}`);
  const row2 = htmlEl("div", "im-row", "", box);
  const g2 = htmlEl("div", "im-group", "", row2);
  htmlEl("span", "", "r =", g2);
  const rVals = [0.4, 0.5, 0.6, 1, 4, 10, 40, 100];
  const rBtns = presets(g2, rVals.map(r => ({ label: String(r), apply: () => { st.lr = Math.log10(r); rS.set(st.lr); draw(); } })));
  const tS = slider(row2, "time t", 0, 120, 0.5, st.t, v => { stop(); st.t = v; draw(); }, v => `${v.toFixed(1)} s`);
  const g3 = htmlEl("div", "im-group", "", row2);
  const playB = button(g3, "▶ play", () => { if (timer) stop(); else play(); }, "im-play");
  button(g3, "reset", () => {
    stop(); st.ic = START.ic; st.lr = START.lr; st.t = TREF[START.ic]; rS.set(st.lr); tS.set(st.t);
    for (const s of SCHEMES) { st.show[s] = true; checks[s].checked = true; }
    draw();
  });
  const row3 = htmlEl("div", "im-row", "", box);
  const checks = {} as Record<Scheme, HTMLInputElement>;
  for (const s of SCHEMES) {
    checks[s] = checkbox(row3, LABEL[s], true, v => { st.show[s] = v; draw(); });
    (checks[s].parentElement as HTMLElement).style.color = COL[s];
    (checks[s].parentElement as HTMLElement).style.fontWeight = "600";
  }

  function setIC(ic: IC): void {
    stop(); st.ic = ic; st.t = TREF[ic]; tS.set(st.t); draw();
  }

  // ---- animation ----
  let timer = 0, last = 0;
  const RATE = 15;   // simulated seconds per real second
  function play(): void {
    if (st.t >= 119.99) { st.t = 0; tS.set(0); }
    playB.textContent = "❚❚ pause";
    last = performance.now();
    const tick = (now: number) => {
      st.t = Math.min(120, st.t + RATE * (now - last) / 1000); last = now;
      tS.set(Math.round(st.t * 2) / 2);
      draw();
      if (st.t >= 120) { stop(); return; }
      timer = requestAnimationFrame(tick);
    };
    timer = requestAnimationFrame(tick);
  }
  function stop(): void {
    if (timer) cancelAnimationFrame(timer);
    timer = 0; playB.textContent = "▶ play";
  }

  // ---- plots ----
  const plots = htmlEl("div", "im-plots", "", box);
  const H = 330, WA = 560, WB = 540, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Concentration profile" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 5px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Error against time step, log-log" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 5px)`;
  const ro = htmlEl("div", "qb-readout im-read", "", box);

  const top = 36, bot = H - 50, la = 62;
  const ax = scale(0, L, la, WA - 14);
  const LT0 = -2.2, LT1 = 1.2, LE0 = -7, LE1 = 0.5;
  const bx = scale(LT0, LT1, 72, WB - 14), by = scale(LE0, LE1, bot, top);

  function draw(): void {
    const s0 = advance();
    const dt = dtNow(), r = rOf(st.lr), n = s0.n, tn = n * dt, ic = st.ic;
    icBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (ic === "step")));
    rBtns.forEach((b, i) => b.classList.toggle("active", Math.abs(Math.log10(rVals[i]) - st.lr) < 1e-6));

    // ---------- (a) profile ----------
    svgA.replaceChildren();
    const y0 = ic === "step" ? -0.4 : -0.4, y1 = ic === "step" ? 1.4 : 1.4;
    const ay = scale(y0, y1, bot, top);
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent = `(a) profile at t = ${tn.toFixed(tn < 10 ? 2 : 1)} s  (${fmtInt(n)} steps)`;
    for (let v = -0.4; v <= 1.4 + 1e-9; v += 0.2) {
      const vv = Math.round(v * 10) / 10;
      svgEl("line", { x1: la, x2: WA - 14, y1: ay(vv), y2: ay(vv), class: vv === 0 ? "axis" : "grid" }, svgA);
      if (Math.abs(vv * 5 - Math.round(vv * 5)) < 1e-9 && Math.round(vv * 10) % 4 === 0)
        svgEl("text", { x: la - 6, y: ay(vv) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = vv.toFixed(1);
    }
    for (let x = 0; x <= L; x += 100) {
      svgEl("line", { x1: ax(x), x2: ax(x), y1: bot, y2: bot + 5, class: "axis" }, svgA);
      svgEl("text", { x: ax(x), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgA).textContent = String(x);
    }
    svgEl("text", { x: (la + WA) / 2, y: H - 6, "text-anchor": "middle" }, svgA).textContent = "position x (µm)";
    const yl = svgEl("text", { x: 18, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 18 ${(top + bot) / 2})` }, svgA);
    yl.textContent = "concentration c";
    const cp = svgEl("clipPath", { id: `im-ca-${uid}` }, svgA);
    svgEl("rect", { x: la, y: top, width: WA - 14 - la, height: bot - top }, cp);
    const pg = svgEl("g", { "clip-path": `url(#im-ca-${uid})` }, svgA);
    // exact
    const ex = Array.from({ length: 201 }, (_, i) => i * L / 200);
    const pts = ex.map(x => `${ax(x)},${ay(tn === 0 && ic === "step" ? (x === 0 ? 1 : 0) : exact(ic, x, tn))}`).join(" ");
    svgEl("polyline", { points: pts, fill: "none", stroke: COL.exact, "stroke-width": 3, "stroke-opacity": 0.55 }, pg);
    const clampV = (v: number) => Number.isFinite(v) ? Math.max(y0 - 1, Math.min(y1 + 1, v)) : (v > 0 ? y1 + 1 : y0 - 1);
    for (const s of SCHEMES) {
      if (!st.show[s]) continue;
      const c = s0.c[s];
      const p = XS.map((x, i) => `${ax(x)},${ay(clampV(c[i]))}`).join(" ");
      svgEl("polyline", { points: p, fill: "none", stroke: COL[s], "stroke-width": s === "ftcs" ? 2 : 2.5 }, pg);
      XS.forEach((x, i) => {
        if (s === "ftcs" && i % 1) return;
        svgEl("circle", { cx: ax(x), cy: ay(clampV(c[i])), r: s === "ftcs" ? 3 : 3.2, fill: COL[s] }, pg);
      });
    }
    // legend
    const lg = (txt: string, i: number, c: string, op = 1) => {
      const x = WA - 210, y = top + 24 + 24 * i;
      svgEl("line", { x1: x, x2: x + 22, y1: y - 6, y2: y - 6, stroke: c, "stroke-width": 3, "stroke-opacity": op }, svgA);
      const t = svgEl("text", { x: x + 28, y, "font-weight": 600 }, svgA);
      t.style.fill = c; t.style.fillOpacity = String(op); t.textContent = txt;
    };
    lg("exact", 0, COL.exact, 0.6);
    if (st.show.ftcs && (s0.blownAt.ftcs >= 0 || Math.abs(s0.c.ftcs[N >> 1]) > 2 || maxAbs(s0.c.ftcs) > 3)) {
      const t = svgEl("text", { x: la + 10, y: bot - 12, "font-weight": 700 }, svgA);
      t.style.fill = COL.warn;
      t.textContent = s0.blownAt.ftcs >= 0 ? `FTCS overflowed after ${s0.blownAt.ftcs} steps` : `FTCS: max |c| = ${sci(maxAbs(s0.c.ftcs), 2)}`;
    }

    // ---------- (b) error vs dt ----------
    const ec = errCurve(ic);
    svgB.replaceChildren();
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgB).textContent = `(b) max error at t = ${TREF[ic]} s against Δt`;
    for (let k = LE0; k <= 0; k++) {
      svgEl("line", { x1: 72, x2: WB - 14, y1: by(k), y2: by(k), class: "grid" }, svgB);
      svgEl("text", { x: 66, y: by(k) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = pow10(k);
    }
    svgEl("line", { x1: 72, x2: WB - 14, y1: bot, y2: bot, class: "axis" }, svgB);
    const xt: [number, string][] = [[-2, "0.01"], [-1, "0.1"], [0, "1"], [1, "10"]];
    for (const [k, lab] of xt) {
      svgEl("line", { x1: bx(k), x2: bx(k), y1: top, y2: bot + 5, class: "grid" }, svgB);
      svgEl("text", { x: bx(k), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgB).textContent = lab;
    }
    svgEl("text", { x: (72 + WB) / 2, y: H - 6, "text-anchor": "middle" }, svgB).textContent = "time step Δt (s)";
    const yl2 = svgEl("text", { x: 18, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 18 ${(top + bot) / 2})` }, svgB);
    yl2.textContent = "max |error|";
    const cpb = svgEl("clipPath", { id: `im-cb-${uid}` }, svgB);
    svgEl("rect", { x: 72, y: top - 2, width: WB - 86, height: bot - top + 2 }, cpb);
    const pgb = svgEl("g", { "clip-path": `url(#im-cb-${uid})` }, svgB);
    const ly = (e: number) => Number.isFinite(e) && e > 0 ? Math.log10(e) : LE1 + 1;
    // FTCS limit
    const lim = Math.log10(dtOf(0.5));
    svgEl("rect", { x: bx(lim), y: top, width: WB - 14 - bx(lim), height: bot - top, fill: COL.ftcs, "fill-opacity": 0.07 }, pgb);
    svgEl("line", { x1: bx(lim), x2: bx(lim), y1: top, y2: bot, stroke: COL.ftcs, "stroke-width": 2, "stroke-dasharray": "6 4" }, pgb);
    // slope guides, anchored below BE and CN at dt = 1 s
    const ld = ec.dt.map(Math.log10);
    const at1 = (s: Scheme) => { let j = 0; for (let i = 0; i < ld.length; i++) if (Math.abs(ld[i]) < Math.abs(ld[j])) j = i; return Math.log10(ec.err[s][j]) - ld[j] * (s === "be" ? 1 : 2); };
    const guide = (p: number, c0: number, xa: number, xb: number, lab: string) => {
      svgEl("line", { x1: bx(xa), x2: bx(xb), y1: by(c0 + p * xa), y2: by(c0 + p * xb), stroke: "#777", "stroke-width": 2, "stroke-dasharray": "8 6" }, pgb);
      svgEl("text", { x: bx(xa) + 4, y: by(c0 + p * xa) + 22, class: "sm" }, pgb).textContent = lab;
    };
    guide(1, at1("be") - 0.5, -1.6, 0.2, "slope 1");
    guide(2, at1("cn") - 0.6, -0.8, 0.9, "slope 2");
    for (const s of SCHEMES) {
      if (!st.show[s]) continue;
      const segs: string[][] = [[]];
      ec.err[s].forEach((e, i) => {
        const y = ly(e);
        if (y > LE1 + 0.3) { segs.push([]); return; }
        segs[segs.length - 1].push(`${bx(ld[i])},${by(y)}`);
      });
      for (const sg of segs) if (sg.length > 1) svgEl("polyline", { points: sg.join(" "), fill: "none", stroke: COL[s], "stroke-width": 2.5 }, pgb);
      ec.err[s].forEach((e, i) => {
        const y = ly(e);
        if (y <= LE1 + 0.3) svgEl("circle", { cx: bx(ld[i]), cy: by(y), r: 2.6, fill: COL[s] }, pgb);
      });
    }
    const tl = svgEl("text", { x: bx(lim) + 6, y: top + 20, "font-weight": 600 }, svgB);
    tl.style.fill = "#a06800"; tl.textContent = "r > ½: FTCS unstable";
    // current dt markers
    const [nRef, dtRef] = refSteps(ic, dt);
    svgEl("line", { x1: bx(Math.log10(dtRef)), x2: bx(Math.log10(dtRef)), y1: top, y2: bot, stroke: COL.ink, "stroke-width": 1.2, "stroke-dasharray": "2 3" }, svgB);
    const curErr = {} as Record<Scheme, number>;
    for (const s of SCHEMES) {
      curErr[s] = errAt(ic, s, nRef, dtRef);
      if (!st.show[s]) continue;
      const y = ly(curErr[s]);
      if (y > LE1 + 0.3) continue;
      svgEl("circle", { cx: bx(Math.log10(dtRef)), cy: by(Math.max(LE0, y)), r: 8, fill: COL[s], stroke: COL.ink, "stroke-width": 2 }, svgB);
    }

    // ---------- readout ----------
    ro.replaceChildren();
    const item = (k: string, v: string, cls = "") => {
      const sp = htmlEl("span", cls, "", ro);
      htmlEl("b", "", k + " ", sp); sp.append(v);
    };
    item("r =", `${fmtR(r)}, Δt = ${fmtDt(dt)}`);
    item("zigzag factor per step:", `FTCS ${sci(1 - 4 * r, 3)}, BE ${sci(1 / (1 + 4 * r), 3)}, CN ${sci((1 - 2 * r) / (1 + 2 * r), 3)}`);
    const errNow = (s: Scheme) => s0.blownAt[s] >= 0 ? Infinity : maxErr(ic, s0.c[s], tn);
    item(`error at t = ${tn.toFixed(1)} s:`, SCHEMES.map(s => `${s === "ftcs" ? "FTCS" : s.toUpperCase()} ${sci(errNow(s), 2)}`).join(", "),
      Number.isFinite(errNow("ftcs")) && errNow("ftcs") < 10 ? "" : "im-warn");
    item(`error at t = ${TREF[ic]} s (Δt = ${fmtDt(dtRef)}):`, `BE ${sci(curErr.be, 3)}, CN ${sci(curErr.cn, 3)}`);
    if (ic === "step") item("c(10 µm):", `BE ${s0.c.be[1].toFixed(4)}, CN ${s0.c.cn[1].toFixed(4)}, exact ${exact(ic, 10, tn).toFixed(4)}`);
    else item("amplitude (x = L/2):", `BE ${s0.c.be[N / 2].toFixed(4)}, CN ${s0.c.cn[N / 2].toFixed(4)}, exact ${exact(ic, L / 2, tn).toFixed(4)}`);
    const ftcsSteps = Math.ceil(tn / dtOf(0.5) - 1e-9);
    item("cost:", `${fmtInt(n)} steps = ${fmtInt(n)} Thomas solves of ${N - 1} unknowns ` +
      `(BE ≈ ${fmtInt(n * OPS.be * (N - 1))}, CN ≈ ${fmtInt(n * OPS.cn * (N - 1))} operations); ` +
      `stable FTCS (r = ½) needs ${fmtInt(ftcsSteps)} steps ≈ ${fmtInt(ftcsSteps * OPS.ftcs * (N - 1))} operations`);
  }
  draw();
  return stop;
}

function maxAbs(c: Float64Array): number {
  let m = 0;
  for (const v of c) { const a = Math.abs(v); if (!(a <= m)) m = a; }
  return m;
}

// ---------------------------------------------------------------- embryo tab
const DT_W = 1.4e5, LE = 500, NE = 200, DXE = LE / NE;
const XE = Array.from({ length: NE + 1 }, (_, i) => i * DXE);
const TOUT = XE.map(x => x < LE / 2 ? 20 : x > LE / 2 ? 27 : 23.5);
const IA = Math.round(125 / DXE), IP = Math.round(375 / DXE);
type Mode = "on" | "off";
const TMAX: Record<Mode, number> = { on: 0.1, off: 1 };

function embStep(u: Float64Array, dt: number, theta: number, kex: number): Float64Array {
  const me = NE + 1, r = DT_W * dt / (DXE * DXE), k = kex * dt;
  const a = new Float64Array(me), b = new Float64Array(me), c = new Float64Array(me), d = new Float64Array(me);
  for (let i = 0; i < me; i++) {
    const lo = i === 0 ? 0 : i === me - 1 ? 2 : 1, up = i === 0 ? 2 : i === me - 1 ? 0 : 1;
    const Ku = -2 * u[i] + lo * (i > 0 ? u[i - 1] : 0) + up * (i < me - 1 ? u[i + 1] : 0);
    d[i] = u[i] + (1 - theta) * (r * Ku + k * (TOUT[i] - u[i])) + theta * k * TOUT[i];
    a[i] = -theta * r * lo; b[i] = 1 + 2 * theta * r + theta * k; c[i] = -theta * r * up;
  }
  return thomas(a, b, c, d);
}
function stopflowExact(x: number, t: number): number {
  if (t <= 0) return x < LE / 2 ? 20 : x > LE / 2 ? 27 : 23.5;
  let s = 23.5;
  for (let k = 1; k <= 8001; k += 2) {
    const e = Math.exp(-DT_W * (k * Math.PI / LE) ** 2 * t);
    if (e < 1e-16) break;
    s -= 14 / (k * Math.PI) * Math.sin(k * Math.PI / 2) * Math.cos(k * Math.PI * x / LE) * e;
  }
  return s;
}
/** 10-90 % crossing positions (20.7 and 26.3 °C), or NaN. */
function width(T: ArrayLike<number>): [number, number] {
  let x10 = NaN, x90 = NaN;
  for (let i = 0; i < NE; i++) if (T[i] < 20.7 && T[i + 1] >= 20.7) { x10 = XE[i] + DXE * (20.7 - T[i]) / (T[i + 1] - T[i]); break; }
  for (let i = NE; i > 0; i--) if (T[i] > 26.3 && T[i - 1] <= 26.3) { x90 = XE[i - 1] + DXE * (26.3 - T[i - 1]) / (T[i] - T[i - 1]); break; }
  return [x10, x90];
}
const tempCol = (T: number): string => {
  const f = Math.max(0, Math.min(1, (T - 20) / 7));
  const c0 = [0x00, 0x72, 0xb2], c1 = [0xd5, 0x5e, 0x00];
  return `rgb(${c0.map((v, i) => Math.round(v + f * (c1[i] - v))).join(",")})`;
};

function buildEmbryo(box: HTMLElement, uid: number): () => void {
  const DTS = [0.1e-3, 1e-3, 5e-3, 50e-3];
  const START = { mode: "on" as Mode, sch: "cn" as "be" | "cn", dt: 1e-3, R: 100, t: 0.1 };
  const st = { ...START };
  const tau = () => st.R ** 2 / (2.405 ** 2 * DT_W);
  const lam = () => Math.sqrt(DT_W * tau());

  // full runs: profiles at every step (at most 1000 steps x 201)
  let run: { key: string; dt: number; prof: Float64Array[] } | null = null;
  function getRun() {
    const key = `${st.mode}|${st.sch}|${st.dt}|${st.R}`;
    if (run && run.key === key) return run;
    const n = Math.round(TMAX[st.mode] / st.dt);
    const stride = Math.max(1, Math.ceil(n / 1000));
    const prof: Float64Array[] = [];
    let u: Float64Array = st.mode === "on" ? new Float64Array(NE + 1).fill(20) : Float64Array.from(TOUT);
    const th = st.sch === "be" ? 1 : 0.5, kex = st.mode === "on" ? 1 / tau() : 0;
    prof.push(u);
    for (let k = 1; k <= n; k++) { u = embStep(u, st.dt, th, kex); if (k % stride === 0) prof.push(u); }
    run = { key, dt: st.dt * stride, prof };
    return run;
  }

  // ---- controls ----
  const row1 = htmlEl("div", "im-row", "", box);
  const g1 = htmlEl("div", "im-group", "", row1);
  const modeBtns = presets(g1, [
    { label: "streams on (from 20 °C)", apply: () => setMode("on") },
    { label: "flow stopped (from the step)", apply: () => setMode("off") },
  ]);
  const g2 = htmlEl("div", "im-group", "", row1);
  const schBtns = presets(g2, [
    { label: "backward Euler", apply: () => { st.sch = "be"; draw(); } },
    { label: "Crank–Nicolson", apply: () => { st.sch = "cn"; draw(); } },
  ]);
  const row2 = htmlEl("div", "im-row", "", box);
  const g3 = htmlEl("div", "im-group", "", row2);
  htmlEl("span", "", "Δt =", g3);
  const dtBtns = presets(g3, DTS.map(d => ({ label: fmtDt(d), apply: () => { st.dt = d; draw(); } })));
  const RS = slider(row2, "radius R", 50, 150, 5, st.R, v => { st.R = v; draw(); }, v => `${v} µm`);
  const row3 = htmlEl("div", "im-row", "", box);
  const tS = slider(row3, "time t", 0, 1, 0.005, st.t / TMAX[st.mode], v => { stop(); st.t = v * TMAX[st.mode]; draw(); },
    v => fmtDt(Math.max(v * TMAX[st.mode], 0)));
  const g4 = htmlEl("div", "im-group", "", row3);
  const playB = button(g4, "▶ play", () => { if (timer) stop(); else play(); }, "im-play");
  button(g4, "reset", () => {
    stop(); Object.assign(st, START); RS.set(st.R); tS.set(st.t / TMAX[st.mode]); draw();
  });

  function setMode(m: Mode): void {
    stop(); st.mode = m; st.t = m === "on" ? 0.1 : 0.05;
    if (m === "off" && st.dt < 1e-3) st.dt = 1e-3;
    tS.set(st.t / TMAX[m]); draw();
  }

  let timer = 0, last = 0;
  function play(): void {
    if (st.t >= TMAX[st.mode] * 0.999) st.t = 0;
    playB.textContent = "❚❚ pause";
    last = performance.now();
    const tick = (now: number) => {
      st.t = Math.min(TMAX[st.mode], st.t + TMAX[st.mode] / 8 * (now - last) / 1000); last = now;
      tS.set(Math.round(st.t / TMAX[st.mode] * 200) / 200);
      draw();
      if (st.t >= TMAX[st.mode]) { stop(); return; }
      timer = requestAnimationFrame(tick);
    };
    timer = requestAnimationFrame(tick);
  }
  function stop(): void {
    if (timer) cancelAnimationFrame(timer);
    timer = 0; playB.textContent = "▶ play";
  }

  const plots = htmlEl("div", "im-plots", "", box);
  const H = 370, WA = 620, WB = 480, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Temperature along the embryo" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 5px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Temperature against time" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 5px)`;
  const ro = htmlEl("div", "qb-readout im-read", "", box);

  const la = 62, top = 118, bot = H - 50;
  const ax = scale(0, LE, la, WA - 14), ay = scale(19.5, 27.5, bot, top);

  function draw(): void {
    modeBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.mode === "on")));
    schBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.sch === "be")));
    dtBtns.forEach((b, i) => b.classList.toggle("active", DTS[i] === st.dt));
    const rn = getRun();
    const k = Math.min(rn.prof.length - 1, Math.floor(st.t / rn.dt + 1e-9));
    const tn = k * rn.dt, T = rn.prof[k];
    const tu = tau(), lm = lam(), on = st.mode === "on";
    const schCol = st.sch === "be" ? COL.be : COL.cn;

    // ---------- (c) ----------
    svgA.replaceChildren();
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgA).textContent =
      `(c) temperature along the embryo, t = ${fmtDt(tn)}`;
    // streams and embryo
    const sy = 34, sh = 66, mid = ax(LE / 2);
    svgEl("rect", { x: la, y: sy, width: mid - la, height: sh, fill: COL.be, "fill-opacity": on ? 0.22 : 0.06 }, svgA);
    svgEl("rect", { x: mid, y: sy, width: WA - 14 - mid, height: sh, fill: COL.cn, "fill-opacity": on ? 0.22 : 0.06 }, svgA);
    const tl = svgEl("text", { x: la + 6, y: sy + 20, "font-weight": 600 }, svgA); tl.style.fill = COL.be;
    tl.textContent = on ? "20 °C stream" : "no flow";
    const tr = svgEl("text", { x: WA - 20, y: sy + 20, "font-weight": 600, "text-anchor": "end" }, svgA); tr.style.fill = COL.cn;
    tr.textContent = on ? "27 °C stream" : "no flow";
    const gid = `im-grad-${uid}`;
    const defs = svgEl("defs", {}, svgA);
    const lgd = svgEl("linearGradient", { id: gid, x1: 0, x2: 1, y1: 0, y2: 0 }, defs);
    for (let i = 0; i <= NE; i += 5) svgEl("stop", { offset: (i / NE).toFixed(3), "stop-color": tempCol(T[i]) }, lgd);
    svgEl("ellipse", { cx: (ax(0) + ax(LE)) / 2, cy: sy + sh / 2 + 6, rx: (ax(LE) - ax(0)) / 2, ry: 20, fill: `url(#${gid})`, stroke: COL.ink, "stroke-width": 2 }, svgA);
    // axes
    for (let v = 20; v <= 27; v += 1) {
      svgEl("line", { x1: la, x2: WA - 14, y1: ay(v), y2: ay(v), class: "grid" }, svgA);
      svgEl("text", { x: la - 6, y: ay(v) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = String(v);
    }
    for (let x = 0; x <= LE; x += 100) {
      svgEl("line", { x1: ax(x), x2: ax(x), y1: bot, y2: bot + 5, class: "axis" }, svgA);
      svgEl("text", { x: ax(x), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgA).textContent = String(x);
    }
    svgEl("line", { x1: la, x2: WA - 14, y1: bot, y2: bot, class: "axis" }, svgA);
    svgEl("text", { x: (la + WA) / 2, y: H - 6, "text-anchor": "middle" }, svgA).textContent = "position along the embryo x (µm)";
    const yl = svgEl("text", { x: 18, y: (top + bot) / 2, "text-anchor": "middle", transform: `rotate(-90 18 ${(top + bot) / 2})` }, svgA);
    yl.textContent = "T (°C)";
    // band
    const [x10, x90] = width(T);
    if (Number.isFinite(x10) && Number.isFinite(x90) && x90 > x10)
      svgEl("rect", { x: ax(x10), y: top, width: ax(x90) - ax(x10), height: bot - top, fill: COL.yellow, "fill-opacity": 0.45 }, svgA);
    // stream (dashed step)
    if (on) svgEl("polyline", { points: `${ax(0)},${ay(20)} ${mid},${ay(20)} ${mid},${ay(27)} ${ax(LE)},${ay(27)}`, fill: "none", stroke: COL.grey, "stroke-width": 2.5, "stroke-dasharray": "7 5" }, svgA);
    // reference: steady formula (on) or exact stop-flow series (off)
    for (let i = 0; i <= NE; i += 8) {
      const ref = on ? TOUT[i] - 3.5 * Math.sign(XE[i] - LE / 2) * Math.exp(-Math.abs(XE[i] - LE / 2) / lm) : stopflowExact(XE[i], tn);
      svgEl("circle", { cx: ax(XE[i]), cy: ay(ref), r: 4, fill: "none", stroke: COL.ink, "stroke-width": 1.5 }, svgA);
    }
    const cpa = svgEl("clipPath", { id: `im-ce-${uid}` }, svgA);
    svgEl("rect", { x: la, y: top - 4, width: WA - 14 - la, height: bot - top + 8 }, cpa);
    const pga = svgEl("g", { "clip-path": `url(#im-ce-${uid})` }, svgA);
    svgEl("polyline", { points: XE.map((x, i) => `${ax(x)},${ay(T[i])}`).join(" "), fill: "none", stroke: schCol, "stroke-width": 3.5 }, pga);
    svgEl("circle", { cx: ax(XE[IA]), cy: ay(T[IA]), r: 7, fill: COL.sky, stroke: COL.ink, "stroke-width": 1.5 }, svgA);
    svgEl("circle", { cx: ax(XE[IP]), cy: ay(T[IP]), r: 7, fill: COL.purple, stroke: COL.ink, "stroke-width": 1.5 }, svgA);
    // legend
    const lg = (txt: string, i: number, c: string, kind: "line" | "dash" | "circ") => {
      const x = la + 10, y = top + 22 + 22 * i;
      if (kind === "circ") svgEl("circle", { cx: x + 11, cy: y - 6, r: 4, fill: "none", stroke: c, "stroke-width": 1.5 }, svgA);
      else svgEl("line", { x1: x, x2: x + 22, y1: y - 6, y2: y - 6, stroke: c, "stroke-width": 3, "stroke-dasharray": kind === "dash" ? "6 4" : "" }, svgA);
      const t = svgEl("text", { x: x + 28, y, class: "sm" }, svgA);
      t.textContent = txt;
    };
    lg(st.sch === "be" ? "backward Euler" : "Crank–Nicolson", 0, schCol, "line");
    lg(on ? "steady formula" : "exact (series)", 1, COL.ink, "circ");
    if (on) lg("stream outside", 2, COL.grey, "dash");

    // ---------- (d) ----------
    svgB.replaceChildren();
    const tmax = TMAX[st.mode], ms = on;
    const bl = 62, bx = scale(0, tmax, bl, WB - 14), by = scale(19.5, 27.5, bot, 40);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgB).textContent = "(d) temperature against time";
    for (let v = 20; v <= 27; v += 1) {
      svgEl("line", { x1: bl, x2: WB - 14, y1: by(v), y2: by(v), class: "grid" }, svgB);
      svgEl("text", { x: bl - 6, y: by(v) + 6, "text-anchor": "end", class: "sm" }, svgB).textContent = String(v);
    }
    const tstep = tmax / 5;
    for (let i = 0; i <= 5; i++) {
      const t = i * tstep;
      svgEl("line", { x1: bx(t), x2: bx(t), y1: bot, y2: bot + 5, class: "axis" }, svgB);
      svgEl("text", { x: bx(t), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgB).textContent = ms ? String(Math.round(t * 1000)) : t.toFixed(1);
    }
    svgEl("line", { x1: bl, x2: WB - 14, y1: bot, y2: bot, class: "axis" }, svgB);
    svgEl("text", { x: (bl + WB) / 2, y: H - 6, "text-anchor": "middle" }, svgB).textContent = ms ? "time after switch-on (ms)" : "time after the flow stops (s)";
    if (on) for (const [m, lab] of [[1, "τ"], [4, "4τ"]] as [number, string][]) {
      if (m * tu > tmax) continue;
      svgEl("line", { x1: bx(m * tu), x2: bx(m * tu), y1: 40, y2: bot, stroke: COL.green, "stroke-width": 2, "stroke-dasharray": "5 4" }, svgB);
      const t = svgEl("text", { x: bx(m * tu) + 4, y: by(23.5) + 6, "font-weight": 600 }, svgB);
      t.style.fill = COL.green; t.textContent = lab;
    }
    for (const [i, c] of [[IA, COL.sky], [IP, COL.purple]] as [number, string][]) {
      const p = rn.prof.map((u, j) => `${bx(j * rn.dt)},${by(Math.max(19, Math.min(28, u[i])))}`).join(" ");
      svgEl("polyline", { points: p, fill: "none", stroke: c, "stroke-width": 3.5 }, svgB);
      svgEl("circle", { cx: bx(tn), cy: by(T[i]), r: 7, fill: c, stroke: COL.ink, "stroke-width": 1.5 }, svgB);
    }
    svgEl("line", { x1: bx(tn), x2: bx(tn), y1: 40, y2: bot, stroke: COL.ink, "stroke-width": 1.2, "stroke-dasharray": "2 3" }, svgB);
    const lgB = (txt: string, y: number, c: string) => {
      svgEl("circle", { cx: WB - 210, cy: y - 6, r: 6, fill: c, stroke: COL.ink, "stroke-width": 1.2 }, svgB);
      svgEl("text", { x: WB - 198, y, class: "sm" }, svgB).textContent = txt;
    };
    lgB("posterior, 375 µm", on ? by(24.6) : by(26.2), COL.purple);
    lgB("anterior, 125 µm", on ? by(23.9) : by(25.5), COL.sky);

    // ---------- readout ----------
    ro.replaceChildren();
    const item = (k: string, v: string, cls = "") => {
      const sp = htmlEl("span", cls, "", ro);
      htmlEl("b", "", k + " ", sp); sp.append(v);
    };
    const r = DT_W * st.dt / (DXE * DXE);
    item("r =", `${sci(r, 3)}`);
    item("steps:", `${fmtInt(tn / st.dt)} Thomas solves of ${NE + 1} unknowns`);
    if (on) {
      item("τ =", `${(tu * 1000).toFixed(1)} ms, λ = √(Dτ) = ${lm.toFixed(1)} µm`);
      item("10–90 % width:", `${Number.isFinite(x10) && Number.isFinite(x90) ? (x90 - x10).toFixed(1) + " µm" : "–"} (formula 2λ ln 5 = ${(2 * lm * Math.log(5)).toFixed(1)} µm)`);
      let dev = 0;
      for (let i = 0; i <= NE; i++) dev = Math.max(dev, Math.abs(T[i] - (TOUT[i] - 3.5 * Math.sign(XE[i] - LE / 2) * Math.exp(-Math.abs(XE[i] - LE / 2) / lm))));
      item("max |T − steady formula|:", `${sci(dev, 2)} °C`);
    } else {
      let dev = 0;
      for (let i = 0; i <= NE; i++) dev = Math.max(dev, Math.abs(T[i] - stopflowExact(XE[i], tn)));
      item("10–90 % width:", Number.isFinite(x10) && Number.isFinite(x90) ? `${(x90 - x10).toFixed(1)} µm` : "–");
      item("max |T − exact|:", `${sci(dev, 2)} °C`);
      let flips = 0;
      if (k >= 2) { const a = rn.prof[k], b = rn.prof[k - 1], c = rn.prof[k - 2]; for (let i = 0; i <= NE; i++) if ((a[i] - b[i]) * (b[i] - c[i]) < -1e-12) flips++; }
      if (flips > 0) item("", `${flips} points change direction between steps`, "im-warn");
    }
    item("T(200 µm) =", `${T[80].toFixed(4)} °C`);
  }
  draw();
  return stop;
}

autoMount(NAME, mount);

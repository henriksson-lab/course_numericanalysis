// Gillespie lab for the gillespie unit ("Stochastic simulation").
// Two genes with the same ODE mean of 10 mRNA (time in hours):
//   constitutive: make at k = 10, degrade at delta m, delta = 1;
//   bursty (telegraph): promoter off -> on at k_on, on -> off at k_off = 9 k_on,
//   make at k_m = 100 while on, degrade at delta m. The promoter is on 10 % of the time
//   for every k_on, so the mean stays 10 while the burst size k_m / k_off changes.
// Tab "One cell": a single trajectory as a step function (promoter-on periods shaded),
//   the ODE mean overlaid. Step = one reaction: (b) shows the exponential waiting time
//   drawn by inverse transform, and the stacked bar of propensities with u2 a0 marked.
//   Run animates, Reset restarts the same cell, "New cell" uses the next seed.
//   With tau-leaping (constitutive only) each step is a leap of tau hours with Poisson
//   numbers of reactions; negative counts are marked, then clipped to 0.
// Tab "Many cells": N cells to time T; (c) histogram against the exact distribution
//   (Poisson, or the telegraph master equation); (d) Fano factor against k_on.
// Auto-mounts on <div data-demo="ssa">. Optional attributes:
//   data-tab="cell|pop"  data-model="const|tele"  data-kon="0.5"  data-tau="0|1.0"
//   data-n="1000"  data-t="30"  data-tend="20"  data-seed="1"
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, rng, autoMount } from "../lib/ui";
import * as M from "./math";

const NAME = "ssa";
// Okabe-Ito, as in make-gillespie-figs.R (constitutive blue, bursty vermillion, promoter on green)
const COL = {
  const: "#0072B2", tele: "#D55E00", on: "#009E73", off: "#CC79A7", make: "#56B4E9", deg: "#E69F00",
  ink: "#222", ref: "#737373", neg: "#b00020",
};
const RCOL: Record<M.Kind, string[]> = { const: [COL.make, COL.deg], tele: [COL.on, COL.off, COL.make, COL.deg] };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .sa-row { display: flex; flex-wrap: wrap; gap: .35em 1.2em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .sa-group { display: flex; gap: .35em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .sa-row .qb-slider { min-width: 200px; flex: 0 1 280px; grid-template-columns: auto 1fr 3em; }
[data-demo="${NAME}"] .sa-nowrap { flex-wrap: nowrap; white-space: nowrap; }
[data-demo="${NAME}"] .sa-nowrap .qb-slider { flex: 0 0 230px; min-width: 0; }
[data-demo="${NAME}"] .sa-nowrap > * { flex-shrink: 0; }
[data-demo="${NAME}"] .sa-plots { display: flex; flex-wrap: wrap; gap: 10px; width: 100%; }
[data-demo="${NAME}"] .sa-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .sa-read { display: flex; flex-wrap: wrap; gap: .15em 1.3em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .sa-read span b { font-weight: 600; }
[data-demo="${NAME}"] .sa-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"] .off { opacity: .4; pointer-events: none; }
[data-demo="${NAME}"] .sa-tabs .qb-btn { font-size: .95em; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

const fx = (v: number, d: number): string => Number.isFinite(v) ? v.toFixed(d).replace(/^-/, "−") : "–";
const f3 = (v: number): string => !Number.isFinite(v) ? "–" : Math.abs(v) >= 100 ? v.toFixed(1) : Math.abs(v) >= 10 ? v.toFixed(2) : v.toFixed(3);
const fk = (v: number): string => v >= 10 ? v.toFixed(0) : v >= 1 ? v.toPrecision(2) : v.toFixed(2);
const pct = (v: number): string => (100 * v).toFixed(1) + " %";
const NVALS = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000];
const KON0 = 0.2, KON1 = 20;

interface Ev { t: number; m: number; g: number; raw: number; }
type Pending = { t: number; ssa: M.SsaStep; leap?: undefined } | { t: number; leap: M.Leap; ssa?: undefined };

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (s: string | undefined, d: number, lo: number, hi: number): number => {
    const v = Number(s);
    return s !== undefined && s !== "" && Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : d;
  };
  const start = {
    tab: ds.tab === "pop" ? "pop" : "cell",
    kind: (ds.model === "const" ? "const" : "tele") as M.Kind,
    kon: num(ds.kon, 0.5, KON0, KON1),
    tauOn: Number(ds.tau) > 0,
    tau: Number(ds.tau) > 0 ? num(ds.tau, 1, 0.05, 1.5) : 1,
    N: num(ds.n, 1000, 10, 5000),
    T: num(ds.t, 30, 0.5, 30),
  };
  const seed = Math.round(num(ds.seed, 1, 0, 1e9));
  const tEnd = num(ds.tend, 20, 2, 40);
  const st = { ...start, cell: 0, speed: 2 };
  const model = (): M.Model => M.makeModel(st.kind, st.kon);
  const leaping = (): boolean => st.tauOn && st.kind === "const";

  // ---------------------------------------------------------------- controls
  const row1 = htmlEl("div", "sa-row", "", host);
  const tabG = htmlEl("div", "sa-group sa-tabs", "", row1);
  const tabBtns = presets(tabG, [
    { label: "One cell, step by step", apply: () => { st.tab = "cell"; showTab(); } },
    { label: "Many cells: distribution", apply: () => { st.tab = "pop"; showTab(); } },
  ]);
  const modG = htmlEl("div", "sa-group", "", row1);
  htmlEl("span", "", "gene:", modG);
  const modBtns = presets(modG, [
    { label: "constitutive", apply: () => { st.kind = "const"; changed(); } },
    { label: "bursty (telegraph)", apply: () => { st.kind = "tele"; changed(); } },
  ]);
  button(row1, "reset all", () => {
    stopRun();
    Object.assign(st, start, { cell: 0 });
    konS.set(Math.log10(st.kon)); tauS.set(st.tau); tauBox.checked = st.tauOn;
    nS.set(NVALS.indexOf(nearestN(st.N))); tS.set(st.T);
    showTab(); changed();
  });

  const row2 = htmlEl("div", "sa-row", "", host);
  const konG = htmlEl("div", "sa-group sa-nowrap", "", row2);
  const konS = slider(konG, "k_on", Math.log10(KON0), Math.log10(KON1), 0.01, Math.log10(st.kon),
    v => { st.kon = 10 ** v; changed(); }, v => fk(10 ** v));
  const konInfo = htmlEl("span", "", "", konG);
  presets(konG, [0.5, 1, 5].map(k => ({ label: String(k), apply: () => { st.kon = k; konS.set(Math.log10(k)); changed(); } })));
  const tauG = htmlEl("div", "sa-group sa-nowrap", "", row2);
  const tauBox = checkbox(tauG, "tau-leaping", st.tauOn, v => { st.tauOn = v; changed(); });
  const tauS = slider(tauG, "τ (h)", 0.05, 1.5, 0.05, st.tau, v => { st.tau = v; changed(); });

  // ---- tab 1 ----
  const cellBox = htmlEl("div", "", "", host);
  const popBox = htmlEl("div", "", "", host);
  for (const b of [cellBox, popBox]) Object.assign(b.style, { display: "flex", flexDirection: "column", gap: ".4em", width: "100%" });
  const row3 = htmlEl("div", "sa-row", "", cellBox);
  const stepG = htmlEl("div", "sa-group", "", row3);
  const stepB = button(stepG, "Step", () => { stopRun(); stepOnce(); });
  const runB = button(stepG, "Run", () => { if (raf) stopRun(); else startRun(); });
  button(stepG, "Reset", () => { stopRun(); resetCell(); drawCell(); });
  button(stepG, "New cell", () => { stopRun(); st.cell++; resetCell(); drawCell(); });
  const spG = htmlEl("div", "sa-group", "", row3);
  htmlEl("span", "", "speed:", spG);
  const spBtns = presets(spG, [
    { label: "slow", apply: () => { st.speed = 0.5; } },
    { label: "medium", apply: () => { st.speed = 2; } },
    { label: "fast", apply: () => { st.speed = 8; } },
  ]);
  spBtns[1].classList.add("active");
  const plots1 = htmlEl("div", "sa-plots", "", cellBox);
  const H = 340, WA = 680, WB = 420, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "One cell: mRNA against time" }, plots1);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 6px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "The last step" }, plots1);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 6px)`;
  const ro1 = htmlEl("div", "qb-readout sa-read", "", cellBox);

  // ---- tab 2 ----
  const row4 = htmlEl("div", "sa-row", "", popBox);
  const nearestN = (n: number): number => NVALS.reduce((b, v) => Math.abs(Math.log(v / n)) < Math.abs(Math.log(b / n)) ? v : b, NVALS[0]);
  st.N = nearestN(st.N);
  const nS = slider(row4, "cells N", 0, NVALS.length - 1, 1, NVALS.indexOf(st.N), v => { st.N = NVALS[v]; schedulePop(); }, v => String(NVALS[v]));
  const tS = slider(row4, "time T (h)", 0.5, 30, 0.5, st.T, v => { st.T = v; schedulePop(); });
  const tG = htmlEl("div", "sa-group", "", row4);
  presets(tG, [1, 3, 30].map(T => ({ label: `T = ${T}`, apply: () => { st.T = T; tS.set(T); schedulePop(); } })));
  const plots2 = htmlEl("div", "sa-plots", "", popBox);
  const WC = 640, WD = 460, WT2 = WC + WD;
  const svgC = svgEl("svg", { viewBox: `0 0 ${WC} ${H}`, role: "img", "aria-label": "Histogram of mRNA counts" }, plots2);
  svgC.style.width = `calc(${(100 * WC / WT2).toFixed(2)}% - 6px)`;
  const svgD = svgEl("svg", { viewBox: `0 0 ${WD} ${H}`, role: "img", "aria-label": "Fano factor against k_on" }, plots2);
  svgD.style.width = `calc(${(100 * WD / WT2).toFixed(2)}% - 6px)`;
  const ro2 = htmlEl("div", "qb-readout sa-read", "", popBox);

  function showTab(): void {
    cellBox.style.display = st.tab === "cell" ? "flex" : "none";
    popBox.style.display = st.tab === "pop" ? "flex" : "none";
    tabBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.tab === "cell")));
    if (st.tab === "pop") schedulePop(); else stopRun();
  }

  function syncControls(): void {
    modBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.kind === "const")));
    konG.classList.toggle("off", st.kind === "const");
    tauG.classList.toggle("off", st.kind === "tele");
    tauS.input.parentElement!.classList.toggle("off", !st.tauOn);
    const koff = M.OFF_PER_ON * st.kon;
    konInfo.textContent = `k_off = ${fk(koff)}, burst ${fk(M.KM / koff)}`;
  }

  function changed(): void {
    stopRun();
    syncControls();
    resetCell();
    drawCell();
    schedulePop();
  }

  // ---------------------------------------------------------------- one cell
  let path: Ev[] = [], tNow = 0, pending: Pending | null = null, last: Pending | null = null;
  let R = rng(1), nNeg = 0, ode = M.odeMean(model(), tEnd);
  let raf = 0, lastTs = 0;

  function resetCell(): void {
    R = rng(seed * 1009 + st.cell + 8);
    path = [{ t: 0, m: 0, g: 0, raw: 0 }];
    tNow = 0; pending = null; last = null; nNeg = 0;
    ode = M.odeMean(model(), tEnd);
  }
  const cur = (): Ev => path[path.length - 1];
  function ensurePending(): Pending {
    if (pending) return pending;
    const c = cur();
    if (leaping()) {
      const leap = M.tauLeap(c.m, c.t, st.tau, R.unif);
      pending = { t: c.t + st.tau, leap };
    } else {
      const ssa = M.ssaStep(model(), { g: c.g, m: c.m }, c.t, R.unif);
      pending = { t: c.t + ssa.tau, ssa };
    }
    return pending;
  }
  function applyPending(): void {
    const p = pending!;
    if (p.leap) {
      path.push({ t: p.t, m: p.leap.m, g: 1, raw: p.leap.raw });
      if (p.leap.raw < 0) nNeg++;
    } else {
      const a = p.ssa.after;
      path.push({ t: p.t, m: a.m, g: a.g, raw: a.m });
    }
    last = p; pending = null; tNow = p.t;
  }
  function stepOnce(): void {
    if (tNow >= tEnd) return;
    const p = ensurePending();
    if (p.t > tEnd) tNow = tEnd; else applyPending();
    drawCell();
  }
  function startRun(): void {
    if (tNow >= tEnd) { resetCell(); }
    runB.textContent = "Pause";
    lastTs = 0;
    raf = window.requestAnimationFrame(tick);
  }
  function stopRun(): void {
    if (raf) window.cancelAnimationFrame(raf);
    raf = 0; runB.textContent = "Run";
  }
  function tick(ts: number): void {
    const dt = lastTs ? Math.min(0.1, (ts - lastTs) / 1000) : 0;
    lastTs = ts;
    const target = Math.min(tEnd, tNow + dt * st.speed);
    for (let guard = 0; guard < 20000; guard++) {
      const p = ensurePending();
      if (p.t > target) break;
      applyPending();
    }
    tNow = target;
    drawCell();
    if (tNow >= tEnd || !host.isConnected) { stopRun(); return; }
    raf = window.requestAnimationFrame(tick);
  }

  // panel geometry
  const L = 62, top = 40, bot = H - 50;

  function drawCell(): void {
    const md = model(), lp = leaping();
    stepB.disabled = tNow >= tEnd;
    // ---------- (a) trajectory ----------
    svgA.replaceChildren();
    let mx = 0, mn = 0;
    for (const e of path) { mx = Math.max(mx, e.m); mn = Math.min(mn, e.raw); }
    const ymax = Math.max(md.kind === "const" ? 25 : 40, Math.ceil((mx + 3) / 10) * 10);
    const ymin = mn < 0 ? Math.floor((mn - 1) / 5) * 5 : 0;
    const ax = scale(0, tEnd, L, WA - 14), ay = scale(ymin, ymax, bot, top);
    const tc = md.kind === "const" ? COL.const : COL.tele;
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgA).textContent =
      `(a) one cell, ${md.kind === "const" ? "constitutive" : "bursty"} gene${lp ? `, tau-leaping τ = ${st.tau} h` : " (SSA)"}`;
    // promoter-on bands
    if (md.kind === "tele") {
      for (let i = 0; i < path.length; i++) {
        if (!path[i].g) continue;
        const t1 = i + 1 < path.length ? path[i + 1].t : tNow;
        if (t1 <= path[i].t) continue;
        svgEl("rect", { x: ax(path[i].t), y: top, width: Math.max(1, ax(t1) - ax(path[i].t)), height: bot - top, fill: COL.on, "fill-opacity": 0.28 }, svgA);
      }
    }
    const ystep = ymax - ymin > 60 ? 20 : 10;
    for (let v = Math.ceil(ymin / ystep) * ystep; v <= ymax; v += ystep) {
      svgEl("line", { x1: L, x2: WA - 14, y1: ay(v), y2: ay(v), class: v === 0 ? "axis" : "grid" }, svgA);
      svgEl("text", { x: L - 6, y: ay(v) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = fx(v, 0);
    }
    if (ymin < 0) {
      svgEl("rect", { x: L, y: ay(0), width: WA - 14 - L, height: ay(ymin) - ay(0), fill: COL.neg, "fill-opacity": 0.08 }, svgA);
    }
    const tstep = tEnd > 20 ? 5 : 2;
    for (let t = 0; t <= tEnd + 1e-9; t += tstep) {
      svgEl("line", { x1: ax(t), x2: ax(t), y1: bot, y2: bot + 5, class: "axis" }, svgA);
      svgEl("text", { x: ax(t), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgA).textContent = String(t);
    }
    svgEl("line", { x1: L, x2: L, y1: top, y2: bot, class: "axis" }, svgA);
    svgEl("text", { x: (L + WA) / 2, y: H - 6, "text-anchor": "middle" }, svgA).textContent = "time t (h)";
    const cy = (top + bot) / 2;
    svgEl("text", { x: 16, y: cy, "text-anchor": "middle", transform: `rotate(-90 16 ${cy})` }, svgA).textContent = "mRNA molecules";
    // ODE mean
    svgEl("polyline", { points: ode.t.map((t, i) => `${ax(t).toFixed(1)},${ay(ode.m[i]).toFixed(1)}`).join(" "), fill: "none", stroke: COL.ink, "stroke-width": 3, "stroke-dasharray": "9 6" }, svgA);
    // trajectory
    const pts: string[] = [];
    for (let i = 0; i < path.length; i++) {
      const e = path[i];
      if (i > 0) pts.push(`${ax(e.t).toFixed(1)},${ay(path[i - 1].m).toFixed(1)}`);
      pts.push(`${ax(e.t).toFixed(1)},${ay(e.m).toFixed(1)}`);
    }
    pts.push(`${ax(tNow).toFixed(1)},${ay(cur().m).toFixed(1)}`);
    svgEl("polyline", { points: pts.join(" "), fill: "none", stroke: tc, "stroke-width": 2.4, "stroke-linejoin": "round" }, svgA);
    for (const e of path) {
      if (e.raw >= 0) continue;
      svgEl("line", { x1: ax(e.t), x2: ax(e.t), y1: ay(0), y2: ay(e.raw), stroke: COL.neg, "stroke-width": 2, "stroke-dasharray": "3 3" }, svgA);
      svgEl("circle", { cx: ax(e.t), cy: ay(e.raw), r: 6, fill: COL.neg }, svgA);
    }
    svgEl("circle", { cx: ax(tNow), cy: ay(cur().m), r: 5, fill: tc }, svgA);
    // legend
    const lgx = WA - 250;
    const leg = (y: number, draw: () => void, txt: string) => { draw(); svgEl("text", { x: lgx + 34, y }, svgA).textContent = txt; };
    leg(top + 20, () => svgEl("line", { x1: lgx, x2: lgx + 26, y1: top + 14, y2: top + 14, stroke: tc, "stroke-width": 3 }, svgA), "one cell");
    leg(top + 44, () => svgEl("line", { x1: lgx, x2: lgx + 26, y1: top + 38, y2: top + 38, stroke: COL.ink, "stroke-width": 3, "stroke-dasharray": "7 4" }, svgA), "ODE mean");
    if (md.kind === "tele") leg(top + 68, () => svgEl("rect", { x: lgx, y: top + 54, width: 26, height: 18, fill: COL.on, "fill-opacity": 0.35 }, svgA), "promoter on");
    if (lp && nNeg > 0) leg(top + 68, () => svgEl("circle", { cx: lgx + 13, cy: top + 62, r: 6, fill: COL.neg }, svgA), "negative count");

    drawStep();

    // ---------- readout ----------
    ro1.replaceChildren();
    const item = (k: string, v: string, cls = "") => { const s = htmlEl("span", cls, "", ro1); htmlEl("b", "", k + " ", s); s.append(v); };
    const c = cur();
    item("t =", fx(tNow, 3) + " h");
    item("m =", String(c.m));
    if (md.kind === "tele") item("promoter", c.g ? "on" : "off");
    item(lp ? "leaps:" : "reactions:", String(path.length - 1));
    if (!lp) {
      const a0 = M.propensities(md, { g: c.g, m: c.m }).reduce((x, y) => x + y, 0);
      item("now a₀ =", `${f3(a0)} /h, mean wait ${fx(60 / a0, 1)} min`);
    } else item("went negative:", `${nNeg} time${nNeg === 1 ? "" : "s"}`, nNeg ? "sa-warn" : "");
    if (tNow >= tEnd) item("", `end of the ${tEnd} h window: Reset or New cell`);
  }

  function drawStep(): void {
    svgB.replaceChildren();
    const md = model();
    const txt = (x: number, y: number, s: string, opt: Record<string, string | number> = {}, fill?: string) => {
      const t = svgEl("text", { x, y, ...opt }, svgB);
      if (fill) t.style.fill = fill;
      t.textContent = s;
      return t;
    };
    if (last && last.leap) {
      const lpp = last.leap;
      txt(8, 24, "(b) the last leap", { class: "ttl" });
      txt(14, 64, `from t = ${fx(lpp.t0, 2)} h to ${fx(lpp.t0 + lpp.tau, 2)} h,`);
      txt(14, 90, `propensities frozen at m = ${lpp.m0}:`);
      txt(14, 134, `make:  Poisson(kτ = ${fx(lpp.lamMake, 2)})`, {}, COL.make);
      txt(WB - 14, 134, `+${lpp.nMake}`, { "text-anchor": "end", "font-weight": 600 });
      txt(14, 166, `degrade:  Poisson(δmτ = ${fx(lpp.lamDeg, 2)})`, {}, COL.deg);
      txt(WB - 14, 166, `−${lpp.nDeg}`, { "text-anchor": "end", "font-weight": 600 });
      svgEl("line", { x1: 14, x2: WB - 14, y1: 182, y2: 182, class: "axis" }, svgB);
      txt(14, 210, `m: ${lpp.m0} + ${lpp.nMake} − ${lpp.nDeg} = ${fx(lpp.raw, 0)}`, { "font-weight": 600 }, lpp.raw < 0 ? COL.neg : undefined);
      if (lpp.raw < 0) {
        txt(14, 244, "a negative number of molecules!", { "font-weight": 600 }, COL.neg);
        txt(14, 270, "clipped to 0 (as in the lecture)", { class: "sm" });
      }
      txt(14, 316, `mean degradations δmτ = ${fx(lpp.lamDeg, 1)} vs m = ${lpp.m0}`, { class: "sm" });
      return;
    }
    txt(8, 24, "(b) the last reaction", { class: "ttl" });
    const s = last?.ssa;
    // propensities to show: of the last step, or of the current state
    const c = cur();
    const a = s ? s.a : M.propensities(md, { g: c.g, m: c.m });
    const a0 = a.reduce((x, y) => x + y, 0);
    const names = M.NAMES[md.kind], cols = RCOL[md.kind];
    // --- when? ---
    txt(14, 54, s ? `when?  a₀ = ${f3(a0)} /h,  u₁ = ${fx(s.u1, 3)}` : `when?  a₀ = ${f3(a0)} /h  (press Step)`, { "font-weight": 600 });
    const cx0 = 58, cx1 = WB - 18, cyT = 70, cyB = 160;
    const tauMax = Math.max(4 / a0, s ? s.tau * 1.15 : 0);
    const sx = scale(0, tauMax, cx0, cx1), sy = scale(0, 1, cyB, cyT);
    svgEl("line", { x1: cx0, x2: cx1, y1: cyB, y2: cyB, class: "axis" }, svgB);
    svgEl("line", { x1: cx0, x2: cx0, y1: cyT, y2: cyB, class: "axis" }, svgB);
    for (const v of [0, 1]) txt(cx0 - 6, sy(v) + 6, String(v), { "text-anchor": "end", class: "sm" });
    const cdf = Array.from({ length: 81 }, (_, i) => tauMax * i / 80).map(t => `${sx(t).toFixed(1)},${sy(1 - Math.exp(-a0 * t)).toFixed(1)}`).join(" ");
    svgEl("polyline", { points: cdf, fill: "none", stroke: COL.const, "stroke-width": 3 }, svgB);
    txt(cx1, cyB + 20, "waiting time τ", { "text-anchor": "end", class: "sm" });
    txt(cx1, sy(0.3), "P(wait ≤ τ) = 1 − e^(−a₀τ)", { "text-anchor": "end", class: "sm" });
    if (s) {
      const yv = sy(1 - s.u1);
      svgEl("line", { x1: cx0, x2: sx(s.tau), y1: yv, y2: yv, stroke: COL.tele, "stroke-width": 2.5 }, svgB);
      svgEl("line", { x1: sx(s.tau), x2: sx(s.tau), y1: yv, y2: cyB, stroke: COL.tele, "stroke-width": 2.5 }, svgB);
      svgEl("path", { d: `M${sx(s.tau) - 6},${cyB - 10} L${sx(s.tau) + 6},${cyB - 10} L${sx(s.tau)},${cyB} Z`, fill: COL.tele }, svgB);
      txt(14, 202, `τ = −ln(u₁)/a₀ = ${fx(s.tau, 4)} h = ${s.tau * 60 >= 1 ? fx(s.tau * 60, 1) + " min" : fx(s.tau * 3600, 0) + " s"}`, {}, COL.tele);
    } else txt(14, 202, `mean waiting time 1/a₀ = ${fx(60 / a0, 1)} min`, { class: "sm" });
    // --- which? ---
    txt(14, 236, s ? `which?  u₂ = ${fx(s.u2, 3)},  u₂·a₀ = ${f3(s.u2 * a0)}` : "which?  reaction j with probability aⱼ/a₀", { "font-weight": 600 });
    const bx0 = 14, bx1 = WB - 14, by0 = 250, bh = 34;
    const bs = scale(0, a0, bx0, bx1);
    let acc = 0;
    a.forEach((aj, j) => {
      const x0 = bs(acc), x1 = bs(acc + aj);
      acc += aj;
      if (aj > 0) {
        const chosen = s && s.j === j;
        svgEl("rect", { x: x0, y: by0, width: Math.max(0.5, x1 - x0), height: bh, fill: cols[j], stroke: chosen ? COL.ink : "#fff", "stroke-width": chosen ? 3.5 : 1 }, svgB);
        if (x1 - x0 > 62) txt((x0 + x1) / 2, by0 + bh / 2 + 6, names[j], { "text-anchor": "middle", "font-weight": 600 }, "#111");
      }
    });
    if (s) {
      const xm = bs(s.u2 * a0);
      svgEl("line", { x1: xm, x2: xm, y1: by0 - 6, y2: by0 + bh + 6, stroke: COL.ink, "stroke-width": 3 }, svgB);
      svgEl("path", { d: `M${xm - 7},${by0 - 12} L${xm + 7},${by0 - 12} L${xm},${by0 - 2} Z`, fill: COL.ink }, svgB);
    }
    // propensity values, coloured
    const tv = svgEl("text", { x: 14, y: by0 + bh + 24, class: "sm" }, svgB);
    a.forEach((aj, j) => {
      const sp = svgEl("tspan", {}, tv);
      sp.style.fill = cols[j]; sp.style.fontWeight = "600";
      sp.textContent = (j ? "  " : "") + `${names[j]} ${Number(aj.toPrecision(3))}`;
    });
    if (s) {
      const b = s.before, af = s.after;
      const what = md.kind === "tele" && s.j < 2 ? `promoter ${b.g ? "on" : "off"} → ${af.g ? "on" : "off"}` : `m ${b.m} → ${af.m}`;
      txt(14, H - 8, `→ ${names[s.j]}:  ${what}`, { "font-weight": 600 }, COL.ink);
    }
  }

  // ---------------------------------------------------------------- many cells
  let popTimer = 0;
  const pmfCache = new Map<string, Float64Array>();
  const curveCache = new Map<number, { k: number[]; F: number[] }>();
  const hist: { kon: number; F: number }[] = [];
  let histKey = "";
  function schedulePop(): void {
    if (st.tab !== "pop") return;
    if (popTimer) window.clearTimeout(popTimer);
    popTimer = window.setTimeout(() => { popTimer = 0; drawPop(); }, 60);
  }
  function getPmf(md: M.Model, T: number): Float64Array {
    const key = `${md.kind}|${md.kind === "tele" ? md.kon.toFixed(3) : ""}|${T}`;
    let p = pmfCache.get(key);
    if (!p) { p = M.pmf(md, T); if (pmfCache.size > 60) pmfCache.clear(); pmfCache.set(key, p); }
    return p;
  }
  function getCurve(T: number): { k: number[]; F: number[] } {
    let c = curveCache.get(T);
    if (!c) {
      const k: number[] = [], F: number[] = [];
      for (let i = 0; i <= 40; i++) {
        const kon = 10 ** (Math.log10(KON0) + i / 40 * (Math.log10(KON1) - Math.log10(KON0)));
        const mo = M.moments(M.makeModel("tele", kon), T);
        k.push(kon); F.push(mo.variance / mo.mean);
      }
      c = { k, F };
      curveCache.set(T, c);
    }
    return c;
  }

  function drawPop(): void {
    const md = model(), lp = leaping();
    // a seed fixed by the settings, so the same settings give the same cells
    const r = rng(seed * 7919 + st.N * 31 + Math.round(st.T * 4) * 7 + Math.round(st.kon * 1000) + (md.kind === "const" ? 13 : 0) + (lp ? Math.round(st.tau * 100) * 101 : 0));
    const pop = M.population(md, st.N, st.T, r.unif, lp ? st.tau : 0);
    const S = M.stats(pop.counts);
    const ex = M.moments(md, st.T);
    const p = getPmf(md, st.T);
    const tc = md.kind === "const" ? COL.const : COL.tele;

    // ---------- (c) histogram ----------
    svgC.replaceChildren();
    let q = 0, n99 = 0;
    for (let n = 0; n < p.length; n++) { q += p[n]; if (q > 0.995) { n99 = n; break; } }
    const xmax0 = Math.max(25, n99, S.max);
    const bw = xmax0 <= 60 ? 1 : xmax0 <= 120 ? 2 : 5;
    const xmax = Math.ceil((Math.min(xmax0, 250) + 1) / bw) * bw;
    const nb = xmax / bw;
    const cnt = new Array(nb).fill(0), ep = new Array(nb).fill(0);
    let over = 0;
    for (const v of pop.counts) { const b = Math.floor(v / bw); if (b < nb) cnt[b]++; else over++; }
    for (let n = 0; n < p.length; n++) { const b = Math.floor(n / bw); if (b < nb) ep[b] += p[n]; }
    const fr = cnt.map(c => c / st.N / bw), epd = ep.map(e => e / bw);
    const ymax = Math.max(...fr, ...epd) * 1.12 || 0.1;
    const Lc = 74, Rc = WC - 14;
    const sx = scale(0, xmax, Lc, Rc), sy = scale(0, ymax, bot, top);
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgC).textContent =
      `(c) mRNA in ${st.N} cells at T = ${st.T} h${lp ? `, tau-leap τ = ${st.tau}` : ""}`;
    const yt = niceStep(ymax / 4);
    for (let v = 0; v <= ymax; v += yt) {
      svgEl("line", { x1: Lc, x2: Rc, y1: sy(v), y2: sy(v), class: v === 0 ? "axis" : "grid" }, svgC);
      svgEl("text", { x: Lc - 6, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svgC).textContent = v === 0 ? "0" : v.toFixed(yt < 0.01 ? 3 : 2);
    }
    const xt = niceStep(xmax / 6);
    for (let v = 0; v <= xmax; v += xt) {
      svgEl("line", { x1: sx(v), x2: sx(v), y1: bot, y2: bot + 5, class: "axis" }, svgC);
      svgEl("text", { x: sx(v), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgC).textContent = String(v);
    }
    svgEl("text", { x: (Lc + Rc) / 2, y: H - 6, "text-anchor": "middle" }, svgC).textContent = "mRNA molecules per cell";
    const cy = (top + bot) / 2;
    svgEl("text", { x: 16, y: cy, "text-anchor": "middle", transform: `rotate(-90 16 ${cy})` }, svgC).textContent = "fraction of cells";
    for (let b = 0; b < nb; b++) {
      if (!cnt[b]) continue;
      svgEl("rect", { x: sx(b * bw) + 0.5, y: sy(fr[b]), width: Math.max(1, sx(bw) - sx(0) - 1), height: sy(0) - sy(fr[b]), fill: tc, "fill-opacity": 0.6 }, svgC);
    }
    for (let b = 0; b < nb; b++) {
      if (epd[b] < ymax * 0.004 && b * bw > ex.mean) continue;
      svgEl("circle", { cx: sx(b * bw + bw / 2), cy: sy(epd[b]), r: bw === 1 ? 3.6 : 3, fill: COL.ink }, svgC);
    }
    svgEl("line", { x1: sx(S.mean), x2: sx(S.mean), y1: top + 8, y2: bot, stroke: tc, "stroke-width": 2.5, "stroke-dasharray": "6 4" }, svgC);
    const lx = Rc - 230;
    svgEl("rect", { x: lx, y: top + 4, width: 24, height: 16, fill: tc, "fill-opacity": 0.6 }, svgC);
    svgEl("text", { x: lx + 32, y: top + 18 }, svgC).textContent = lp ? "tau-leaping" : "SSA";
    svgEl("circle", { cx: lx + 12, cy: top + 38, r: 4, fill: COL.ink }, svgC);
    svgEl("text", { x: lx + 32, y: top + 44 }, svgC).textContent = md.kind === "const" ? "exact: Poisson" : "exact (master eq.)";
    svgEl("line", { x1: lx, x2: lx + 24, y1: top + 64, y2: top + 64, stroke: tc, "stroke-width": 2.5, "stroke-dasharray": "6 4" }, svgC);
    svgEl("text", { x: lx + 32, y: top + 70 }, svgC).textContent = "sample mean";
    if (over) svgEl("text", { x: Rc, y: top + 96, "text-anchor": "end", class: "sm" }, svgC).textContent = `${over} cells above ${xmax}`;

    // ---------- (d) Fano against k_on ----------
    const key = `${st.N}|${st.T}`;
    if (key !== histKey) { hist.length = 0; histKey = key; }
    if (md.kind === "tele" && !hist.some(h => Math.abs(h.kon - st.kon) < 1e-9)) hist.push({ kon: st.kon, F: S.fano });
    svgD.replaceChildren();
    const curve = getCurve(st.T);
    const Ld = 56, Rd = WD - 92, xc = WD - 34;
    const dx = scale(Math.log10(KON0), Math.log10(KON1), Ld, Rd);
    const Fmax = Math.max(35, ...hist.map(h => h.F), S.fano) * 1.05;
    const dy = scale(0, Fmax, bot, top);
    svgEl("text", { x: 8, y: 24, class: "ttl" }, svgD).textContent = `(d) Fano factor at T = ${st.T} h`;
    for (let v = 0; v <= Fmax; v += 10) {
      svgEl("line", { x1: Ld, x2: WD - 10, y1: dy(v), y2: dy(v), class: v === 0 ? "axis" : "grid" }, svgD);
      svgEl("text", { x: Ld - 6, y: dy(v) + 6, "text-anchor": "end", class: "sm" }, svgD).textContent = String(v);
    }
    for (const k of [0.2, 0.5, 1, 2, 5, 10, 20]) {
      svgEl("line", { x1: dx(Math.log10(k)), x2: dx(Math.log10(k)), y1: bot, y2: bot + 5, class: "axis" }, svgD);
      svgEl("text", { x: dx(Math.log10(k)), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgD).textContent = String(k);
    }
    svgEl("text", { x: xc, y: bot + 22, "text-anchor": "middle", class: "sm" }, svgD).textContent = "const.";
    svgEl("text", { x: (Ld + Rd) / 2, y: H - 6, "text-anchor": "middle" }, svgD).textContent = "k_on (per h), k_off = 9 k_on";
    svgEl("text", { x: 14, y: cy, "text-anchor": "middle", transform: `rotate(-90 14 ${cy})` }, svgD).textContent = "variance / mean";
    svgEl("line", { x1: Ld, x2: WD - 10, y1: dy(1), y2: dy(1), stroke: COL.const, "stroke-width": 2, "stroke-dasharray": "6 4" }, svgD);
    svgEl("polyline", { points: curve.k.map((k, i) => `${dx(Math.log10(k)).toFixed(1)},${dy(curve.F[i]).toFixed(1)}`).join(" "), fill: "none", stroke: COL.tele, "stroke-width": 3 }, svgD);
    for (const h of hist) svgEl("circle", { cx: dx(Math.log10(h.kon)), cy: dy(h.F), r: 5, fill: "#fff", stroke: COL.tele, "stroke-width": 2.5 }, svgD);
    if (md.kind === "tele") svgEl("circle", { cx: dx(Math.log10(st.kon)), cy: dy(S.fano), r: 8, fill: COL.tele, stroke: COL.ink, "stroke-width": 2 }, svgD);
    else svgEl("circle", { cx: xc, cy: dy(S.fano), r: 8, fill: COL.const, stroke: COL.ink, "stroke-width": 2 }, svgD);
    const t1 = svgEl("text", { x: Rd, y: top + 22, "text-anchor": "end" }, svgD);
    t1.style.fill = COL.tele; t1.textContent = "telegraph, exact";
    const t2 = svgEl("text", { x: Rd, y: top + 46, "text-anchor": "end", class: "sm" }, svgD);
    t2.textContent = "dots: simulated";
    const t3 = svgEl("text", { x: Ld + 8, y: dy(1) - 8 }, svgD);
    t3.style.fill = COL.const; t3.textContent = "Poisson, F = 1";

    // ---------- readout ----------
    ro2.replaceChildren();
    const item = (k: string, v: string, cls = "") => { const s = htmlEl("span", cls, "", ro2); htmlEl("b", "", k + " ", s); s.append(v); };
    const pz = p[0];
    item("mean", `${f3(S.mean)} ± ${f3(S.se)} (exact ${f3(ex.mean)})`);
    item("variance", `${f3(S.variance)} (exact ${f3(ex.variance)})`);
    item("Fano", `${f3(S.fano)} (exact ${f3(ex.variance / ex.mean)})`);
    item("cells with 0 mRNA", `${pct(S.zero)} (exact ${pct(pz)})`);
    if (lp) {
      const tm = M.tauMoments(st.T, st.tau);
      item("tau-leap theory: variance", f3(tm.variance));
      item("runs that went negative", pct(pop.negFrac), pop.negFrac > 0 ? "sa-warn" : "");
    } else item("reactions simulated", pop.events.toLocaleString("en-US"));
  }

  // ---------------------------------------------------------------- start
  syncControls();
  resetCell();
  showTab();
  drawCell();
}

function niceStep(raw: number): number {
  const p = 10 ** Math.floor(Math.log10(raw));
  const m = raw / p;
  return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p;
}

autoMount(NAME, mount);

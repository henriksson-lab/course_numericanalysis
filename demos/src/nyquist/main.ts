// Sampling-and-spectrum lab for the fourier-fft unit.
// True signal: 8 + A cos(2 pi (t - 6)/P) + B cos(2 pi t / Pu) (a circadian clock gene with
// an optional fast "ultradian" rhythm), sampled every dt hours over a chosen duration,
// with seeded normal noise of sd s added to each sample.
//  (top)    the true signal (grey), the samples (blue) and the Fourier curve through the
//           samples (orange dashed: the trigonometric interpolant), or, with the low-pass
//           filter on, the smoothed curve built from the kept coefficients (green);
//  (bottom) the one-sided amplitude spectrum 2|X_k|/n (or 2|X_k|/sum(w) with the Hann
//           window), on a fixed axis in cycles per day, with the Nyquist frequency, the
//           true frequencies, and the fold of a fast rhythm onto its alias.
// With the defaults (dt = 2 h, 48 h, P = 24, B = 0, sd 0.3) the samples are exactly the
// 24 values of the lecture's circadian dataset: amplitude 1.4645, peak 6.45 h.
// The DFT is computed directly (n <= 240), O(n^2), in doubles.
// Auto-mounts on <div data-demo="nyquist">. Optional attributes:
//   data-dt="2" data-dur="48" data-period="24" data-fast-period="4.8" data-fast-amp="0"
//   data-noise="0.3" data-cut="8" data-window="0|1" data-filter="0|1" data-log="0|1"
import { svgEl, htmlEl, slider, presets, checkbox, scale, rng, autoMount, Slider } from "../lib/ui";

const NAME = "nyquist";
// Okabe-Ito, as in make-fourier-fft-figs.R
const OI = { orange: "#E69F00", sky: "#56B4E9", green: "#009E73", blue: "#0072B2", red: "#D55E00", purple: "#CC79A7", grey: "#999999", ink: "#222" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-wrap: nowrap; gap: .8em; align-items: flex-start; }
[data-demo="${NAME}"] .nq-ctrl { flex: 0 0 360px; display: flex; flex-direction: column; gap: .3em; }
[data-demo="${NAME}"] .nq-ctrl .qb-slider { grid-template-columns: 7.6em 1fr 4.3em; gap: .35em; }
[data-demo="${NAME}"] .nq-row { display: flex; flex-wrap: wrap; gap: .2em 1em; align-items: center; }
[data-demo="${NAME}"] .nq-plot { flex: 1 1 600px; min-width: 0; }
[data-demo="${NAME}"] .nq-plot svg { width: 100%; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 19px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .qb-readout { font-size: .85em; line-height: 1.3; }
[data-demo="${NAME}"] .nq-warn { color: #b00020; font-weight: 600; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 18px; }
.reveal [data-demo="${NAME}"] .nq-ctrl { flex-basis: 330px; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 16px; }
.reveal [data-demo="${NAME}"] .nq-plot svg { max-height: 530px; }
.reveal [data-demo="${NAME}"] .qb-btn { font-size: 16px; padding: .15em .5em; }
@media (max-width: 900px) {
  [data-demo="${NAME}"].qb-demo { flex-wrap: wrap; }
  [data-demo="${NAME}"] .nq-ctrl { flex-basis: 100%; }
}
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- the lecture's circadian dataset (t = 0, 2, ..., 46) ----
const DATA = [7.81, 8.81, 9.05, 9.98, 9.40, 8.50, 8.15, 7.47, 6.87, 6.41, 7.15, 7.37,
  7.81, 8.09, 9.64, 9.49, 9.29, 9.03, 8.25, 7.43, 6.98, 6.73, 6.72, 6.65];
const MEAN = 8, A24 = 1.5, PEAK = 6, SD0 = 0.3;
const TMAX = 120, GRIDH = 0.5;           // noise lives on a 0.5-h grid, 0..120 h
// standard-normal noise per grid time; at the 2-h times of the dataset it is the
// dataset's residual / 0.3, so sd = 0.3 reproduces the data exactly
const Z: number[] = (() => {
  const r = rng(17);
  const z = Array.from({ length: Math.round(TMAX / GRIDH) + 1 }, () => r.norm());
  DATA.forEach((x, i) => {
    const t = 2 * i;
    z[Math.round(t / GRIDH)] = (x - (MEAN + A24 * Math.cos(2 * Math.PI * (t - PEAK) / 24))) / SD0;
  });
  return z;
})();

interface State {
  dt: number; dur: number; P: number; Pu: number; Bu: number; sd: number; cut: number;
  win: boolean; filt: boolean; log: boolean;
}

const num = (s: string | undefined, d: number, lo: number, hi: number): number => {
  const v = Number(s);
  return s !== undefined && s !== "" && Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : d;
};

function dft(x: number[]): { re: number[]; im: number[] } {
  const n = x.length, re = new Array(n).fill(0), im = new Array(n).fill(0);
  for (let k = 0; k <= n / 2; k++) {
    let a = 0, b = 0;
    for (let j = 0; j < n; j++) {
      const th = 2 * Math.PI * ((j * k) % n) / n;
      a += x[j] * Math.cos(th); b -= x[j] * Math.sin(th);
    }
    re[k] = a; im[k] = b;
    if (k > 0 && k < n - k) { re[n - k] = a; im[n - k] = -b; }
  }
  return { re, im };
}

const fmt = (v: number, d = 2) => v.toFixed(d);
const fmtP = (p: number) => Math.abs(p - Math.round(p)) < 0.005 ? String(Math.round(p)) : p.toFixed(1);

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const init: State = {
    dt: num(ds.dt, 2, 0.5, 12), dur: num(ds.dur, 48, 24, TMAX),
    P: num(ds.period, 24, 20, 28), Pu: num(ds.fastPeriod, 4.8, 2, 12), Bu: num(ds.fastAmp, 0, 0, 1.5),
    sd: num(ds.noise, SD0, 0, 1.5), cut: num(ds.cut, 8, 2, 48),
    win: ds.window === "1", filt: ds.filter === "1", log: ds.log === "1",
  };
  const st: State = { ...init };

  // ---- controls ----
  const ctrl = htmlEl("div", "nq-ctrl", "", host);
  const sl: Record<string, Slider> = {};
  sl.dt = slider(ctrl, "sample every", 0.5, 12, 0.5, st.dt, v => { st.dt = v; draw(); }, v => v.toFixed(1) + " h");
  sl.dur = slider(ctrl, "duration", 24, TMAX, 1, st.dur, v => { st.dur = v; draw(); }, v => v + " h");
  sl.P = slider(ctrl, "clock period", 20, 28, 0.5, st.P, v => { st.P = v; draw(); }, v => v.toFixed(1) + " h");
  sl.Pu = slider(ctrl, "fast period", 2, 12, 0.1, st.Pu, v => { st.Pu = v; draw(); }, v => v.toFixed(1) + " h");
  sl.Bu = slider(ctrl, "fast amplitude", 0, 1.5, 0.1, st.Bu, v => { st.Bu = v; draw(); }, v => v.toFixed(1));
  sl.sd = slider(ctrl, "noise sd", 0, 1.5, 0.05, st.sd, v => { st.sd = v; draw(); }, v => v.toFixed(2));
  sl.cut = slider(ctrl, "keep periods ≥", 2, 48, 0.5, st.cut, v => { st.cut = v; draw(); }, v => v.toFixed(1) + " h");
  const cr = htmlEl("div", "nq-row", "", ctrl);
  const cbWin = checkbox(cr, "Hann window", st.win, v => { st.win = v; draw(); });
  const cbFilt = checkbox(cr, "low-pass", st.filt, v => { st.filt = v; draw(); });
  const cbLog = checkbox(cr, "log scale", st.log, v => { st.log = v; draw(); });

  const apply = (p: Partial<State>) => {
    Object.assign(st, init, { win: false, filt: false, log: false }, p);
    for (const k of Object.keys(sl)) sl[k].set((st as unknown as Record<string, number>)[k]);
    cbWin.checked = st.win; cbFilt.checked = st.filt; cbLog.checked = st.log;
    draw();
  };
  const base = { dt: 2, dur: 48, P: 24, Pu: 4.8, Bu: 0, sd: SD0, cut: 8 };
  const pb = presets(ctrl, [
    { label: "course data", apply: () => apply({ ...base }) },
    { label: "+ 4.8 h, every 2 h", apply: () => apply({ ...base, Bu: 1, sd: 0 }) },
    { label: "every 4 h", apply: () => apply({ ...base, Bu: 1, sd: 0, dt: 4 }) },
    { label: "leakage", apply: () => apply({ ...base, dt: 0.5, dur: 120, P: 25.5, sd: 0, log: true }) },
    { label: "low-pass", apply: () => apply({ ...base, dt: 0.5, dur: 96, sd: 1, Bu: 0, filt: true, cut: 8 }) },
    { label: "reset", apply: () => { apply({ ...init }); pb.forEach(b => b.classList.remove("active")); } },
  ]);
  const ro = htmlEl("div", "qb-readout", "", ctrl);

  // ---- plot ----
  const W = 800, H = 560;
  const plot = htmlEl("div", "nq-plot", "", host);
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Sampled signal and its amplitude spectrum" }, plot);
  const L = 66, R = W - 14;
  const tTop = 62, tBot = 250;            // time panel
  const sTop = 352, sBot = H - 48;        // spectrum panel
  const FMAX = 12;                        // cycles per day
  const fx = scale(0, FMAX, L, R);

  function draw(): void {
    const { dt, dur, P, Pu, Bu, sd } = st;
    const n = Math.max(2, Math.floor(dur / dt + 1e-9));
    const T = n * dt;
    const truth = (t: number) => MEAN + A24 * Math.cos(2 * Math.PI * (t - PEAK) / P) + Bu * Math.cos(2 * Math.PI * t / Pu);
    const ts = Array.from({ length: n }, (_, j) => j * dt);
    const xs = ts.map(t => truth(t) + sd * Z[Math.round(t / GRIDH)]);
    const w = ts.map((_, j) => st.win ? 0.5 - 0.5 * Math.cos(2 * Math.PI * j / n) : 1);
    const sw = w.reduce((a, b) => a + b, 0);
    const X = dft(xs);                                  // for the time-domain curve
    const xbar = X.re[0] / n;                         // remove the mean before windowing,
    const Xw = st.win ? dft(xs.map((x, j) => (x - xbar) * w[j])) : X;   // or it leaks into k = 1
    const half = Math.floor(n / 2);
    const amp: number[] = [];
    for (let k = 0; k <= half; k++) {
      const m = Math.hypot(Xw.re[k], Xw.im[k]) / sw;
      amp.push(k === 0 || (n % 2 === 0 && k === half) ? m : 2 * m);
    }
    const kmax = st.filt ? Math.floor(T / st.cut + 1e-9) : half;
    // trigonometric interpolant / low-pass curve: real form with k = 0..half
    const curve = (t: number): number => {
      let s = X.re[0] / n;
      for (let k = 1; k <= Math.min(half, kmax); k++) {
        const th = 2 * Math.PI * k * t / T;
        const c = (n % 2 === 0 && k === half) ? 1 : 2;
        s += c * (X.re[k] * Math.cos(th) - X.im[k] * Math.sin(th)) / n;
      }
      return s;
    };

    svg.replaceChildren();
    // ---------- time panel ----------
    const tx = scale(0, dur, L, R);
    const hspan = Math.max(2, Math.ceil(A24 + Bu + 2.5 * sd));
    const ty = scale(MEAN - hspan, MEAN + hspan, tBot, tTop);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svg).textContent =
      `Time: n = ${n} samples every ${fmtP(dt)} h, record T = n·Δt = ${fmtP(T)} h`;
    const ystep = hspan <= 2 ? 1 : 2;
    for (let v = MEAN - hspan; v <= MEAN + hspan + 1e-9; v += ystep) {
      svgEl("line", { x1: L, x2: R, y1: ty(v), y2: ty(v), class: "grid" }, svg);
      svgEl("text", { x: L - 6, y: ty(v) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = String(v);
    }
    const tstep = dur > 72 ? 24 : 12;
    for (let t = 0; t <= dur + 1e-9; t += tstep) {
      svgEl("line", { x1: tx(t), x2: tx(t), y1: tBot, y2: tBot + 5, class: "axis" }, svg);
      svgEl("text", { x: tx(t), y: tBot + 22, "text-anchor": "middle", class: "sm" }, svg).textContent = String(t);
    }
    svgEl("line", { x1: L, x2: R, y1: tBot, y2: tBot, class: "axis" }, svg);
    svgEl("text", { x: R, y: tBot + 44, "text-anchor": "end" }, svg).textContent = "time (h)";
    const yl = svgEl("text", { x: 18, y: (tTop + tBot) / 2, "text-anchor": "middle", transform: `rotate(-90 18 ${(tTop + tBot) / 2})` }, svg);
    yl.textContent = "log2 expr.";
    const NP = 700;
    const line = (f: (t: number) => number, t1: number, attrs: Record<string, string | number>) => {
      const pts: string[] = [];
      for (let i = 0; i <= NP; i++) {
        const t = t1 * i / NP;
        const y = Math.max(MEAN - hspan - 0.3, Math.min(MEAN + hspan + 0.3, f(t)));
        pts.push(`${tx(t).toFixed(1)},${ty(y).toFixed(1)}`);
      }
      svgEl("polyline", { points: pts.join(" "), fill: "none", ...attrs }, svg);
    };
    line(truth, dur, { stroke: OI.grey, "stroke-width": 2 });
    line(curve, T, st.filt ? { stroke: OI.green, "stroke-width": 3.5 } : { stroke: OI.orange, "stroke-width": 3, "stroke-dasharray": "9 5" });
    const rmax = n > 120 ? 3 : n > 60 ? 4 : 6;
    ts.forEach((t, j) => svgEl("circle", { cx: tx(t), cy: ty(Math.max(MEAN - hspan - 0.3, Math.min(MEAN + hspan + 0.3, xs[j]))), r: rmax, fill: OI.blue }, svg));
    // legend
    const lg = svgEl("g", { transform: `translate(${L + 8}, 44)` }, svg);
        const items: [string, Record<string, string | number>][] = [
      ["true signal", { stroke: OI.grey, "stroke-width": 2 }],
      ["samples", {}],
      [st.filt ? "low-pass curve" : "Fourier curve through samples", st.filt ? { stroke: OI.green, "stroke-width": 3.5 } : { stroke: OI.orange, "stroke-width": 3, "stroke-dasharray": "9 5" }],
    ];
    let lx = 0;
    items.forEach(([lab, a], i) => {
      if (i === 1) svgEl("circle", { cx: lx + 12, cy: 0, r: 6, fill: OI.blue }, lg);
      else svgEl("line", { x1: lx, x2: lx + 26, y1: 0, y2: 0, ...a }, lg);
      svgEl("text", { x: lx + 32, y: 6 }, lg).textContent = lab;
      lx += 52 + lab.length * 9.5;
    });

    // ---------- spectrum panel ----------
    const fN = 24 / (2 * dt), df = 24 / T;           // cycles per day
    let kpk = 1;
    for (let k = 1; k <= half; k++) if (amp[k] > amp[kpk]) kpk = k;
    svgEl("text", { x: 8, y: sTop - 34, class: "ttl" }, svg).textContent =
      `Spectrum${st.win ? " (Hann)" : ""}: Nyquist period 2Δt = ${fmtP(2 * dt)} h; bins at ${fmtP(T)}/k h`;
    const amax = Math.max(...amp.slice(1), 0.1);
    const yTopLin = amax > 2 ? Math.ceil(amax * 1.1 / 0.5) * 0.5 : amax > 1 ? 2 : 1;
    const LY0 = -4, LY1 = 1;
    const sy = st.log
      ? (a: number) => scale(LY0, LY1, sBot, sTop)(Math.max(LY0, Math.log10(Math.max(a, 1e-12))))
      : (a: number) => scale(0, yTopLin, sBot, sTop)(Math.min(a, yTopLin));
    // beyond Nyquist
    if (fN < FMAX) {
      svgEl("rect", { x: fx(fN), y: sTop, width: R - fx(fN), height: sBot - sTop, fill: "#000", opacity: 0.06 }, svg);
      svgEl("text", { x: (fx(fN) + R) / 2, y: sBot - 10, "text-anchor": "middle", class: "sm" }, svg).textContent =
        R - fx(fN) > 200 ? "beyond Nyquist: folds back" : R - fx(fN) > 90 ? "folds back" : "";
    }
    // y grid
    if (st.log) {
      for (let e = LY0; e <= LY1; e++) {
        svgEl("line", { x1: L, x2: R, y1: sy(10 ** e), y2: sy(10 ** e), class: "grid" }, svg);
        svgEl("text", { x: L - 6, y: sy(10 ** e) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = e >= -2 ? String(10 ** e) : `1e${e}`;
      }
    } else {
      const ys = yTopLin <= 1 ? 0.25 : 0.5;
      for (let v = 0; v <= yTopLin + 1e-9; v += ys) {
        svgEl("line", { x1: L, x2: R, y1: sy(v), y2: sy(v), class: v === 0 ? "axis" : "grid" }, svg);
        svgEl("text", { x: L - 6, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svg).textContent = String(v);
      }
    }
    for (let f = 0; f <= FMAX; f++) {
      svgEl("line", { x1: fx(f), x2: fx(f), y1: sBot, y2: sBot + 5, class: "axis" }, svg);
      svgEl("text", { x: fx(f), y: sBot + 22, "text-anchor": "middle", class: "sm" }, svg).textContent = String(f);
    }
    svgEl("line", { x1: L, x2: R, y1: sBot, y2: sBot, class: "axis" }, svg);
    svgEl("text", { x: R, y: sBot + 44, "text-anchor": "end" }, svg).textContent = "frequency (cycles per day)";
    // period labels on top
    for (const p of [24, 12, 8, 6, 4.8, 4, 3, 2]) {
      svgEl("line", { x1: fx(24 / p), x2: fx(24 / p), y1: sTop, y2: sTop - 5, class: "axis" }, svg);
      svgEl("text", { x: fx(24 / p), y: sTop - 9, "text-anchor": "middle", class: "sm" }, svg).textContent = p + " h";
    }
    svgEl("line", { x1: L, x2: R, y1: sTop, y2: sTop, class: "grid" }, svg);
    const ylab = svgEl("text", { x: 18, y: (sTop + sBot) / 2, "text-anchor": "middle", transform: `rotate(-90 18 ${(sTop + sBot) / 2})` }, svg);
    ylab.textContent = "amplitude";
    // true frequencies
    const fc = 24 / P, fu = 24 / Pu;
    const trueMark = (f: number, col: string, lab: string) => {
      svgEl("line", { x1: fx(f), x2: fx(f), y1: sTop, y2: sBot, stroke: col, "stroke-width": 2.5, "stroke-dasharray": "3 4" }, svg);
      svgEl("text", { x: fx(f) + 5, y: sTop + 20, fill: col, style: `fill:${col}` }, svg).textContent = lab;
    };
    trueMark(fc, OI.orange, "clock");
    if (Bu > 0) {
      trueMark(fu, OI.purple, "fast");
      if (fu > fN + 1e-9) {
        const fs = 24 / dt;
        const fa = Math.abs(fu - Math.round(fu / fs) * fs);
        const x1 = fx(fu), x2 = fx(fa), yA = sBot - 50;
        const mid = (x1 + x2) / 2;
        svgEl("path", { d: `M ${x1} ${yA} Q ${mid} ${yA - 30} ${x2} ${yA}`, fill: "none", stroke: OI.purple, "stroke-width": 3 }, svg);
        svgEl("path", { d: `M ${x2} ${yA} l ${x1 > x2 ? 10 : -10} -9 m 0 9 l ${x1 > x2 ? 2 : -2} -13`, fill: "none", stroke: OI.purple, "stroke-width": 2.5 }, svg);
        svgEl("text", { x: mid, y: yA - 24, "text-anchor": "middle", style: `fill:${OI.purple}`, "font-weight": 600 }, svg)
          .textContent = fa < 1e-9 ? "alias: looks constant" : `alias: ${fmtP(24 / fa)} h`;
      }
    }
    // Nyquist line
    if (fN <= FMAX) {
      svgEl("line", { x1: fx(fN), x2: fx(fN), y1: sTop, y2: sBot, stroke: OI.ink, "stroke-width": 2.5 }, svg);
      svgEl("text", { x: fx(fN) - 5, y: sTop + 44, "text-anchor": "end", "font-weight": 600 }, svg).textContent = "Nyquist";
    } else {
      svgEl("text", { x: R - 4, y: sTop + 44, "text-anchor": "end", "font-weight": 600 }, svg).textContent = `Nyquist at ${fmtP(fN)}/day →`;
    }
    // stems (k >= 1)
    const col = st.win ? OI.blue : OI.red;
    const y0 = st.log ? sBot : sy(0);
    for (let k = 1; k <= half; k++) {
      const f = k * df;
      if (f > FMAX + 1e-9) break;
      const kept = !st.filt || k <= kmax;
      const a = amp[k];
      const c = kept ? col : "#bbb";
      svgEl("line", { x1: fx(f), x2: fx(f), y1: y0, y2: sy(a), stroke: c, "stroke-width": n > 120 ? 2 : 3 }, svg);
      svgEl("circle", { cx: fx(f), cy: sy(a), r: n > 120 ? 3 : 5, fill: c }, svg);
    }
    if (st.filt) {
      const fcut = 24 / st.cut;
      if (fcut <= FMAX) {
        svgEl("line", { x1: fx(fcut), x2: fx(fcut), y1: sTop + 50, y2: sBot, stroke: OI.green, "stroke-width": 3 }, svg);
        svgEl("text", { x: fx(fcut) + 5, y: sTop + 66, style: `fill:${OI.green}`, "font-weight": 600 }, svg).textContent = "cut-off";
      }
    }

    // ---------- readout ----------
    const amp1 = (k: number) => k < 1 || k > half ? NaN : amp[k];
    const peakTime = (k: number) => {
      const per = T / k;
      const ph = -Math.atan2(Xw.im[k], Xw.re[k]);
      return ((ph / (2 * Math.PI) * per) % per + per) % per;
    };
    const lines: string[] = [];
    lines.push(`mean ${fmt(X.re[0] / n, 3)}; strongest k = ${kpk}: period ${fmtP(T / kpk)} h, amplitude ${fmt(amp[kpk], 3)}`);
    const k24 = Math.round(T / 24);
    if (k24 >= 1 && k24 <= half) {
      const onGrid = Math.abs(T / k24 - 24) < 1e-9;
      lines.push(`${onGrid ? "24-h bin" : `nearest bin to 24 h (${fmtP(T / k24)} h)`}, k = ${k24}: amplitude ${fmt(amp1(k24), 3)}, peak at ${fmt(peakTime(k24), 2)} h`);
    }
    let warn = "";
    if (Bu > 0 && fu > fN + 1e-9) {
      const fs = 24 / dt, fa = Math.abs(fu - Math.round(fu / fs) * fs);
      warn = `${fmtP(Pu)}-h rhythm is faster than Nyquist (${fmtP(2 * dt)} h): it shows up at ${fa < 1e-9 ? "period ∞ (the mean)" : fmtP(24 / fa) + " h"}`;
    } else if (Bu > 0 && Math.abs(fu - fN) < 1e-9) warn = `${fmtP(Pu)}-h rhythm sits exactly at Nyquist`;
    if (Math.abs(T / P - Math.round(T / P)) > 1e-6) lines.push(`T/P = ${(T / P).toFixed(2)} cycles of the clock: not whole → leakage`);
    else lines.push(`T/P = ${Math.round(T / P)} whole cycles of the clock`);
    if (st.filt) {
      const rms = (f: (j: number) => number) => Math.sqrt(ts.reduce((s, t, j) => s + f(j) ** 2, 0) / n);
      lines.push(`low-pass keeps |k| ≤ ${Math.min(kmax, half)} (${Math.min(n, 2 * Math.min(kmax, half) + 1)} of ${n} coefficients); rms error vs truth: samples ${fmt(rms(j => xs[j] - truth(ts[j])))}, smoothed ${fmt(rms(j => curve(ts[j]) - truth(ts[j])))}`);
    }
    ro.replaceChildren();
    lines.forEach(s => htmlEl("div", "", s, ro));
    if (warn) htmlEl("div", "nq-warn", warn, ro);
  }
  draw();
}

autoMount(NAME, mount);

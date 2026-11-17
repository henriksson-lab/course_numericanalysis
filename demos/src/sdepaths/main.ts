// Sample paths of an SDE for the sde unit ("Stochastic differential equations").
// Ornstein-Uhlenbeck  dX = θ(μ − X) dt + σ dW   (additive noise), or
// geometric Brownian motion  dX = r X dt + σ X dW   (multiplicative noise).
// (a) ~20 paths from x0 by Euler-Maruyama or Milstein with step h. Every path's
//     Brownian motion lives on a fine grid of 2^16 steps; a coarse increment ΔW is the
//     sum of the fine ones inside the step, so "same path, finer step" keeps the path.
//     Path 1 is highlighted together with the exact solution for the same W (GBM: the
//     closed form; OU: the exact OU recursion on the fine grid), and the largest gap.
// (b) a histogram of N paths at time t (sideways, sharing the y axis of (a)) against
//     the exact density (the Fokker-Planck solution): normal for OU, lognormal for GBM.
// For OU with θh ≥ 2 Euler-Maruyama is not mean-square stable and the paths blow up.
// Auto-mounts on <div data-demo="sdepaths">. Optional attributes:
//   data-model="ou|gbm"  data-scheme="em|milstein"  data-theta="1"  data-r="1"
//   data-sigma="1"  data-mu="0"  data-x0="3"  data-steps="32" (a power of 2, 4..1024)
//   data-n="2000" (histogram paths)  data-paths="20"  data-t="5" (histogram time)
//   data-seed="1"
import { svgEl, htmlEl, slider, button, presets, scale, rng, autoMount } from "../lib/ui";

const NAME = "sdepaths";
// Okabe-Ito: numerical paths sky blue / blue, exact vermillion, time cursor green
const COL = {
  path: "#56B4E9", hi: "#0072B2", exact: "#D55E00", bar: "#56B4E9", cursor: "#009E73",
  ink: "#222", grey: "#737373", warn: "#b00020",
};
const KF = 16, NF = 1 << KF;           // fine Brownian grid
const KMIN = 2, KMAX = 10, NMAX = 1 << KMAX;
const FPC = NF / NMAX;                 // fine steps per finest coarse step
const NVALS = [200, 500, 1000, 2000, 5000, 10000];
const NSNAP = 100;                     // histogram snapshots in [0, T]
const TEND = { ou: 5, gbm: 1 };
type Model = "ou" | "gbm";
type Scheme = "em" | "mil";

const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .35em; }
[data-demo="${NAME}"] .sp-row { display: flex; flex-wrap: wrap; gap: .3em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .sp-group { display: flex; gap: .35em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .sp-row .qb-slider { flex: 1 1 200px; min-width: 0; grid-template-columns: auto minmax(3em, 1fr) 4.2em; gap: .35em; }
[data-demo="${NAME}"] .sp-row .qb-slider input[type=range] { min-width: 3em; }
[data-demo="${NAME}"] .sp-plots { display: flex; gap: 6px; width: 100%; }
[data-demo="${NAME}"] .sp-plots svg { min-width: 0; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 20px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
[data-demo="${NAME}"] .sp-read { display: flex; flex-direction: column; gap: .1em; width: 100%; box-sizing: border-box; }
[data-demo="${NAME}"] .sp-read { font-family: system-ui, -apple-system, "Segoe UI", sans-serif; white-space: normal; }
[data-demo="${NAME}"] .sp-read b { font-weight: 600; }
[data-demo="${NAME}"] .sp-warn { color: ${COL.warn}; font-weight: 600; font-family: system-ui, sans-serif; }
[data-demo="${NAME}"] .sp-note { font-size: .85em; color: #555; }
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

const minus = (s: string): string => s.replace(/^-/, "−");
const fx = (v: number, d: number): string => Number.isFinite(v) ? minus(v.toFixed(d)) : "∞";
const fe = (v: number): string => {
  if (!Number.isFinite(v)) return "∞";
  const a = Math.abs(v);
  if (a < 1e-12) return "0";
  if (a !== 0 && (a >= 1e4 || a < 1e-3)) return minus(v.toExponential(2)).replace("e+", "e");
  return minus(v.toPrecision(3));
};
const fh = (v: number): string => v >= 0.1 ? v.toFixed(3).replace(/0+$/, "").replace(/\.$/, "") : v.toPrecision(3);

function niceStep(range: number, target: number): number {
  const raw = range / target, p = 10 ** Math.floor(Math.log10(raw)), m = raw / p;
  return (m < 1.5 ? 1 : m < 3 ? 2 : m < 7 ? 5 : 10) * p;
}

/** Normal density with mean m and variance v. */
const npdf = (x: number, m: number, v: number): number => Math.exp(-((x - m) ** 2) / (2 * v)) / Math.sqrt(2 * Math.PI * v);

let clipCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const num = (s: string | undefined, d: number, lo: number, hi: number): number => {
    const v = Number(s);
    return s !== undefined && s !== "" && Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : d;
  };
  const model0: Model = ds.model === "gbm" ? "gbm" : "ou";
  const nearestN = (n: number): number => NVALS.reduce((b, v) => Math.abs(Math.log(v / n)) < Math.abs(Math.log(b / n)) ? v : b, NVALS[0]);
  const start = {
    model: model0,
    scheme: (ds.scheme === "milstein" || ds.scheme === "mil" ? "mil" : "em") as Scheme,
    theta: num(ds.theta, 1, 0.1, 10),
    r: num(ds.r, 1, -4, 3),
    sigma: num(ds.sigma, 1, 0.1, 2),
    mu: num(ds.mu, 0, -5, 5),
    x0: { ou: model0 === "ou" ? num(ds.x0, 3, -5, 5) : 3, gbm: model0 === "gbm" ? num(ds.x0, 1, 0.1, 5) : 1 },
    k: Math.round(Math.log2(num(ds.steps, 32, 1 << KMIN, NMAX))),
    N: nearestN(num(ds.n, 2000, 200, 10000)),
    j: NSNAP,
    noise: 0,
  };
  if (ds.t !== undefined && ds.t !== "") start.j = Math.round(NSNAP * num(ds.t, TEND[model0], 0, TEND[model0]) / TEND[model0]);
  const P = Math.round(num(ds.paths, 20, 2, 50));
  const seed = Math.round(num(ds.seed, 1, 0, 1e9));
  const st = { ...start, x0: { ...start.x0 } };
  const T = (): number => TEND[st.model];
  const nSteps = (): number => 1 << st.k;
  const hh = (): number => T() / nSteps();
  const x0 = (): number => st.x0[st.model];

  // ---------------------------------------------------------------- controls
  const row1 = htmlEl("div", "sp-row", "", host);
  const modG = htmlEl("div", "sp-group", "", row1);
  htmlEl("span", "", "SDE:", modG);
  const modBtns = presets(modG, [
    { label: "Ornstein–Uhlenbeck", apply: () => { setModel("ou"); } },
    { label: "geometric BM", apply: () => { setModel("gbm"); } },
  ]);
  const schG = htmlEl("div", "sp-group", "", row1);
  htmlEl("span", "", "scheme:", schG);
  const schBtns = presets(schG, [
    { label: "Euler–Maruyama", apply: () => { st.scheme = "em"; changed(); } },
    { label: "Milstein", apply: () => { st.scheme = "mil"; changed(); } },
  ]);
  const schNote = htmlEl("span", "sp-note", "", row1);

  const row2 = htmlEl("div", "sp-row", "", host);
  const hS = slider(row2, "step h", KMIN, KMAX, 1, st.k, v => { st.k = v; changed(); }, v => fh(T() / (1 << v)));
  const sS = slider(row2, "noise σ", 0.1, 2, 0.05, st.sigma, v => { st.sigma = v; changed(); });
  const dS = slider(row2, "θ", 0.1, 10, 0.1, st.theta, v => { if (st.model === "ou") st.theta = v; else st.r = v; changed(); });
  const dName = dS.input.parentElement!.querySelector(".qb-slider-name") as HTMLElement;
  const nS = slider(row2, "paths in (b)", 0, NVALS.length - 1, 1, NVALS.indexOf(st.N), v => { st.N = NVALS[v]; changed(); }, v => String(NVALS[v]));

  const row3 = htmlEl("div", "sp-row", "", host);
  const bG = htmlEl("div", "sp-group", "", row3);
  const finerB = button(bG, "same path, finer step h/2", () => { if (st.k < KMAX) { st.k++; hS.set(st.k); changed(); } });
  const coarserB = button(bG, "coarser 2h", () => { if (st.k > KMIN) { st.k--; hS.set(st.k); changed(); } });
  button(bG, "new noise", () => { st.noise++; changed(); });
  const tS = slider(row3, "time t", 0, NSNAP, 1, st.j, v => { stopPlay(); st.j = v; draw(); }, v => fx(tAt(v), 2));
  const playB = button(row3, "▶ play", () => { if (raf) stopPlay(); else startPlay(); });

  const row4 = htmlEl("div", "sp-row", "", host);
  const pG = htmlEl("div", "sp-group", "", row4);
  htmlEl("span", "", "presets:", pG);
  presets(pG, [
    { label: "OU relaxes", apply: () => applyPreset({ model: "ou", scheme: "em", theta: 1, sigma: 1, k: 5, j: NSNAP }) },
    { label: "OU, θh = 2.5 (unstable)", apply: () => applyPreset({ model: "ou", scheme: "em", theta: 8, sigma: 1, k: 4, j: NSNAP }) },
    { label: "GBM: EM vs Milstein", apply: () => applyPreset({ model: "gbm", scheme: "em", r: 1, sigma: 1, k: 4, j: NSNAP }) },
  ]);
  button(row4, "reset", () => {
    stopPlay();
    Object.assign(st, start, { x0: { ...start.x0 } });
    syncSliders(); changed();
  });

  const plots = htmlEl("div", "sp-plots", "", host);
  const H = 288, WA = 720, WB = 380, WT = WA + WB;
  const svgA = svgEl("svg", { viewBox: `0 0 ${WA} ${H}`, role: "img", "aria-label": "Sample paths against time" }, plots);
  svgA.style.width = `calc(${(100 * WA / WT).toFixed(2)}% - 3px)`;
  const svgB = svgEl("svg", { viewBox: `0 0 ${WB} ${H}`, role: "img", "aria-label": "Histogram of X at time t with the exact density" }, plots);
  svgB.style.width = `calc(${(100 * WB / WT).toFixed(2)}% - 3px)`;
  const ro = htmlEl("div", "qb-readout sp-read", "", host);
  const clipId = `sdp-clip-${++clipCount}`;

  function applyPreset(p: Partial<typeof st>): void {
    stopPlay();
    Object.assign(st, p);
    if (p.model) st.x0[p.model] = start.x0[p.model] ?? st.x0[p.model];
    st.mu = start.mu;
    syncSliders(); changed();
  }
  function setModel(m: Model): void {
    stopPlay();
    st.model = m;
    syncSliders(); changed();
  }
  function syncSliders(): void {
    hS.set(st.k); sS.set(st.sigma); nS.set(NVALS.indexOf(st.N)); tS.set(st.j);
    if (st.model === "ou") {
      dName.textContent = "drift θ";
      Object.assign(dS.input, { min: "0.1", max: "10", step: "0.1" });
      dS.set(st.theta);
    } else {
      dName.textContent = "growth r";
      Object.assign(dS.input, { min: "-4", max: "3", step: "0.1" });
      dS.set(st.r);
    }
  }
  function syncControls(): void {
    modBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.model === "ou")));
    schBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (st.scheme === "em")));
    finerB.disabled = st.k >= KMAX;
    coarserB.disabled = st.k <= KMIN;
    schNote.textContent = st.model === "ou"
      ? "additive noise: Milstein = EM here (order 1)"
      : "Milstein adds ½σ²X(ΔW² − h): order 1 vs ½";
  }

  // ---------------------------------------------------------------- numerics
  /** One step of the scheme. */
  function step(x: number, h: number, dW: number): number {
    if (st.model === "ou") return x + st.theta * (st.mu - x) * h + st.sigma * dW;
    let y = x + st.r * x * h + st.sigma * x * dW;
    if (st.scheme === "mil") y += 0.5 * st.sigma * st.sigma * x * (dW * dW - h);
    return y;
  }

  // fine noise of the shown paths (standard normals), and its sums on the finest coarse grid
  let Zf: Float32Array[] = [], S: Float64Array[] = [], zKey = "";
  function ensureNoise(): void {
    const key = String(st.noise);
    if (key === zKey) return;
    zKey = key;
    const R = rng(seed * 7919 + st.noise * 101 + 3);
    Zf = []; S = [];
    for (let p = 0; p < P; p++) {
      const z = new Float32Array(NF), s = new Float64Array(NMAX + 1);
      let acc = 0;
      for (let i = 0; i < NF; i++) {
        z[i] = R.norm(); acc += z[i];
        if ((i + 1) % FPC === 0) s[(i + 1) / FPC] = acc;
      }
      Zf.push(z); S.push(s);
    }
  }

  // exact solution on the finest coarse grid (NMAX + 1 points) for every shown path
  let EX: Float64Array[] = [], exKey = "", yLo = 0, yHi = 1;
  function ensureExact(): void {
    const key = [zKey, st.model, st.theta, st.r, st.sigma, st.mu, x0()].join("|");
    if (key === exKey) return;
    exKey = key;
    const dt = T() / NF, sq = Math.sqrt(dt);
    EX = [];
    let lo = Infinity, hi = -Infinity;
    for (let p = 0; p < P; p++) {
      const z = Zf[p], ex = new Float64Array(NMAX + 1);
      ex[0] = x0();
      if (st.model === "ou") {
        // exact OU transition over dt; given ΔW its best estimate of ∫e^{−θ(t−s)}dW is ΔW(1−e)/(θdt)
        const e = Math.exp(-st.theta * dt), c = st.sigma * sq * (1 - e) / (st.theta * dt);
        let x = x0();
        for (let i = 0; i < NF; i++) {
          x = st.mu + (x - st.mu) * e + c * z[i];
          if ((i + 1) % FPC === 0) ex[(i + 1) / FPC] = x;
        }
      } else {
        const s = S[p], a = st.r - st.sigma * st.sigma / 2, H0 = T() / NMAX;
        for (let i = 1; i <= NMAX; i++) ex[i] = x0() * Math.exp(a * i * H0 + st.sigma * sq * s[i]);
      }
      for (const v of ex) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
      EX.push(ex);
    }
    if (st.model === "ou") {
      const sInf = st.sigma / Math.sqrt(2 * st.theta);
      lo = Math.min(lo, st.mu - 3 * sInf, x0()); hi = Math.max(hi, st.mu + 3 * sInf, x0());
      const pad = 0.06 * (hi - lo);
      yLo = lo - pad; yHi = hi + pad;
    } else {
      const a = Math.log(x0()) + (st.r - st.sigma ** 2 / 2) * T(), s = st.sigma * Math.sqrt(T());
      const q995 = Math.exp(a + 2.576 * s);
      yLo = 0; yHi = Math.max(1.15 * x0(), 1.05 * Math.min(hi, Math.max(q995, x0())));
    }
  }

  // numerical shown paths at the current step
  let NUM: Float64Array[] = [];
  function computeNum(): void {
    const n = nSteps(), h = hh(), m = NMAX / n, sq = Math.sqrt(T() / NF);
    NUM = [];
    for (let p = 0; p < P; p++) {
      const s = S[p], x = new Float64Array(n + 1);
      x[0] = x0();
      for (let i = 0; i < n; i++) x[i + 1] = step(x[i], h, sq * (s[(i + 1) * m] - s[i * m]));
      NUM.push(x);
    }
  }

  // histogram paths: their own noise, drawn directly on the coarse grid; NSNAP+1 snapshots
  let HS = new Float64Array(0), hsKey = "";
  function snapIdx(j: number): number { return Math.floor(j * nSteps() / NSNAP); }
  function tAt(j: number): number { return snapIdx(j) * hh(); }
  function computeHist(): void {
    const key = [st.noise, st.model, st.scheme, st.theta, st.r, st.sigma, st.mu, x0(), st.k, st.N].join("|");
    if (key === hsKey) return;
    hsKey = key;
    const R = rng(seed * 104729 + st.noise * 13 + 5);
    const N = st.N, h = hh(), sq = Math.sqrt(h);
    HS = new Float64Array(N * (NSNAP + 1));
    for (let p = 0; p < N; p++) {
      let x = x0(), i = 0;
      for (let j = 0; j <= NSNAP; j++) {
        const target = snapIdx(j);
        while (i < target) { x = step(x, h, sq * R.norm()); i++; }
        HS[j * N + p] = x;
      }
    }
  }

  // strong errors for the same Brownian paths at different h
  let hist = new Map<number, number>(), histKey = "";
  let maxErr = 0, maxAt = 0, meanMax = 0;
  function errors(): void {
    const n = nSteps(), m = NMAX / n;
    maxErr = 0; maxAt = 0;
    const x = NUM[0], ex = EX[0];
    for (let i = 0; i <= n; i++) {
      const e = Math.abs(x[i] - ex[i * m]);
      if (!(e <= maxErr)) { maxErr = e; maxAt = i; if (!Number.isFinite(e)) break; }
    }
    // strong error estimate: mean over the shown paths of each path's largest gap
    let s1 = 0;
    for (let p = 0; p < P; p++) {
      let mp = 0;
      for (let i = 0; i <= n; i++) { const e = Math.abs(NUM[p][i] - EX[p][i * m]); if (!(e <= mp)) mp = e; }
      s1 += mp;
    }
    meanMax = s1 / P;
    const key = [exKey, st.scheme].join("|");
    if (key !== histKey) { histKey = key; hist = new Map(); }
    hist.set(st.k, meanMax);
  }

  function changed(): void {
    syncControls();
    ensureNoise(); ensureExact(); computeNum(); errors(); computeHist();
    draw();
  }

  // ---------------------------------------------------------------- play
  let raf = 0, lastTs = 0;
  function startPlay(): void {
    if (st.j >= NSNAP) st.j = 0;
    playB.textContent = "❚❚ pause";
    lastTs = 0;
    raf = window.requestAnimationFrame(tick);
  }
  function stopPlay(): void {
    if (raf) window.cancelAnimationFrame(raf);
    raf = 0; playB.textContent = "▶ play";
  }
  let jf = 0;
  function tick(ts: number): void {
    if (!lastTs) jf = st.j;
    const dt = lastTs ? Math.min(0.1, (ts - lastTs) / 1000) : 0;
    lastTs = ts;
    jf = Math.min(NSNAP, jf + dt * 20);
    st.j = Math.floor(jf); tS.set(st.j);
    draw();
    if (st.j >= NSNAP || !host.isConnected) { stopPlay(); return; }
    raf = window.requestAnimationFrame(tick);
  }

  // ---------------------------------------------------------------- drawing
  const L = 66, R0 = 22, top = 34, bot = H - 46;
  const clampPx = (v: number): number => Number.isFinite(v) ? Math.max(-3000, Math.min(H + 3000, v)) : (v > 0 ? -3000 : H + 3000);

  function draw(): void {
    const ay = scale(yLo, yHi, bot, top);
    const py = (v: number): number => Number.isNaN(v) ? -3000 : clampPx(ay(v));
    const n = nSteps(), h = hh(), Tn = T();
    const tj = tAt(st.j);

    // ---------- (a) paths ----------
    svgA.replaceChildren();
    const ax = scale(0, Tn, L, WA - R0);
    const defs = svgEl("defs", {}, svgA);
    const cp = svgEl("clipPath", { id: clipId }, defs);
    svgEl("rect", { x: L, y: top, width: WA - R0 - L, height: bot - top }, cp);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgA).textContent =
      `(a) ${P} paths, ${st.scheme === "em" ? "Euler–Maruyama" : "Milstein"}, h = ${fh(h)}`;
    const ys = niceStep(yHi - yLo, 5);
    for (let v = Math.ceil(yLo / ys) * ys; v <= yHi + 1e-9; v += ys) {
      svgEl("line", { x1: L, x2: WA - R0, y1: ay(v), y2: ay(v), class: Math.abs(v) < 1e-9 ? "axis" : "grid" }, svgA);
      svgEl("text", { x: L - 6, y: ay(v) + 6, "text-anchor": "end", class: "sm" }, svgA).textContent = fx(v, ys < 1 ? 1 : 0);
    }
    const ts = niceStep(Tn, 5);
    for (let t = 0; t <= Tn + 1e-9; t += ts) {
      svgEl("line", { x1: ax(t), x2: ax(t), y1: bot, y2: bot + 5, class: "axis" }, svgA);
      svgEl("text", { x: ax(t), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgA).textContent = fx(t, ts < 1 ? 1 : 0);
    }
    svgEl("line", { x1: L, x2: L, y1: top, y2: bot, class: "axis" }, svgA);
    svgEl("line", { x1: L, x2: WA - R0, y1: bot, y2: bot, class: "axis" }, svgA);
    svgEl("text", { x: (L + WA) / 2, y: H - 4, "text-anchor": "middle" }, svgA).textContent = "time t";
    const cy = (top + bot) / 2;
    svgEl("text", { x: 18, y: cy, "text-anchor": "middle", transform: `rotate(-90 18 ${cy})` }, svgA).textContent = "X(t)";
    if (st.model === "ou") {
      svgEl("line", { x1: L, x2: WA - R0, y1: ay(st.mu), y2: ay(st.mu), stroke: COL.grey, "stroke-width": 1.5, "stroke-dasharray": "6 5" }, svgA);
    }
    const g = svgEl("g", { "clip-path": `url(#${clipId})` }, svgA);
    const line = (xs: Float64Array, nn: number, attrs: Record<string, string | number>) => {
      const pts: string[] = [];
      for (let i = 0; i <= nn; i++) pts.push(`${ax(i * Tn / nn).toFixed(1)},${py(xs[i]).toFixed(1)}`);
      svgEl("polyline", { points: pts.join(" "), fill: "none", "stroke-linejoin": "round", ...attrs }, g);
    };
    for (let p = P - 1; p >= 1; p--) line(NUM[p], n, { stroke: COL.path, "stroke-width": 1.6, "stroke-opacity": 0.7 });
    line(EX[0], NMAX, { stroke: COL.exact, "stroke-width": 2.6 });
    line(NUM[0], n, { stroke: COL.hi, "stroke-width": 3 });
    if (n <= 64) for (let i = 0; i <= n; i++) svgEl("circle", { cx: ax(i * Tn / n), cy: py(NUM[0][i]), r: 3.5, fill: COL.hi }, g);
    // largest gap on the highlighted path
    if (maxErr > 0 && Number.isFinite(maxErr)) {
      const xm = ax(maxAt * Tn / n), y1 = py(NUM[0][maxAt]), y2 = py(EX[0][maxAt * NMAX / n]);
      svgEl("line", { x1: xm, x2: xm, y1, y2, stroke: COL.ink, "stroke-width": 2.5 }, g);
      for (const yy of [y1, y2]) svgEl("line", { x1: xm - 6, x2: xm + 6, y1: yy, y2: yy, stroke: COL.ink, "stroke-width": 2.5 }, g);
    }
    // time cursor
    svgEl("line", { x1: ax(tj), x2: ax(tj), y1: top, y2: bot, stroke: COL.cursor, "stroke-width": 2.5, "stroke-dasharray": "7 4" }, svgA);
    // legend
    const lgx = WA - 268, lgy = top + 4;
    svgEl("rect", { x: lgx - 8, y: lgy, width: 268, height: 76, fill: "#fff", "fill-opacity": 0.85, rx: 4 }, svgA);
    const leg = (i: number, col: string, w: number, txt: string) => {
      const y = lgy + 18 + 23 * i;
      svgEl("line", { x1: lgx, x2: lgx + 26, y1: y - 6, y2: y - 6, stroke: col, "stroke-width": w }, svgA);
      svgEl("text", { x: lgx + 34, y }, svgA).textContent = txt;
    };
    leg(0, COL.path, 2.5, "other paths");
    leg(1, COL.hi, 3, "path 1, numerical");
    leg(2, COL.exact, 3, "path 1, exact (same W)");

    // ---------- (b) histogram ----------
    svgB.replaceChildren();
    const N = st.N, base = st.j * N;
    const nb = 34, bw = (yHi - yLo) / nb, counts = new Array<number>(nb).fill(0);
    let off = 0, s1 = 0, s2 = 0, nf = 0;
    for (let p = 0; p < N; p++) {
      const v = HS[base + p];
      if (Number.isFinite(v)) { s1 += v; s2 += v * v; nf++; }
      const b = Math.floor((v - yLo) / bw);
      if (b >= 0 && b < nb) counts[b]++; else off++;
    }
    const dens = counts.map(c => c / (N * bw));
    // exact density at time tj
    let eMean: number, eVar: number, pdf: ((y: number) => number) | null = null;
    if (st.model === "ou") {
      const e = Math.exp(-st.theta * tj);
      eMean = st.mu + (x0() - st.mu) * e;
      eVar = st.sigma ** 2 / (2 * st.theta) * (1 - e * e);
      if (eVar > 0) { const m = eMean, v = eVar; pdf = y => npdf(y, m, v); }
    } else {
      const ml = Math.log(x0()) + (st.r - st.sigma ** 2 / 2) * tj, vl = st.sigma ** 2 * tj;
      eMean = x0() * Math.exp(st.r * tj);
      eVar = x0() ** 2 * Math.exp(2 * st.r * tj) * (Math.exp(vl) - 1);
      if (vl > 0) pdf = y => y > 0 ? npdf(Math.log(y), ml, vl) / y : 0;
    }
    const ny = 240, curve: [number, number][] = [];
    let pk = 0;
    if (pdf) for (let i = 0; i <= ny; i++) { const y = yLo + (yHi - yLo) * i / ny, d = pdf(y); curve.push([y, d]); pk = Math.max(pk, d); }
    const dMax = Math.max(1e-9, ...dens, pk) * 1.08;
    const bL = 14, bR = WB - 16;
    const dx = scale(0, dMax, bL, bR);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgB).textContent = `(b) ${N} paths at t = ${fx(tj, 2)}`;
    const ds2 = niceStep(dMax, 3);
    for (let d = 0; d <= dMax + 1e-12; d += ds2) {
      svgEl("line", { x1: dx(d), x2: dx(d), y1: top, y2: bot, class: d === 0 ? "axis" : "grid" }, svgB);
      svgEl("text", { x: dx(d), y: bot + 22, "text-anchor": "middle", class: "sm" }, svgB).textContent = fx(d, ds2 < 0.1 ? 2 : ds2 < 1 ? 1 : 0);
    }
    svgEl("line", { x1: bL, x2: bR, y1: bot, y2: bot, class: "axis" }, svgB);
    svgEl("text", { x: (bL + bR) / 2, y: H - 4, "text-anchor": "middle" }, svgB).textContent = "density";
    for (let b = 0; b < nb; b++) {
      if (!dens[b]) continue;
      const y0 = ay(yLo + (b + 1) * bw), y1 = ay(yLo + b * bw);
      svgEl("rect", { x: bL, y: y0, width: Math.max(0.5, dx(dens[b]) - bL), height: Math.max(0.5, y1 - y0 - 0.6), fill: COL.bar, "fill-opacity": 0.85 }, svgB);
    }
    if (pdf) {
      svgEl("polyline", { points: curve.map(([y, d]) => `${dx(d).toFixed(1)},${ay(y).toFixed(1)}`).join(" "), fill: "none", stroke: COL.exact, "stroke-width": 3 }, svgB);
      const lx = bR - 190;
      svgEl("line", { x1: lx, x2: lx + 24, y1: top + 12, y2: top + 12, stroke: COL.exact, "stroke-width": 3 }, svgB);
      svgEl("text", { x: lx + 30, y: top + 18 }, svgB).textContent = st.model === "ou" ? "exact normal" : "exact lognormal";
    } else {
      svgEl("text", { x: bR, y: top + 18, "text-anchor": "end", class: "sm" }, svgB).textContent = "t = 0: all paths at x₀";
    }
    if (off > 0) {
      const t = svgEl("text", { x: bR, y: top + 46, "text-anchor": "end", "font-weight": 600 }, svgB);
      t.style.fill = COL.warn;
      t.textContent = `${off} off the scale`;
    }

    // ---------- readout ----------
    ro.replaceChildren();
    const line1 = htmlEl("div", "", "", ro);
    const item = (el: HTMLElement, k: string, v: string) => { const s = htmlEl("span", "", "", el); s.style.marginRight = "1.1em"; htmlEl("b", "", k + " ", s); s.append(v); };
    item(line1, "path 1: max |X_h − X| =", fe(maxErr));
    item(line1, `mean over ${P} paths =`, fe(meanMax));
    item(line1, st.model === "ou" ? "θh =" : "rh =", fx((st.model === "ou" ? st.theta : st.r) * h, 3));
    const line2 = htmlEl("div", "", "", ro);
    htmlEl("b", "", `same W, mean max error at h = `, line2);
    const ks = [...hist.keys()].sort((a, b) => a - b).slice(-4);
    line2.append(ks.map((k, i) => {
      const e = hist.get(k)!, prev = i > 0 && ks[i - 1] === k - 1 ? hist.get(k - 1)! : NaN;
      const r = Number.isFinite(prev) && e > 0 && Number.isFinite(e) ? ` (÷${(prev / e).toFixed(2)})` : "";
      return `${fh(Tn / (1 << k))}: ${fe(e)}${r}`;
    }).join(" · "));
    const line3 = htmlEl("div", "", "", ro);
    const sm = nf ? s1 / nf : NaN, sv = nf > 1 ? (s2 - nf * sm * sm) / (nf - 1) : NaN;
    item(line3, `(b) at t = ${fx(tj, 2)}: mean`, `${fe(sm)} (exact ${fe(eMean)}),`);
    let vtxt = `${fe(sv)} (exact ${fe(eVar)}`;
    if (st.model === "ou" && st.j === NSNAP) {
      const th = st.theta * h;
      vtxt += th < 2 ? `; EM long-run σ²/(θ(2−θh)) = ${fe(st.sigma ** 2 / (st.theta * (2 - th)))}` : "; EM long-run: ∞";
    }
    item(line3, "variance", vtxt + ")");
    // mean-square stability
    let warn = "";
    if (st.model === "ou" && st.theta * h >= 2) {
      warn = `θh = ${fx(st.theta * h, 2)} ≥ 2: EM is not mean-square stable (|1 − θh| > 1); paths blow up. Need h < 2/θ = ${fx(2 / st.theta, 3)}.`;
    } else if (st.model === "gbm" && 2 * st.r + st.sigma ** 2 < 0) {
      const a = (1 + st.r * h) ** 2 + st.sigma ** 2 * h + (st.scheme === "mil" ? st.sigma ** 4 * h * h / 2 : 0);
      if (a >= 1) warn = `2r + σ² < 0: exact E[X²] → 0, but ${st.scheme === "em" ? "EM" : "Milstein"} has E[(X_{n+1}/X_n)²] = ${fx(a, 2)} ≥ 1: not mean-square stable.`;
    }
    if (!warn && NUM.some(x => !Number.isFinite(x[n]) || Math.abs(x[n]) > 1e6)) warn = "Some numerical paths overflow.";
    if (warn) htmlEl("div", "sp-warn", "⚠ " + warn, ro);
  }

  syncSliders();
  changed();
}

autoMount(NAME, mount);

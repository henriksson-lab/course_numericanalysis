// Phase-plane explorer for the phase-plane unit: the DYNAMICS of a two-variable ODE.
// (Newton on the same toggle-switch plane is the newton2d demo; this one is about the flow.)
// Models:
//   toggle switch  du/dt = 3 / (1 + v^2) - u,  dv/dt = alpha_v / (1 + u^2) - v
//                  (alpha_v = 3 symmetric, 2.9 biased; slider 2.0 ... 3.5)
//   Sel'kov        du/dt = -u + a v + u^2 v,   dv/dt = b - a v - u^2 v   (a = 0.08, b slider)
// Left: the (u, v) plane with normalised velocity arrows, the nullclines (blue du/dt = 0,
// red dv/dt = 0) and the fixed points (black = stable, white = unstable) with their type
// from the eigenvalues of the Jacobian. Click anywhere to launch a trajectory (RK4,
// h = 0.01, animated); it stops when it is within 1e-3 of a fixed point. Right: u(t) and
// v(t) of the latest trajectory, the fixed points with their eigenvalues, and where the
// latest trajectory ended (Sel'kov: range of u and period after t = 100).
// "basins of the flow" colours every start in [0, 4]^2 by the stable state its trajectory
// reaches (RK4, h = 0.05) and reports the shares; "separatrix" draws the saddle's stable
// manifold, computed by integrating backwards in time from saddle +- 1e-6 v_stable.
// Auto-mounts on <div data-demo="phaseplane">. Optional attributes:
//   data-model="toggle|selkov"  data-av="3"  data-b="0.6"  data-basins="0|1"
//   data-sep="0|1"  data-field="1|0"
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, svgPoint, autoMount } from "../lib/ui";
import { Par, Vec, Model, Fixed, F, fixedPoints, trajectory, separatrix, fate, Traj } from "./math";

const NAME = "phaseplane";
// Okabe-Ito, as in make-phase-plane-figs.R
const COL = {
  f1: "#0072B2", f2: "#D55E00", sep: "#E69F00", field: "#9a9a9a", ink: "#111111",
  hi: "#009E73", mid: "#bdbdbd", lo: "#CC79A7", old: "#555555", start: "#F0E442",
};
const RGB: Record<string, [number, number, number]> = {
  hi: [0x00, 0x9e, 0x73], mid: [0x99, 0x99, 0x99], lo: [0xcc, 0x79, 0xa7],
};

const STYLE = `
[data-demo="${NAME}"].qb-demo { gap: 1em; flex-wrap: wrap; }
[data-demo="${NAME}"].qb-demo .qb-plot { flex: 1 1 420px; max-width: 560px; }
[data-demo="${NAME}"].qb-demo .qb-panel { flex: 1 1 380px; max-width: 600px; gap: .4em; }
[data-demo="${NAME}"].qb-demo svg { width: 100%; height: auto; }
[data-demo="${NAME}"].qb-demo svg.pp-main { touch-action: none; cursor: crosshair; }
[data-demo="${NAME}"].qb-demo svg text { font-size: 19px; }
[data-demo="${NAME}"].qb-demo svg .sm { font-size: 18px; }
[data-demo="${NAME}"].qb-demo .pp-row { display: flex; flex-wrap: wrap; gap: .3em .8em; align-items: center; }
[data-demo="${NAME}"].qb-demo .pp-row .qb-slider { flex: 1 1 300px; grid-template-columns: auto 1fr 3.4em; }
[data-demo="${NAME}"].qb-demo .qb-readout { line-height: 1.4; }
[data-demo="${NAME}"].qb-demo .pp-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"].qb-demo .pp-sw { display: inline-block; width: .9em; height: .9em; border: 1px solid #555; vertical-align: -.1em; margin-right: .2em; }
.reveal [data-demo="${NAME}"].qb-demo { flex-wrap: nowrap; font-size: 18px; gap: .8em; margin: 0; }
.reveal [data-demo="${NAME}"].qb-demo .qb-plot { flex: 0 0 545px; max-width: 545px; }
.reveal [data-demo="${NAME}"].qb-demo .qb-panel { flex: 1 1 560px; max-width: none; gap: .25em; }
.reveal [data-demo="${NAME}"].qb-demo .qb-readout { font-size: 17px; padding: .25em .5em; line-height: 1.3; }
.reveal [data-demo="${NAME}"].qb-demo .qb-btn { font-size: 17px; padding: .12em .5em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- formatting ----
const f3 = (v: number) => v.toFixed(3);
const fl = (v: number) => (v < 0 ? "−" : "") + Math.abs(v).toFixed(3);
const lamStr = (f: Fixed): string => {
  const e = f.e;
  if (e.im > 0) return `λ = ${fl(e.re[0])} ± ${e.im.toFixed(3)}i`;
  return `λ = ${fl(e.re[0])}, ${fl(e.re[1])}`;
};

// ---- plot geometry ----
const W = 560, H = 560;
const L = 58, R = W - 12, T = 12, B = H - 54;
const VIEW: Record<Model, number> = { toggle: 4, selkov: 3 };
const H_STEP = 0.01;
const TMAX: Record<Model, number> = { toggle: 100, selkov: 200 };
const SPEED: Record<Model, number> = { toggle: 10, selkov: 25 };   // time units per second
const MAXTRAJ = 20;
const GRID = 160;          // basin cells per side on [0, 4]^2
// time-course panel
const TW = 560, TH = 190, TL = 52, TR = TW - 24, TT = 10, TB = TH - 44;

interface Run { x0: Vec; tr: Traj; shown: number; els?: [SVGPolylineElement, SVGPolylineElement, SVGCircleElement] }

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const uid = ++mountCount;
  const d = host.dataset;
  const num = (s: string | undefined, def: number) => (s !== undefined && s !== "" && Number.isFinite(Number(s)) ? Number(s) : def);
  const init = {
    model: (d.model === "selkov" ? "selkov" : "toggle") as Model,
    av: Math.min(3.5, Math.max(2, num(d.av, 3))),
    b: Math.min(1.2, Math.max(0.1, num(d.b, 0.6))),
    basins: d.basins === "1",
    sep: d.sep === "1",
    field: d.field !== "0",
  };
  const st = { ...init };
  const par = (): Par => ({ model: st.model, au: 3, av: st.av, a: 0.08, b: st.b });
  let vmax = VIEW[st.model];
  let sx = scale(0, vmax, L, R), sy = scale(0, vmax, B, T);

  // ---- layout ----
  const plot = htmlEl("div", "qb-plot", "", host);
  const panel = htmlEl("div", "qb-panel", "", host);
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, class: "pp-main", role: "img",
    "aria-label": "The (u, v) phase plane with velocity arrows, nullclines, fixed points and trajectories started by clicking." }, plot);
  const defs = svgEl("defs", {}, svg);
  const clip = svgEl("clipPath", { id: `pp-clip-${uid}` }, defs);
  svgEl("rect", { x: L, y: T, width: R - L, height: B - T }, clip);
  const mkF = svgEl("marker", { id: `pp-af-${uid}`, viewBox: "0 0 10 10", refX: 8, refY: 5, markerWidth: 5, markerHeight: 5, orient: "auto" }, defs);
  svgEl("path", { d: "M0,0 L10,5 L0,10 z", fill: COL.field }, mkF);
  const mkT = svgEl("marker", { id: `pp-at-${uid}`, viewBox: "0 0 10 10", refX: 5, refY: 5, markerWidth: 4, markerHeight: 4, orient: "auto" }, defs);
  svgEl("path", { d: "M0,0 L10,5 L0,10 z", fill: COL.ink }, mkT);

  const basinImg = svgEl("image", { x: L, y: T, width: R - L, height: B - T, preserveAspectRatio: "none", "clip-path": `url(#pp-clip-${uid})` }, svg);
  basinImg.style.imageRendering = "pixelated";
  const axes = svgEl("g", {}, svg);
  const fieldG = svgEl("g", { "clip-path": `url(#pp-clip-${uid})` }, svg);
  const ncG = svgEl("g", { "clip-path": `url(#pp-clip-${uid})` }, svg);
  const trajG = svgEl("g", { "clip-path": `url(#pp-clip-${uid})` }, svg);
  const fpG = svgEl("g", {}, svg);

  const tsvg = svgEl("svg", { viewBox: `0 0 ${TW} ${TH}`, role: "img", "aria-label": "u(t) and v(t) of the latest trajectory." });

  // ---- controls ----
  const r1 = htmlEl("div", "pp-row", "", panel);
  const modeBtns = presets(r1, [
    { label: "symmetric toggle", apply: () => setModel("toggle", 3) },
    { label: "biased toggle", apply: () => setModel("toggle", 2.9) },
    { label: "Sel'kov", apply: () => setModel("selkov") },
  ]);
  const r2 = htmlEl("div", "pp-row", "", panel);
  const avS = slider(r2, "promoter αᵥ", 2, 3.5, 0.01, st.av, v => { st.av = v; changedPar(); });
  const bS = slider(r2, "feed b", 0.1, 1.2, 0.01, st.b, v => { st.b = v; changedPar(); });
  const r3 = htmlEl("div", "pp-row", "", panel);
  const fieldC = checkbox(r3, "arrows", st.field, v => { st.field = v; drawStatic(); });
  const basC = checkbox(r3, "basins of the flow", st.basins, v => { st.basins = v; changedBasins(); });
  const sepC = checkbox(r3, "separatrix", st.sep, v => { st.sep = v; drawStatic(); });
  const basLab = basC.parentElement as HTMLElement, sepLab = sepC.parentElement as HTMLElement;
  const r4 = htmlEl("div", "pp-row", "", panel);
  htmlEl("span", "", "start:", r4);
  const startsBox = htmlEl("span", "pp-row", "", r4);
  button(r4, "clear", () => { runs = []; drawRuns(); });
  button(r4, "reset", () => {
    st.model = init.model; st.av = init.av; st.b = init.b; st.basins = init.basins; st.sep = init.sep; st.field = init.field;
    avS.set(st.av); bS.set(st.b); basC.checked = st.basins; sepC.checked = st.sep; fieldC.checked = st.field;
    modeBtns.forEach(b => b.classList.remove("active"));
    runs = []; changedModel();
  });
  panel.appendChild(tsvg);
  const fpRo = htmlEl("div", "qb-readout", "", panel);
  const trRo = htmlEl("div", "qb-readout", "", panel);

  // ---- state ----
  let fps: Fixed[] = fixedPoints(par());
  let runs: Run[] = [];
  let shares = "";
  let anim = 0, lastTs = 0;
  let basinJob = 0;
  let basinKey = "";

  function setModel(m: Model, av?: number): void {
    const was = st.model;
    st.model = m;
    if (av !== undefined) { st.av = av; avS.set(av); }
    if (m === "selkov") { st.b = 0.6; bS.set(0.6); }
    if (m !== was) runs = [];
    changedModel();
  }

  function changedModel(): void {
    vmax = VIEW[st.model];
    sx = scale(0, vmax, L, R); sy = scale(0, vmax, B, T);
    const tg = st.model === "toggle";
    (avS.input.parentElement as HTMLElement).style.display = tg ? "" : "none";
    (bS.input.parentElement as HTMLElement).style.display = tg ? "none" : "";
    basLab.style.display = tg ? "" : "none";
    sepLab.style.display = tg ? "" : "none";
    startsBox.replaceChildren();
    const S: { label: string; x: Vec }[] = tg
      ? [{ label: "(0.5, 0.3)", x: [0.5, 0.3] }, { label: "(0.3, 0.5)", x: [0.3, 0.5] },
         { label: "(2, 2)", x: [2, 2] }, { label: "(2, 2.0001)", x: [2, 2.0001] }]
      : [{ label: "inside (0.61, 1.36)", x: [0.61, 1.36] }, { label: "outside (2.5, 0.2)", x: [2.5, 0.2] }];
    for (const s of S) button(startsBox, s.label, () => launch(s.x));
    if (tg) button(startsBox, "4 × 4 grid", () => {
      for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) launch([0.5 + i, 0.5 + j]);
    });
    drawAxes();
    changedPar();
  }

  function changedPar(): void {
    fps = fixedPoints(par());
    runs = runs.map(r => ({ x0: r.x0, tr: run(r.x0), shown: Infinity }));
    drawStatic();
    changedBasins();
  }

  function run(x0: Vec): Traj { return trajectory(par(), x0, H_STEP, TMAX[st.model], fps); }

  function launch(x0: Vec): void {
    runs.push({ x0, tr: run(x0), shown: 0 });
    if (runs.length > MAXTRAJ) runs.splice(0, runs.length - MAXTRAJ);
    if (!anim) { lastTs = 0; anim = window.requestAnimationFrame(tick); }
    drawRuns();
  }

  function tick(ts: number): void {
    const dt = lastTs ? Math.min(0.1, (ts - lastTs) / 1000) : 0;
    lastTs = ts;
    const adv = dt * SPEED[st.model] / H_STEP;
    let more = false;
    for (const r of runs) {
      const end = lastIdx(r);
      if (r.shown < end) { r.shown = Math.min(end, r.shown + adv); more = r.shown < end || more; updateRun(r); }
    }
    drawTime();
    writeTraj();
    anim = more ? window.requestAnimationFrame(tick) : 0;
  }

  const lastIdx = (r: Run) => (r.tr.settled >= 0 ? r.tr.settled : r.tr.x.length - 1);

  // ---- basins of the flow ----
  function kindOf(f: Fixed): "hi" | "mid" | "lo" {
    if (f.e.saddle) return "mid";
    return f.x[0] > f.x[1] ? "hi" : "lo";
  }

  function changedBasins(): void {
    const show = st.basins && st.model === "toggle";
    if (!show) { basinImg.style.display = "none"; shares = ""; basinJob++; basinKey = ""; writeFp(); return; }
    basinImg.style.display = "";
    const key = `${st.av}`;
    if (key === basinKey) return;
    basinKey = key;
    const job = ++basinJob;
    const p = par(), fs = fps;
    const res = new Int8Array(GRID * GRID);
    let row = 0;
    shares = "<b>basins:</b> computing …";
    writeFp();
    const step = () => {
      if (job !== basinJob) return;
      const t0 = performance.now();
      while (row < GRID && performance.now() - t0 < 30) {
        const v = (GRID - row - 0.5) * 4 / GRID;
        for (let c = 0; c < GRID; c++) res[row * GRID + c] = fate(p, [(c + 0.5) * 4 / GRID, v], fs);
        row++;
      }
      if (row < GRID) { window.requestAnimationFrame(step); return; }
      const cv = document.createElement("canvas");
      cv.width = GRID; cv.height = GRID;
      const ctx = cv.getContext("2d");
      if (!ctx) return;
      const img = ctx.createImageData(GRID, GRID);
      const cnt = new Array(fs.length).fill(0);
      for (let i = 0; i < GRID * GRID; i++) {
        const k = res[i];
        if (k < 0) continue;
        cnt[k]++;
        const [r0, g0, b0] = RGB[kindOf(fs[k])];
        img.data[4 * i] = r0; img.data[4 * i + 1] = g0; img.data[4 * i + 2] = b0; img.data[4 * i + 3] = 120;
      }
      ctx.putImageData(img, 0, 0);
      basinImg.setAttribute("href", cv.toDataURL());
      const tot = GRID * GRID;
      const pc = (c: number) => { const x = 100 * c / tot; return (x > 0 && x < 0.1 ? "< 0.1" : x > 0 && x < 1 ? x.toFixed(1) : x.toFixed(0)) + " %"; };
      // order: U on, middle, V on
      const order = fs.map((f, i) => i).sort((i, j) => fs[j].x[0] - fs[i].x[0]);
      shares = `<b>basins of the flow, starts in [0, 4]²:</b> ` + order.map(i =>
        `<span class="pp-sw" style="background:${COL[kindOf(fs[i])]}"></span>${fs[i].name} ${pc(cnt[i])}`).join("  ");
      writeFp();
    };
    window.requestAnimationFrame(step);
  }

  // ---- pointer: click launches a trajectory ----
  svg.addEventListener("pointerdown", ev => {
    ev.preventDefault();
    const [px, py] = svgPoint(svg, ev);
    if (px < L || px > R || py < T || py > B) return;
    launch([Math.round(sx.invert(px) * 100) / 100, Math.round(sy.invert(py) * 100) / 100]);
  });

  // ---- drawing ----
  function drawAxes(): void {
    axes.replaceChildren();
    for (let t = 0; t <= vmax; t += 0.5) {
      const major = Number.isInteger(t);
      svgEl("line", { x1: sx(t), x2: sx(t), y1: T, y2: B, class: "grid", "stroke-opacity": major ? 1 : 0.5 }, axes);
      svgEl("line", { x1: L, x2: R, y1: sy(t), y2: sy(t), class: "grid", "stroke-opacity": major ? 1 : 0.5 }, axes);
      if (!major && vmax > 3) continue;
      svgEl("text", { x: sx(t), y: B + 22, "text-anchor": "middle", class: "sm" }, axes).textContent = String(t);
      svgEl("text", { x: L - 8, y: sy(t) + 6, "text-anchor": "end", class: "sm" }, axes).textContent = String(t);
    }
    svgEl("rect", { x: L, y: T, width: R - L, height: B - T, fill: "none", class: "axis" }, axes);
    const tg = st.model === "toggle";
    svgEl("text", { x: (L + R) / 2, y: H - 8, "text-anchor": "middle" }, axes).textContent = tg ? "u (protein U)" : "u (ADP)";
    const yl = svgEl("text", { x: 16, y: (T + B) / 2, "text-anchor": "middle", transform: `rotate(-90 16 ${(T + B) / 2})` }, axes);
    yl.textContent = tg ? "v (protein V)" : "v (fructose-6-phosphate)";
  }

  const pts = (xs: Vec[], every = 1): string => {
    const out: string[] = [];
    for (let i = 0; i < xs.length; i += every) out.push(`${sx(xs[i][0]).toFixed(1)},${sy(xs[i][1]).toFixed(1)}`);
    if ((xs.length - 1) % every) { const l = xs[xs.length - 1]; out.push(`${sx(l[0]).toFixed(1)},${sy(l[1]).toFixed(1)}`); }
    return out.join(" ");
  };

  function drawStatic(): void {
    const p = par();
    fieldG.replaceChildren(); ncG.replaceChildren(); fpG.replaceChildren();
    if (st.field) {
      const n = 15, len = 0.42 * (R - L) / n;
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
        const u = (i + 0.5) * vmax / n, v = (j + 0.5) * vmax / n;
        const f = F(p, [u, v]);
        const dx = f[0], dy = -f[1], m = Math.hypot(dx, dy);
        if (!(m > 1e-9)) continue;
        const cx = sx(u), cy = sy(v), ex = dx / m * len, ey = dy / m * len;
        svgEl("line", { x1: cx - ex / 2, y1: cy - ey / 2, x2: cx + ex / 2, y2: cy + ey / 2, stroke: COL.field, "stroke-width": 1.6,
          "marker-end": `url(#pp-af-${uid})` }, fieldG);
      }
    }
    // nullclines
    const S = Array.from({ length: 401 }, (_, i) => i * vmax / 400);
    let n1: Vec[], n2: Vec[];
    if (st.model === "toggle") {
      n1 = S.map(v => [p.au / (1 + v * v), v]);
      n2 = S.map(u => [u, p.av / (1 + u * u)]);
    } else {
      const S2 = Array.from({ length: 801 }, (_, i) => i * vmax / 800);
      n1 = S2.map(u => [u, u / (p.a + u * u)]);
      n2 = S2.map(u => [u, Math.min(50, p.b / (p.a + u * u))]);
    }
    svgEl("polyline", { points: pts(n1), fill: "none", stroke: COL.f1, "stroke-width": 3.5 }, ncG);
    svgEl("polyline", { points: pts(n2), fill: "none", stroke: COL.f2, "stroke-width": 3.5 }, ncG);
    // separatrix
    if (st.sep && st.model === "toggle") {
      for (const f of fps) if (f.e.saddle) {
        const s = separatrix(p, f, [-0.5, vmax + 0.5, -0.5, vmax + 0.5]);
        svgEl("polyline", { points: pts(s, 4), fill: "none", stroke: "#fff", "stroke-width": 7, "stroke-opacity": 0.7 }, ncG);
        svgEl("polyline", { points: pts(s, 4), fill: "none", stroke: COL.sep, "stroke-width": 4 }, ncG);
      }
    }
    // fixed points with their type
    for (const f of fps) {
      const cx = sx(f.x[0]), cy = sy(f.x[1]);
      svgEl("circle", { cx, cy, r: 9, fill: f.e.stable ? COL.ink : "#fff", stroke: COL.ink, "stroke-width": 2.5 }, fpG);
      const right = cx < (L + R) / 2;
      const t = svgEl("text", { x: cx + (right ? 14 : -14), y: cy + (f.e.saddle ? 26 : -12), "text-anchor": right ? "start" : "end",
        "font-weight": 700, class: "sm", stroke: "#fff", "stroke-width": 5, "paint-order": "stroke" }, fpG);
      t.textContent = f.e.type;
    }
    writeFp();
    drawRuns();
  }

  function writeFp(): void {
    const lines = fps.map(f => {
      const nm = st.model === "toggle" ? f.name.padEnd(8) + " " : "";
      return `<b>${nm}</b>(${f3(f.x[0])}, ${f3(f.x[1])})  ${lamStr(f)}  <b>${f.e.type}</b>`;
    });
    if (st.model === "selkov") {
      const f = fps[0];
      lines.push(`trace ${fl(f.e.tr)}, det ${f3(f.e.det)}`);
    }
    fpRo.innerHTML = lines.join("\n") + (shares ? "\n" + shares : "");
  }

  function drawRuns(): void {
    trajG.replaceChildren();
    runs.forEach((r, i) => {
      const latest = i === runs.length - 1;
      const col = latest ? COL.ink : COL.old;
      const halo = svgEl("polyline", { fill: "none", stroke: "#fff", "stroke-width": latest ? 6 : 4, "stroke-opacity": 0.6 }, trajG);
      const line = svgEl("polyline", { fill: "none", stroke: col, "stroke-width": latest ? 3 : 1.8 }, trajG);
      const s = r.x0;
      svgEl("rect", { x: sx(s[0]) - 6, y: sy(s[1]) - 6, width: 12, height: 12, fill: COL.start, stroke: COL.ink, "stroke-width": 1.5 }, trajG);
      const head = svgEl("circle", { r: latest ? 5 : 3.5, fill: col }, trajG);
      r.els = [halo, line, head];
      updateRun(r);
    });
    drawTime();
    writeTraj();
  }

  /** Redraw one trajectory up to its shown index (called every animation frame). */
  function updateRun(r: Run): void {
    if (!r.els) return;
    const k = Math.floor(Math.min(r.shown, lastIdx(r)));
    const xs = r.tr.x.slice(0, k + 1);
    const every = st.model === "toggle" ? 3 : 5;
    const ps = xs.length > 1 ? pts(xs, every) : "";
    r.els[0].setAttribute("points", ps);
    r.els[1].setAttribute("points", ps);
    const h = xs[xs.length - 1];
    r.els[2].setAttribute("cx", String(sx(h[0])));
    r.els[2].setAttribute("cy", String(sy(h[1])));
  }

  function drawTime(): void {
    tsvg.replaceChildren();
    const r = runs[runs.length - 1];
    const tEnd = r ? Math.max(10, Math.ceil(r.tr.t[lastIdx(r)] / 10) * 10) : TMAX[st.model];
    const ymax = st.model === "toggle" ? 4 : 3;
    const tx = scale(0, tEnd, TL, TR), ty = scale(0, ymax, TB, TT);
    const stepT = tEnd <= 30 ? 5 : tEnd <= 100 ? 10 : 20;
    for (let t = 0; t <= tEnd + 1e-9; t += stepT) {
      svgEl("line", { x1: tx(t), x2: tx(t), y1: TT, y2: TB, class: "grid" }, tsvg);
      svgEl("text", { x: tx(t), y: TB + 20, "text-anchor": "middle", class: "sm" }, tsvg).textContent = String(t);
    }
    for (let y = 0; y <= ymax; y++) {
      svgEl("line", { x1: TL, x2: TR, y1: ty(y), y2: ty(y), class: "grid" }, tsvg);
      svgEl("text", { x: TL - 8, y: ty(y) + 6, "text-anchor": "end", class: "sm" }, tsvg).textContent = String(y);
    }
    svgEl("rect", { x: TL, y: TT, width: TR - TL, height: TB - TT, fill: "none", class: "axis" }, tsvg);
    svgEl("text", { x: (TL + TR) / 2, y: TH - 4, "text-anchor": "middle", class: "sm" }, tsvg).textContent = "time t";
    if (!r) {
      svgEl("text", { x: (TL + TR) / 2, y: (TT + TB) / 2 + 6, "text-anchor": "middle", class: "sm" }, tsvg).textContent = "click in the plane to start a trajectory";
      return;
    }
    const k = Math.floor(Math.min(r.shown, lastIdx(r)));
    const every = Math.max(1, Math.round(k / 800));
    const line = (c: 0 | 1) => {
      const out: string[] = [];
      for (let i = 0; i <= k; i += every) out.push(`${tx(r.tr.t[i]).toFixed(1)},${ty(Math.min(ymax * 1.05, r.tr.x[i][c])).toFixed(1)}`);
      return out.join(" ");
    };
    svgEl("polyline", { points: line(0), fill: "none", stroke: COL.f1, "stroke-width": 2.5 }, tsvg);
    svgEl("polyline", { points: line(1), fill: "none", stroke: COL.f2, "stroke-width": 2.5, "stroke-dasharray": "7 4" }, tsvg);
    const lg = (x: number, txt: string, col: string, dash: string) => {
      svgEl("line", { x1: x, x2: x + 28, y1: TT + 16, y2: TT + 16, stroke: col, "stroke-width": 3, "stroke-dasharray": dash }, tsvg);
      const t = svgEl("text", { x: x + 34, y: TT + 22, class: "sm", stroke: "#fff", "stroke-width": 4, "paint-order": "stroke" }, tsvg);
      t.textContent = txt;
    };
    lg(TR - 170, "u(t)", COL.f1, "");
    lg(TR - 82, "v(t)", COL.f2, "7 4");
  }

  function writeTraj(): void {
    const r = runs[runs.length - 1];
    if (!r) { trRo.innerHTML = "<b>trajectory:</b> click in the plane, or a start button"; return; }
    const tr = r.tr, li = lastIdx(r), k = Math.floor(Math.min(r.shown, li));
    const x = tr.x[k];
    let s = `<b>start</b> (${r.x0[0]}, ${r.x0[1]})   <b>t = ${tr.t[k].toFixed(1)}</b>   (u, v) = (${f3(x[0])}, ${f3(x[1])})`;
    if (k < li) { trRo.innerHTML = s + "\n…"; return; }
    if (tr.settled >= 0) {
      const f = fps[tr.end];
      const nm = st.model === "toggle" ? f.name : "the fixed point";
      if (f.e.stable) s += `\n<b>ends:</b> ${nm}, within 10⁻³ at t = ${tr.t[tr.settled].toFixed(1)}`;
      else s += `\n<b>stuck:</b> within 10⁻³ of ${nm} from t = ${tr.t[tr.settled].toFixed(1)} to ${TMAX[st.model]}` +
        `\n<span class="pp-warn">an UNSTABLE state (${f.e.type}): the start lies on (or very near) its stable manifold</span>`;
    } else if (st.model === "toggle") {
      s += `\n<b>not settled</b> by t = ${TMAX.toggle}`;
    } else {
      const b = st.b, t0 = TMAX.selkov / 2;
      let umin = Infinity, umax = -Infinity;
      const ups: number[] = [];
      for (let i = 1; i < tr.x.length; i++) {
        if (tr.t[i] <= t0) continue;
        const u = tr.x[i][0];
        umin = Math.min(umin, u); umax = Math.max(umax, u);
        const u0 = tr.x[i - 1][0];
        if (u0 < b && u >= b) ups.push(tr.t[i - 1] + (b - u0) / (u - u0) * H_STEP);
      }
      s += `\n<b>after t = ${t0}:</b> u from ${f3(umin)} to ${f3(umax)} (range ${(umax - umin).toPrecision(2)})`;
      if (ups.length > 1) s += `\n<b>period</b> ${((ups[ups.length - 1] - ups[0]) / (ups.length - 1)).toFixed(2)} (upward crossings of u = b)`;
    }
    trRo.innerHTML = s;
  }

  changedModel();
}

autoMount(NAME, mount);

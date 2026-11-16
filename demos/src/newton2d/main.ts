// Newton for a system, on the toggle switch of the newton-systems unit:
//   du/dt = alpha / (1 + v^n) - u,   dv/dt = alpha / (1 + u^n) - v   (alpha = 3, n = 2)
// The (u, v) plane with the two nullclines F1 = 0 (blue) and F2 = 0 (red) and the steady
// states (their crossings). Click or drag to set the start x0; the Newton path is drawn.
// Step through it: at x_k the two dashed lines are the *linearised* nullclines,
// F_i(x_k) + J_i(x_k) . delta = 0; they cross at x_{k+1} = x_k + solve(J, -F).
// The readout gives F, J, det J and the step. "Basins" colours every start by the steady
// state Newton reaches (green U high, grey middle, purple V high, black: no convergence),
// computed here on a grid, and reports the shares on the square [0, 4]^2.
// "det J = 0" draws the curve where the Jacobian is singular.
// Auto-mounts on <div data-demo="newton2d">. Optional attributes:
//   data-alpha="3"  data-n="2"  data-start="3 0.5"  data-basins="0|1"  data-detcurve="0|1"
import { svgEl, htmlEl, slider, button, presets, checkbox, scale, svgPoint, numList, autoMount } from "../lib/ui";
import { Par, Vec, F, step, newton, roots, basinOf, norm, Root, Path, MAXIT } from "./math";

const NAME = "newton2d";
// Okabe-Ito, as in make-newton-systems-figs.R
const COL = {
  f1: "#0072B2", f2: "#D55E00", start: "#E69F00", ink: "#111111",
  hi: "#009E73", mid: "#d9d9d9", lo: "#CC79A7", fail: "#222222",
};
const RGB: Record<string, [number, number, number]> = {
  hi: [0x00, 0x9e, 0x73], mid: [0xd9, 0xd9, 0xd9], lo: [0xcc, 0x79, 0xa7], fail: [0x22, 0x22, 0x22],
};
const KIND_LABEL = { hi: "U high", mid: "middle", lo: "V high" };

const STYLE = `
[data-demo="${NAME}"].qb-demo { gap: 1em; flex-wrap: wrap; }
[data-demo="${NAME}"].qb-demo .qb-plot { flex: 1 1 420px; max-width: 560px; }
[data-demo="${NAME}"].qb-demo .qb-panel { flex: 1 1 380px; max-width: 580px; gap: .45em; }
[data-demo="${NAME}"].qb-demo svg { width: 100%; height: auto; touch-action: none; cursor: crosshair; }
[data-demo="${NAME}"].qb-demo svg text { font-size: 19px; }
[data-demo="${NAME}"].qb-demo svg .sm { font-size: 18px; }
[data-demo="${NAME}"].qb-demo .nw-row { display: flex; flex-wrap: wrap; gap: .35em .8em; align-items: center; }
[data-demo="${NAME}"].qb-demo .qb-readout { line-height: 1.4; }
[data-demo="${NAME}"].qb-demo .qb-readout b { font-weight: 700; }
[data-demo="${NAME}"].qb-demo .nw-warn { color: #b00020; font-weight: 600; }
[data-demo="${NAME}"].qb-demo .nw-sw { display: inline-block; width: .9em; height: .9em; border: 1px solid #555; vertical-align: -.1em; margin-right: .2em; }
[data-demo="${NAME}"].qb-demo .nw-legend { font-size: .85em; display: flex; flex-wrap: wrap; gap: .1em 1em; }
[data-demo="${NAME}"].qb-demo .nw-legend .ln { display: inline-block; width: 1.6em; height: 0; border-top: 3px solid; vertical-align: .3em; margin-right: .25em; }
[data-demo="${NAME}"].qb-demo .nw-legend .ln.dash { border-top-style: dashed; }
.reveal [data-demo="${NAME}"].qb-demo { flex-wrap: nowrap; font-size: 19px; gap: .8em; margin: 0; }
.reveal [data-demo="${NAME}"].qb-demo .qb-plot { flex: 0 0 540px; max-width: 540px; }
.reveal [data-demo="${NAME}"].qb-demo .qb-panel { flex: 1 1 560px; max-width: none; gap: .3em; }
.reveal [data-demo="${NAME}"].qb-demo .qb-readout { font-size: 18px; padding: .3em .5em; }
.reveal [data-demo="${NAME}"].qb-demo .qb-btn { font-size: 17px; padding: .15em .55em; }
.reveal [data-demo="${NAME}"].qb-demo .nw-legend { font-size: 17px; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---- formatting ----
const fx = (v: number, d = 6): string => (!Number.isFinite(v) ? String(v) : (v < 0 ? "" : " ") + v.toFixed(d));
const fe = (v: number): string => {
  if (!Number.isFinite(v)) return String(v);
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e-2 && a < 1e4) return v.toPrecision(3);
  return v.toExponential(1).replace("e+", "e");
};
const fv = (v: number): string => (!Number.isFinite(v) ? String(v) : Math.abs(v) >= 1e-4 || v === 0 ? fx(v) : (v < 0 ? "" : " ") + v.toExponential(2));
const sub = (k: number): string => String(k).split("").map(c => "₀₁₂₃₄₅₆₇₈₉"[Number(c)]).join("");

// ---- plot geometry ----
const W = 560, H = 560;
const L = 62, R = W - 12, T = 14, B = H - 56;
const XMIN = -0.5, XMAX = 4.5;
const sx = scale(XMIN, XMAX, L, R), sy = scale(XMIN, XMAX, B, T);
const GRID = 250;          // basin pixels per side over the visible square
const SHARE = 200;         // grid for the shares on [0, 4]^2 (cell centres)

let mountCount = 0;

function mount(host: HTMLElement): void {
  injectStyle();
  const uid = ++mountCount;
  const d = host.dataset;
  const num = (s: string | undefined, def: number) => (s !== undefined && s !== "" && Number.isFinite(Number(s)) ? Number(s) : def);
  const s0 = numList(d.start, [3, 0.5]);
  const init = {
    alpha: Math.min(5, Math.max(0.5, num(d.alpha, 3))),
    n: Math.min(4, Math.max(1, num(d.n, 2))),
    x0: [s0[0] ?? 3, s0[1] ?? 0.5] as Vec,
    basins: d.basins === "1",
    detc: d.detcurve === "1",
  };
  const st = { alpha: init.alpha, n: init.n, x0: init.x0, k: 0, basins: init.basins, detc: init.detc, lin: true };
  const par = (): Par => ({ alpha: st.alpha, n: st.n });

  // ---- layout ----
  const plot = htmlEl("div", "qb-plot", "", host);
  const panel = htmlEl("div", "qb-panel", "", host);
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
    "aria-label": "The (u, v) plane with the nullclines of the toggle switch, its steady states and the Newton path from the chosen start." }, plot);
  const defs = svgEl("defs", {}, svg);
  const clip = svgEl("clipPath", { id: `nw-clip-${uid}` }, defs);
  svgEl("rect", { x: L, y: T, width: R - L, height: B - T }, clip);
  const mk = svgEl("marker", { id: `nw-arr-${uid}`, viewBox: "0 0 10 10", refX: 9, refY: 5, markerWidth: 7, markerHeight: 7, orient: "auto-start-reverse" }, defs);
  svgEl("path", { d: "M0,0 L10,5 L0,10 z", fill: COL.ink }, mk);

  const basinImg = svgEl("image", { x: sx(XMIN), y: sy(XMAX), width: sx(XMAX) - sx(XMIN), height: sy(XMIN) - sy(XMAX),
    preserveAspectRatio: "none", "clip-path": `url(#nw-clip-${uid})` }, svg);
  basinImg.style.imageRendering = "pixelated";
  const axes = svgEl("g", {}, svg);
  const layer = svgEl("g", { "clip-path": `url(#nw-clip-${uid})` }, svg);
  const over = svgEl("g", {}, svg);
  drawAxes(axes);

  // ---- controls ----
  const r1 = htmlEl("div", "nw-row", "", panel);
  const aS = slider(r1, "promoter α", 0.5, 5, 0.05, st.alpha, v => { st.alpha = v; changedPar(); });
  aS.input.parentElement!.style.flex = "1 1 240px";
  const nS = slider(r1, "Hill n", 1, 4, 0.1, st.n, v => { st.n = v; changedPar(); });
  nS.input.parentElement!.style.flex = "1 1 200px";

  const r2 = htmlEl("div", "nw-row", "", panel);
  htmlEl("span", "", "step:", r2);
  button(r2, "x₀", () => { st.k = 0; draw(); });
  button(r2, "− 1", () => { st.k = Math.max(0, Math.min(st.k, path.xs.length - 1) - 1); draw(); });
  button(r2, "+ 1", () => { st.k = Math.min(path.xs.length - 1, st.k + 1); draw(); });
  button(r2, "all", () => { st.k = Infinity; draw(); });

  const r3 = htmlEl("div", "nw-row", "", panel);
  const linC = checkbox(r3, "linearisation", st.lin, v => { st.lin = v; draw(); });
  const basC = checkbox(r3, "basins", st.basins, v => { st.basins = v; changedBasins(); });
  const detC = checkbox(r3, "det J = 0", st.detc, v => { st.detc = v; draw(); });

  const r4 = htmlEl("div", "nw-row", "", panel);
  htmlEl("span", "", "start:", r4);
  const setStart = (x: Vec, keepAll = false) => {
    st.x0 = x; path = newton(par(), x);
    if (!(keepAll && st.k === Infinity)) st.k = 0;
    draw();
  };
  const pre = presets(r4, [
    { label: "(3, 0.5)", apply: () => setStart([3, 0.5]) },
    { label: "(1, 0.2)", apply: () => setStart([1, 0.2]) },
    { label: "(1, 1)", apply: () => setStart([1, 1]) },
    { label: "(1.9, 0.5)", apply: () => setStart([1.9, 0.5]) },
    { label: "(2, 0.55)", apply: () => setStart([2, 0.55]) },
  ]);
  pre.forEach(b => b.addEventListener("click", () => window.setTimeout(() => b.classList.remove("active"), 600)));
  button(r4, "reset", () => {
    st.alpha = init.alpha; st.n = init.n; aS.set(st.alpha); nS.set(st.n);
    st.basins = init.basins; st.detc = init.detc; st.lin = true;
    basC.checked = st.basins; detC.checked = st.detc; linC.checked = true;
    rs = roots(par()); basinKey = ""; setStart(init.x0); changedBasins();
  });

  const ro = htmlEl("div", "qb-readout", "", panel);
  const fin = htmlEl("div", "qb-readout", "", panel);
  const legend = htmlEl("div", "nw-legend", "", panel);
  legend.innerHTML =
    `<span><span class="ln" style="border-color:${COL.f1}"></span>F₁ = 0</span>` +
    `<span><span class="ln" style="border-color:${COL.f2}"></span>F₂ = 0</span>` +
    `<span><span class="ln dash" style="border-color:${COL.f1}"></span><span class="ln dash" style="border-color:${COL.f2}"></span>linearised at x<sub>k</sub>: they cross at x<sub>k+1</sub></span>` +
    `<span>click or drag in the plot to move the start</span>`;

  // ---- state ----
  let rs: Root[] = roots(par());
  let path: Path = newton(par(), st.x0);
  let basinKey = "";
  let shares = "";
  let pending = 0;

  function changedPar(): void {
    rs = roots(par());
    path = newton(par(), st.x0);
    if (st.k !== Infinity) st.k = Math.min(st.k, path.xs.length - 1);
    draw();
    changedBasins();
  }
  function changedBasins(): void {
    if (!st.basins) { basinImg.removeAttribute("href"); basinImg.style.display = "none"; shares = ""; draw(); return; }
    basinImg.style.display = "";
    if (pending) return;
    pending = window.requestAnimationFrame(() => { pending = 0; computeBasins(); draw(); });
  }

  function computeBasins(): void {
    const key = `${st.alpha}|${st.n}`;
    if (key === basinKey) return;
    basinKey = key;
    const p = par();
    const cv = document.createElement("canvas");
    cv.width = GRID; cv.height = GRID;
    const ctx = cv.getContext("2d");
    if (!ctx) return;
    const img = ctx.createImageData(GRID, GRID);
    const colOf = (i: number) => RGB[i < 0 ? "fail" : rs[i].kind];
    for (let r = 0; r < GRID; r++) {
      const v = XMAX - (r + 0.5) * (XMAX - XMIN) / GRID;
      for (let c = 0; c < GRID; c++) {
        const u = XMIN + (c + 0.5) * (XMAX - XMIN) / GRID;
        const [R0, G0, B0] = colOf(basinOf(p, rs, [u, v]));
        const o = 4 * (r * GRID + c);
        img.data[o] = R0; img.data[o + 1] = G0; img.data[o + 2] = B0; img.data[o + 3] = 150;
      }
    }
    ctx.putImageData(img, 0, 0);
    basinImg.setAttribute("href", cv.toDataURL());
    // shares on [0, 4]^2
    const cnt = new Array(rs.length + 1).fill(0);
    for (let i = 0; i < SHARE; i++) for (let j = 0; j < SHARE; j++) {
      const b = basinOf(p, rs, [(i + 0.5) * 4 / SHARE, (j + 0.5) * 4 / SHARE]);
      cnt[b < 0 ? rs.length : b]++;
    }
    const tot = SHARE * SHARE, pc = (c: number) => (100 * c / tot).toFixed(0) + " %";
    const parts = rs.map((r, i) => `<span class="nw-sw" style="background:${COL[r.kind]}"></span>${rs.length === 1 ? "the only state" : KIND_LABEL[r.kind]} ${pc(cnt[i])}`);
    parts.push(`<span class="nw-sw" style="background:${COL.fail}"></span>no convergence ${pc(cnt[rs.length])}`);
    shares = `<b>basins, starts in [0, 4]²:</b> ` + parts.join(" &nbsp;");
  }

  // ---- pointer: click or drag moves the start ----
  const toWorld = (ev: PointerEvent): Vec => {
    const [px, py] = svgPoint(svg, ev);
    const c = (v: number) => Math.min(XMAX, Math.max(XMIN, v));
    return [Math.round(c(sx.invert(px)) * 100) / 100, Math.round(c(sy.invert(py)) * 100) / 100];
  };
  svg.addEventListener("pointerdown", ev => {
    ev.preventDefault();
    svg.setPointerCapture(ev.pointerId);
    setStart(toWorld(ev), true);
    const move = (e: PointerEvent) => setStart(toWorld(e), true);
    const up = () => {
      svg.removeEventListener("pointermove", move);
      svg.removeEventListener("pointerup", up);
      svg.removeEventListener("pointercancel", up);
    };
    svg.addEventListener("pointermove", move);
    svg.addEventListener("pointerup", up);
    svg.addEventListener("pointercancel", up);
  });

  // ---- drawing ----
  function drawAxes(g: SVGGElement): void {
    for (let t = 0; t <= 4; t++) {
      svgEl("line", { x1: sx(t), x2: sx(t), y1: T, y2: B, class: "grid" }, g);
      svgEl("line", { x1: L, x2: R, y1: sy(t), y2: sy(t), class: "grid" }, g);
      svgEl("text", { x: sx(t), y: B + 22, "text-anchor": "middle", class: "sm" }, g).textContent = String(t);
      svgEl("text", { x: L - 8, y: sy(t) + 6, "text-anchor": "end", class: "sm" }, g).textContent = String(t);
    }
    svgEl("rect", { x: L, y: T, width: R - L, height: B - T, fill: "none", class: "axis" }, g);
    svgEl("text", { x: (L + R) / 2, y: H - 8, "text-anchor": "middle" }, g).textContent = "u (protein U)";
    const yl = svgEl("text", { x: 18, y: (T + B) / 2, "text-anchor": "middle", transform: `rotate(-90 18 ${(T + B) / 2})` }, g);
    yl.textContent = "v (protein V)";
  }

  function polyline(pts: Vec[], attrs: Record<string, string | number>, parent: Element): void {
    const segs: string[][] = [[]];
    for (const [u, v] of pts) {
      if (!Number.isFinite(u) || !Number.isFinite(v)) { segs.push([]); continue; }
      segs[segs.length - 1].push(`${sx(u).toFixed(1)},${sy(v).toFixed(1)}`);
    }
    for (const s of segs) if (s.length > 1) svgEl("polyline", { points: s.join(" "), fill: "none", ...attrs }, parent);
  }

  /** Line {x : a . (x - x0) = -f} across the view, i.e. f + a . delta = 0. */
  function linLine(a: [number, number], f: number, x0: Vec, col: string): void {
    // points x = x0 + delta with a0 d0 + a1 d1 = -f
    const pts: Vec[] = [];
    const lo = XMIN - 50, hi = XMAX + 50;
    if (Math.abs(a[1]) > Math.abs(a[0])) {
      for (const u of [lo, hi]) pts.push([u, x0[1] + (-f - a[0] * (u - x0[0])) / a[1]]);
    } else {
      for (const v of [lo, hi]) pts.push([x0[0] + (-f - a[1] * (v - x0[1])) / a[0], v]);
    }
    if (pts.some(p => !Number.isFinite(p[0]) || !Number.isFinite(p[1]))) return;
    svgEl("line", { x1: sx(pts[0][0]), y1: sy(pts[0][1]), x2: sx(pts[1][0]), y2: sy(pts[1][1]),
      stroke: col, "stroke-width": 2.5, "stroke-dasharray": "9 6" }, layer);
  }

  function detCurve(): void {
    // marching squares on det J(u, v) = 1 - b(v) c(u)
    const p = par(), Nn = 160, h = (XMAX - XMIN) / Nn;
    const a = p.alpha, n = p.n;
    const bf = (v: number) => -a * n * v ** (n - 1) / (1 + v ** n) ** 2;
    const B_ = Array.from({ length: Nn + 1 }, (_, i) => bf(XMIN + i * h));
    const D = (i: number, j: number) => 1 - B_[j] * B_[i];   // c(u) = b(u) by symmetry
    let dstr = "";
    for (let i = 0; i < Nn; i++) for (let j = 0; j < Nn; j++) {
      const c = [D(i, j), D(i + 1, j), D(i + 1, j + 1), D(i, j + 1)];
      if (c.some(x => !Number.isFinite(x))) continue;
      const P: Vec[] = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]];
      const cr: Vec[] = [];
      for (let e = 0; e < 4; e++) {
        const c0 = c[e], c1 = c[(e + 1) % 4];
        if ((c0 < 0) !== (c1 < 0)) {
          const t = c0 / (c0 - c1), p0 = P[e], p1 = P[(e + 1) % 4];
          cr.push([XMIN + (p0[0] + t * (p1[0] - p0[0])) * h, XMIN + (p0[1] + t * (p1[1] - p0[1])) * h]);
        }
      }
      for (let q = 0; q + 1 < cr.length; q += 2)
        dstr += `M${sx(cr[q][0]).toFixed(1)},${sy(cr[q][1]).toFixed(1)}L${sx(cr[q + 1][0]).toFixed(1)},${sy(cr[q + 1][1]).toFixed(1)}`;
    }
    if (dstr) svgEl("path", { d: dstr, fill: "none", stroke: "#000", "stroke-width": 2, "stroke-dasharray": "2 4", "stroke-linecap": "round" }, layer);
  }

  function draw(): void {
    const p = par();
    layer.replaceChildren();
    over.replaceChildren();
    // nullclines: F1 = 0 is u = alpha/(1 + v^n), F2 = 0 is v = alpha/(1 + u^n)
    const S = Array.from({ length: 401 }, (_, i) => XMIN + i * (XMAX - XMIN) / 400);
    polyline(S.map(v => [p.alpha / (1 + v ** p.n), v] as Vec), { stroke: COL.f1, "stroke-width": 3.5 }, layer);
    polyline(S.map(u => [u, p.alpha / (1 + u ** p.n)] as Vec), { stroke: COL.f2, "stroke-width": 3.5 }, layer);
    if (st.detc) detCurve();

    const xs = path.xs, K = xs.length - 1, k = Math.min(st.k, K);
    // whole path, faint
    polyline(xs, { stroke: "#fff", "stroke-width": 5, "stroke-opacity": 0.6 }, layer);
    polyline(xs, { stroke: "#555", "stroke-width": 1.5, "stroke-dasharray": "4 4" }, layer);
    // linearisation at x_k
    let s: ReturnType<typeof step> | null = null;
    if (k < K || path.status !== "converged") {
      const xk = xs[k];
      s = step(p, xk);
      if (st.lin && Number.isFinite(s.d[0]) && Number.isFinite(s.d[1])) {
        linLine([s.j[0], s.j[1]], s.f[0], xk, COL.f1);
        linLine([s.j[2], s.j[3]], s.f[1], xk, COL.f2);
      }
    }
    // path up to x_k (and the step to x_{k+1})
    const kk = s && k < K ? k + 1 : k;
    for (let i = 0; i < kk; i++) {
      const a = xs[i], b = xs[i + 1];
      if (![...a, ...b].every(Number.isFinite)) break;
      const len = Math.hypot(sx(b[0]) - sx(a[0]), sy(b[1]) - sy(a[1]));
      svgEl("line", { x1: sx(a[0]), y1: sy(a[1]), x2: sx(b[0]), y2: sy(b[1]), stroke: "#fff", "stroke-width": 6, "stroke-opacity": 0.8 }, layer);
      svgEl("line", { x1: sx(a[0]), y1: sy(a[1]), x2: sx(b[0]), y2: sy(b[1]), stroke: COL.ink, "stroke-width": 2.5,
        "marker-end": len > 14 ? `url(#nw-arr-${uid})` : "" }, layer);
    }
    // steady states
    for (const r of rs) svgEl("circle", { cx: sx(r.x[0]), cy: sy(r.x[1]), r: 9, fill: r.kind === "mid" ? "#fff" : COL[r.kind], stroke: COL.ink, "stroke-width": 2 }, layer);
    // iterates with labels
    let lastLab: [number, number] | null = null;
    for (let i = 0; i <= kk; i++) {
      const x = xs[i];
      if (!x.every(Number.isFinite)) break;
      const cx = sx(x[0]), cy = sy(x[1]);
      if (i > 0) svgEl("circle", { cx, cy, r: i === kk && kk > k ? 7 : 5, fill: i === kk && kk > k ? "#fff" : COL.ink, stroke: COL.ink, "stroke-width": 2 }, layer);
      if (i > 0 && (!lastLab || Math.hypot(cx - lastLab[0], cy - lastLab[1]) > 26)) {
        const t = svgEl("text", { x: cx + 10, y: cy - 10, "font-weight": 700, stroke: "#fff", "stroke-width": 4, "paint-order": "stroke" }, layer);
        t.textContent = "x" + sub(i);
        lastLab = [cx, cy];
      }
    }
    // start handle (always on top)
    const [u0, v0] = st.x0;
    svgEl("circle", { cx: sx(u0), cy: sy(v0), r: 11, fill: COL.start, stroke: COL.ink, "stroke-width": 2.5 }, over);
    const t0 = svgEl("text", { x: sx(u0) + 14, y: sy(v0) + 22, "font-weight": 700, stroke: "#fff", "stroke-width": 4, "paint-order": "stroke" }, over);
    t0.textContent = "x₀";
    // off-screen next iterate
    if (kk > k) {
      const nx = xs[kk];
      const off = nx.every(Number.isFinite) && (nx[0] < XMIN || nx[0] > XMAX || nx[1] < XMIN || nx[1] > XMAX);
      if (off) {
        const t = svgEl("text", { x: (L + R) / 2, y: T + 26, "text-anchor": "middle", "font-weight": 700, stroke: "#fff", "stroke-width": 4, "paint-order": "stroke" }, over);
        t.style.fill = "#b00020";
        t.textContent = `x${sub(kk)} = (${nx[0].toFixed(2)}, ${nx[1].toFixed(2)}) is off the plot`;
      }
    }

    // ---- readout for step k ----
    const xk = xs[k], fk = F(p, xk);
    const lines: string[] = [];
    lines.push(`<b>step k = ${k}</b>   x${sub(k)} = (${fx(xk[0])}, ${fx(xk[1])})`);
    lines.push(`F(x${sub(k)}) = (${fv(fk[0])}, ${fv(fk[1])})   ‖F‖ = ${fe(norm(fk))}`);
    if (s) {
      const j = s.j;
      lines.push(`J(x${sub(k)}) = [${fx(j[0], 4)}  ${fx(j[1], 4)} ]   det J = ${fe(s.det)}`);
      lines.push(`        [${fx(j[2], 4)}  ${fx(j[3], 4)} ]` + (Math.abs(s.det) < 0.05 ? `   <span class="nw-warn">nearly singular</span>` : ""));
      lines.push(`δ = solve(J, −F) = (${fv(s.d[0])}, ${fv(s.d[1])})   ‖δ‖ = ${fe(norm(s.d))}`);
      if (k < K) lines.push(`x${sub(k + 1)} = x${sub(k)} + δ = (${fx(xs[k + 1][0])}, ${fx(xs[k + 1][1])})`);
    } else {
      lines.push(`‖F‖ < 10⁻¹²: stop. x${sub(k)} is a steady state.`);
    }
    ro.innerHTML = lines.join("\n");

    // ---- where the run ends, and the errors per step ----
    const end = xs[K];
    let res = "";
    let target: Vec | null = null;
    if (path.status === "converged") {
      let best: Root | null = null, bd = 1e-3;
      for (const r of rs) { const dd = Math.hypot(end[0] - r.x[0], end[1] - r.x[1]); if (dd < bd) { bd = dd; best = r; } }
      const name = best ? (rs.length === 1 ? "the only steady state" : KIND_LABEL[best.kind]) : "a root";
      res = `<b>ends:</b> ${name} (${end[0].toFixed(6)}, ${end[1].toFixed(6)}) after <b>${K}</b> step${K === 1 ? "" : "s"}`;
      target = best ? best.x : end;
    } else if (path.status === "nan") {
      res = `<b>ends:</b> <span class="nw-warn">NaN after ${K} step${K === 1 ? "" : "s"}</span> (a negative concentration raised to a non-integer power n)`;
    } else if (path.status === "escaped") {
      res = `<b>ends:</b> <span class="nw-warn">ran off to infinity after ${K} step${K === 1 ? "" : "s"}</span>`;
    } else {
      res = `<b>ends:</b> <span class="nw-warn">no convergence in ${MAXIT} steps</span> (‖F‖ = ${fe(norm(F(p, end)))})`;
    }
    let errs = "";
    if (target) {
      const t = target;
      const es = xs.slice(0, Math.min(K + 1, 12)).map((x, i) => {
        const e = Math.hypot(x[0] - t[0], x[1] - t[1]);
        const txt = e < 1e-16 * 4 ? "≈0" : fe(e);
        return i === k ? `<b>[${txt}]</b>` : txt;
      });
      errs = `\n<b>error ‖xₖ − x*‖:</b> ` + es.join(", ") + (K + 1 > 12 ? ", …" : "");
    }
    fin.innerHTML = res + errs + (shares ? "\n" + shares : "");
  }

  if (st.basins) changedBasins();
  draw();
}

autoMount(NAME, mount);

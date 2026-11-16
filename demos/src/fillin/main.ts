// Fill-in viewer for the structured-sparse unit. A symmetric sparse matrix is factorised
// by Cholesky (symbolically: which entries of L are non-zero), one column at a time.
//  (a) the non-zero pattern, rows and columns in elimination order: blue = non-zero of A,
//      vermillion = fill-in (a zero of A that becomes non-zero in L), outlined = created
//      by the last step; eliminated rows/columns are shaded, the next pivot is boxed;
//  (b) the same matrix as a network (one node per unknown, one link per non-zero pair):
//      eliminating a node links all its remaining neighbours to each other; new links
//      are the fill-in.
// Matrices: tridiagonal (morphogen in a row of cells), band of width 2, the lecture's
// arrow matrix (hub gene linked to all others; the hub is unknown 1 in natural order)
// and the five-point 2D grid Laplacian. Orderings: natural, reversed (hub last),
// shuffled (fixed seed), minimum degree, reverse Cuthill-McKee.
// Operations per column: (c + 1)^2, c = non-zeros below the pivot (1 square root,
// c divisions, c(c + 1)/2 multiply-subtracts); a dense n x n matrix gives about n^3/3.
// The table also gives nnz(L) for the same matrix at n = 1000 (grid 30 x 30), as in the
// lecture (arrow: 500 500 hub first, 1999 hub last).
// Auto-mounts on <div data-demo="fillin">. Optional attributes:
//   data-matrix="tri|band|arrow|grid"  data-n="8" (grid: cells per side)
//   data-order="natural|reversed|shuffled|mindeg|rcm"  data-step="0"
import { svgEl, setAttrs, htmlEl, slider, presets, button, rng, autoMount } from "../lib/ui";

const NAME = "fillin";
// Okabe-Ito, as in make-structured-sparse-figs.R (blue = A, vermillion = fill-in)
const COL = { a: "#0072B2", fill: "#D55E00", grn: "#009E73", done: "#e9e9e9", grey: "#9a9a9a", ink: "#222", zero: "#fff", line: "#d4d4d4" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .fi-row { display: flex; flex-wrap: wrap; gap: .35em 1.1em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .fi-group { display: flex; gap: .4em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .fi-row .qb-slider { min-width: 230px; flex: 1 1 230px; max-width: 340px; grid-template-columns: auto 1fr 2.6em; }
[data-demo="${NAME}"] .fi-main { display: flex; flex-wrap: wrap; gap: 10px; width: 100%; align-items: flex-start; }
[data-demo="${NAME}"] .fi-main svg { flex: 1 1 560px; min-width: 0; max-width: 760px; height: auto; }
[data-demo="${NAME}"] .fi-side { flex: 1 1 300px; max-width: 420px; display: flex; flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] svg text { font-size: 20px; }
[data-demo="${NAME}"] svg .ttl { font-size: 21px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 19px; fill: #555; }
[data-demo="${NAME}"] .fi-count { display: grid; grid-template-columns: auto auto; gap: .05em .8em; }
[data-demo="${NAME}"] .fi-count b { font-weight: 600; }
[data-demo="${NAME}"] .fi-count .v { text-align: right; }
[data-demo="${NAME}"] .qb-readout.fi-info { font-family: inherit; }
[data-demo="${NAME}"] .fi-fill { color: ${COL.fill}; font-weight: 700; }
[data-demo="${NAME}"] table.fi-tab { border-collapse: collapse; font-family: ui-monospace, monospace; font-size: .85em; width: 100%; }
[data-demo="${NAME}"] table.fi-tab td, [data-demo="${NAME}"] table.fi-tab th { padding: .08em .35em; text-align: right; white-space: nowrap; }
[data-demo="${NAME}"] table.fi-tab th { font-family: system-ui, sans-serif; font-weight: 600; border-bottom: 1px solid #bbb; }
[data-demo="${NAME}"] table.fi-tab td:first-child, [data-demo="${NAME}"] table.fi-tab th:first-child { text-align: left; }
[data-demo="${NAME}"] table.fi-tab tr.cur td { font-weight: 700; background: #fde7d8; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; gap: .25em; }
.reveal [data-demo="${NAME}"] .qb-readout { font-size: 17px; padding: .25em .45em; }
.reveal [data-demo="${NAME}"] .qb-btn { font-size: .85em; padding: .15em .55em; }
.reveal [data-demo="${NAME}"] .fi-main svg { flex-basis: 640px; max-width: 750px; }
.reveal [data-demo="${NAME}"] .fi-side { flex-basis: 360px; }
.reveal [data-demo="${NAME}"] .qb-note { display: none; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---------------------------------------------------------------- matrices
type MType = "tri" | "band" | "arrow" | "grid";
type Order = "natural" | "reversed" | "shuffled" | "mindeg" | "rcm";
const ORDERS: { key: Order; label: string }[] = [
  { key: "natural", label: "natural" },
  { key: "reversed", label: "reversed" },
  { key: "shuffled", label: "shuffled" },
  { key: "mindeg", label: "min degree" },
  { key: "rcm", label: "RCM" },
];
const TYPES: { key: MType; label: string }[] = [
  { key: "tri", label: "tridiagonal" },
  { key: "band", label: "band (width 2)" },
  { key: "arrow", label: "arrow (hub gene)" },
  { key: "grid", label: "2D grid" },
];
const SIZE: Record<MType, { min: number; max: number; def: number; big: number }> = {
  tri: { min: 4, max: 30, def: 12, big: 1000 },
  band: { min: 5, max: 30, def: 12, big: 1000 },
  arrow: { min: 4, max: 30, def: 8, big: 1000 },
  grid: { min: 3, max: 7, def: 5, big: 30 },          // cells per side
};

interface Mat { n: number; adj: Set<number>[]; edges: number }

/** Off-diagonal pattern of the matrix (symmetric), as adjacency sets on the unknowns. */
function build(t: MType, size: number): Mat {
  const n = t === "grid" ? size * size : size;
  const adj = Array.from({ length: n }, () => new Set<number>());
  const link = (i: number, j: number) => { adj[i].add(j); adj[j].add(i); };
  if (t === "tri") for (let i = 0; i + 1 < n; i++) link(i, i + 1);
  if (t === "band") for (let i = 0; i < n; i++) for (let d = 1; d <= 2; d++) if (i + d < n) link(i, i + d);
  if (t === "arrow") for (let i = 1; i < n; i++) link(0, i);
  if (t === "grid") for (let r = 0; r < size; r++) for (let c = 0; c < size; c++) {
    const i = r * size + c;
    if (c + 1 < size) link(i, i + 1);
    if (r + 1 < size) link(i, i + size);
  }
  let edges = 0;
  adj.forEach(s => { edges += s.size; });
  return { n, adj, edges: edges / 2 };
}

const copyAdj = (m: Mat) => m.adj.map(s => new Set(s));

interface Elim {
  perm: number[];              // perm[k] = unknown eliminated at step k
  pos: number[];               // pos[u] = its row/column in the reordered matrix
  c: number[];                 // non-zeros below the pivot in column k of L
  fill: Map<number, number>;   // key min*n+max -> step that created it
  nnzL: number; ops: number;
}

/** Symbolic Cholesky: eliminate in the given order (or by minimum degree if order is null). */
function eliminate(m: Mat, order: number[] | null): Elim {
  const n = m.n, adj = copyAdj(m), done = new Uint8Array(n);
  const perm: number[] = [], c: number[] = [], fill = new Map<number, number>();
  for (let k = 0; k < n; k++) {
    let p: number;
    if (order) p = order[k];
    else {
      p = -1;
      for (let u = 0; u < n; u++) if (!done[u] && (p < 0 || adj[u].size < adj[p].size)) p = u;
    }
    const nb = [...adj[p]];
    for (let a = 0; a < nb.length; a++) for (let b = a + 1; b < nb.length; b++) {
      const u = nb[a], v = nb[b];
      if (!adj[u].has(v)) {
        adj[u].add(v); adj[v].add(u);
        fill.set(Math.min(u, v) * n + Math.max(u, v), k);
      }
    }
    for (const q of nb) adj[q].delete(p);
    adj[p].clear();
    done[p] = 1;
    perm.push(p); c.push(nb.length);
  }
  const pos = new Array<number>(n);
  perm.forEach((u, k) => { pos[u] = k; });
  const nnzL = n + c.reduce((s, v) => s + v, 0);
  const ops = c.reduce((s, v) => s + (v + 1) * (v + 1), 0);
  return { perm, pos, c, fill, nnzL, ops };
}

function rcm(m: Mat): number[] {
  const n = m.n, seen = new Uint8Array(n), out: number[] = [];
  const deg = (u: number) => m.adj[u].size;
  while (out.length < n) {
    let s = -1;
    for (let u = 0; u < n; u++) if (!seen[u] && (s < 0 || deg(u) < deg(s))) s = u;
    const q = [s]; seen[s] = 1;
    while (q.length) {
      const u = q.shift() as number;
      out.push(u);
      const nb = [...m.adj[u]].filter(v => !seen[v]).sort((a, b) => deg(a) - deg(b) || a - b);
      for (const v of nb) { seen[v] = 1; q.push(v); }
    }
  }
  return out.reverse();
}

function elimFor(m: Mat, o: Order): Elim {
  const n = m.n, id = Array.from({ length: n }, (_, i) => i);
  if (o === "natural") return eliminate(m, id);
  if (o === "reversed") return eliminate(m, id.reverse());
  if (o === "shuffled") {
    const r = rng(7 + n);
    for (let i = n - 1; i > 0; i--) { const j = Math.floor(r.unif() * (i + 1)); [id[i], id[j]] = [id[j], id[i]]; }
    return eliminate(m, id);
  }
  if (o === "rcm") return eliminate(m, rcm(m));
  return eliminate(m, null);
}

const bigCache = new Map<string, number>();
function bigNnz(t: MType, o: Order): number {
  const key = t + o;
  if (!bigCache.has(key)) bigCache.set(key, elimFor(build(t, SIZE[t].big), o).nnzL);
  return bigCache.get(key) as number;
}

const fmtInt = (v: number): string => Math.round(v).toString().replace(/\B(?=(\d{3})+(?!\d))/g, " ");

const BIO: Record<MType, (size: number) => string> = {
  tri: n => `Morphogen in a row of ${n} cells: each cell exchanges with its two neighbours (Thomas algorithm).`,
  band: n => `${n} cells exchanging with neighbours up to two cells away: bandwidth 2.`,
  arrow: n => `Gene network of ${n} genes: one hub (unknown 1 in natural order) linked to all others.`,
  grid: s => `A sheet of ${s} × ${s} cells, diffusion to the 4 neighbours (2D five-point Laplacian), numbered row by row.`,
};

// ---------------------------------------------------------------- the demo
function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const t0 = (TYPES.some(x => x.key === ds.matrix) ? ds.matrix : "arrow") as MType;
  const n0raw = Number(ds.n);
  const clampSize = (t: MType, v: number) => Math.max(SIZE[t].min, Math.min(SIZE[t].max, Math.round(v)));
  const size0 = Number.isFinite(n0raw) && ds.n !== undefined ? clampSize(t0, n0raw) : SIZE[t0].def;
  const o0 = (ORDERS.some(x => x.key === ds.order) ? ds.order : "natural") as Order;
  const s0raw = Number(ds.step);
  const st = { t: t0, size: size0, o: o0, s: 0 };
  let mat = build(st.t, st.size);
  let el = elimFor(mat, st.o);
  st.s = Number.isFinite(s0raw) ? Math.max(0, Math.min(mat.n, Math.round(s0raw))) : 0;
  let timer: number | null = null;

  // ---- controls ----
  const r1 = htmlEl("div", "fi-row", "", host);
  const g1 = htmlEl("div", "fi-group", "", r1);
  htmlEl("span", "", "matrix:", g1);
  const tBtns = presets(g1, TYPES.map(x => ({ label: x.label, apply: () => setType(x.key) })));
  const sz = slider(r1, "size n", SIZE[st.t].min, SIZE[st.t].max, 1, st.size, v => { st.size = v; rebuild(); });
  const r2 = htmlEl("div", "fi-row", "", host);
  const g2 = htmlEl("div", "fi-group", "", r2);
  htmlEl("span", "", "order:", g2);
  const oBtns = presets(g2, ORDERS.map(x => ({ label: x.label, apply: () => { st.o = x.key; recompute(); } })));
  const g3 = htmlEl("div", "fi-group", "", r2);
  button(g3, "⏮ restart", () => { stop(); st.s = 0; draw(); });
  button(g3, "◀ back", () => { stop(); st.s = Math.max(0, st.s - 1); draw(); });
  button(g3, "step ▶", () => { stop(); st.s = Math.min(mat.n, st.s + 1); draw(); });
  const playB = button(g3, "play", () => { if (timer !== null) stop(); else play(); });
  button(g3, "all ⏭", () => { stop(); st.s = mat.n; draw(); });
  button(g3, "reset", () => {
    stop(); st.t = t0; st.size = size0; st.o = o0; setSliderRange();
    mat = build(st.t, st.size); el = elimFor(mat, st.o);
    st.s = Number.isFinite(s0raw) ? Math.max(0, Math.min(mat.n, Math.round(s0raw))) : 0; draw();
  });

  const main = htmlEl("div", "fi-main", "", host);
  const W = 900, H = 540;
  const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Non-zero pattern and network during Cholesky elimination" }, main);
  const side = htmlEl("div", "fi-side", "", main);
  const count = htmlEl("div", "qb-readout fi-count", "", side);
  const stepInfo = htmlEl("div", "qb-readout fi-info", "", side);
  const tabBox = htmlEl("div", "", "", side);
  const note = htmlEl("div", "qb-note", "", host);

  function stop(): void { if (timer !== null) { clearInterval(timer); timer = null; } playB.textContent = "play"; }
  function play(): void {
    if (st.s >= mat.n) st.s = 0;
    playB.textContent = "pause";
    const delay = Math.max(150, Math.min(700, 9000 / mat.n));
    timer = window.setInterval(() => {
      if (!host.isConnected) { stop(); return; }
      st.s++; draw();
      if (st.s >= mat.n) stop();
    }, delay);
    draw();
  }
  function setSliderRange(): void {
    const r = SIZE[st.t];
    sz.input.min = String(r.min); sz.input.max = String(r.max);
    sz.set(st.size);
    const name = sz.input.parentElement?.querySelector(".qb-slider-name");
    if (name) name.textContent = st.t === "grid" ? "grid side" : "size n";
  }
  function setType(t: MType): void {
    stop();
    st.t = t; st.size = SIZE[t].def;
    setSliderRange(); rebuild();
  }
  function rebuild(): void { stop(); mat = build(st.t, st.size); recompute(); }
  function recompute(): void { stop(); el = elimFor(mat, st.o); st.s = 0; draw(); }

  // ---- layout of the network ----
  const GX0 = 560, GX1 = 890, GY0 = 50, GY1 = 380;
  function nodePos(u: number): [number, number] {
    const cx = (GX0 + GX1) / 2, cy = (GY0 + GY1) / 2, R = (GX1 - GX0) / 2 - 18;
    if (st.t === "grid") {
      const sd = st.size, step = (GX1 - GX0 - 40) / (sd - 1);
      return [GX0 + 20 + (u % sd) * step, GY0 + 20 + Math.floor(u / sd) * step];
    }
    if (st.t === "arrow") {
      if (u === 0) return [cx, cy];
      const a = -Math.PI / 2 + 2 * Math.PI * (u - 1) / (mat.n - 1);
      return [cx + R * Math.cos(a), cy + R * Math.sin(a)];
    }
    // a row of cells drawn around a ring (gap between the last and the first cell)
    const a = -Math.PI / 2 + 2 * Math.PI * (u + 0.5) / (mat.n + 1);
    return [cx + R * Math.cos(a), cy + R * Math.sin(a)];
  }

  function draw(): void {
    const n = mat.n, s = st.s;
    tBtns.forEach((b, i) => b.classList.toggle("active", TYPES[i].key === st.t));
    oBtns.forEach((b, i) => b.classList.toggle("active", ORDERS[i].key === st.o));
    const isA = (u: number, v: number) => mat.adj[u].has(v);
    const fillStep = (u: number, v: number) => el.fill.get(Math.min(u, v) * n + Math.max(u, v));

    svg.replaceChildren();
    // ---------- (a) matrix ----------
    const MX = 40, MY = 52, MS = 470, cs = MS / n;
    svgEl("text", { x: 8, y: 26, class: "ttl" }, svg).textContent = "(a) pattern, in elimination order";
    svgEl("rect", { x: MX, y: MY, width: MS, height: MS, fill: COL.zero, stroke: COL.grey, "stroke-width": 1.5 }, svg);
    if (s > 0) {   // eliminated rows and columns
      svgEl("rect", { x: MX, y: MY, width: s * cs, height: MS, fill: COL.done }, svg);
      svgEl("rect", { x: MX, y: MY, width: MS, height: s * cs, fill: COL.done }, svg);
    }
    if (n <= 30) for (let i = 1; i < n; i++) {
      svgEl("line", { x1: MX + i * cs, x2: MX + i * cs, y1: MY, y2: MY + MS, stroke: COL.line, "stroke-width": 0.8 }, svg);
      svgEl("line", { y1: MY + i * cs, y2: MY + i * cs, x1: MX, x2: MX + MS, stroke: COL.line, "stroke-width": 0.8 }, svg);
    }
    const pad = Math.min(1.5, cs * 0.08);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
      const u = el.perm[i], v = el.perm[j];
      let col = "", fresh = false;
      if (i === j || isA(u, v)) col = COL.a;
      else {
        const f = fillStep(u, v);
        if (f !== undefined && f < s) { col = COL.fill; fresh = f === s - 1; }
      }
      if (!col) continue;
      const r = svgEl("rect", { x: MX + j * cs + pad, y: MY + i * cs + pad, width: cs - 2 * pad, height: cs - 2 * pad, fill: col }, svg);
      if (Math.min(i, j) < s) r.setAttribute("fill-opacity", "0.45");   // finished row/column
      if (fresh) setStroke(r, COL.ink, Math.max(1.5, Math.min(3, cs * 0.12)));
    }
    if (s < n) {   // next pivot and its column
      svgEl("rect", { x: MX + s * cs, y: MY + s * cs, width: cs, height: MS - s * cs, fill: "none", stroke: COL.ink, "stroke-width": 1.5, "stroke-dasharray": "5 4" }, svg);
      svgEl("rect", { x: MX + s * cs, y: MY + s * cs, width: cs, height: cs, fill: "none", stroke: COL.ink, "stroke-width": 3 }, svg);
    }
    // index labels
    const lab = (txt: string, x: number, y: number, anchor: string) => {
      svgEl("text", { x, y, "text-anchor": anchor, class: "sm" }, svg).textContent = txt;
    };
    lab("1", MX + cs / 2, MY - 6, "middle");
    lab(String(n), MX + MS - cs / 2, MY - 6, "middle");
    lab("1", MX - 6, MY + cs / 2 + 6, "end");
    lab(String(n), MX - 6, MY + MS - cs / 2 + 6, "end");

    // ---------- (b) network ----------
    svgEl("text", { x: GX0 - 8, y: 26, class: "ttl" }, svg).textContent = "(b) the network";
    const done = (u: number) => el.pos[u] < s;
    const next = s < n ? el.perm[s] : -1;
    const P = Array.from({ length: n }, (_, u) => nodePos(u));
    // links of A are straight; fill-in links bend a little, so they do not hide behind
    // the links of A that run along the same row of nodes
    const edge = (u: number, v: number, stroke: string, w: number, op: number, dash = "", bend = 0) => {
      const [x1, y1] = P[u], [x2, y2] = P[v];
      const mx = (x1 + x2) / 2 - bend * (y2 - y1), my = (y1 + y2) / 2 + bend * (x2 - x1);
      svgEl("path", { d: `M${x1},${y1} Q${mx},${my} ${x2},${y2}`, fill: "none", stroke, "stroke-width": w, "stroke-opacity": op, "stroke-dasharray": dash }, svg);
    };
    for (let u = 0; u < n; u++) for (const v of mat.adj[u]) if (u < v) {
      const gone = done(u) || done(v);
      edge(u, v, gone ? COL.grey : COL.a, gone ? 1.5 : 3, gone ? 0.5 : 1, gone ? "4 4" : "");
    }
    for (const [key, f] of el.fill) {
      if (f >= s) continue;
      const u = Math.floor(key / n), v = key % n;
      const gone = done(u) || done(v);
      edge(u, v, COL.fill, gone ? 1.5 : (f === s - 1 ? 4.5 : 3), gone ? 0.45 : 1, gone ? "4 4" : "", st.t === "arrow" ? 0 : 0.18);
    }
    const nr = Math.max(6, Math.min(16, 150 / Math.sqrt(n) / 1.6));
    const showNum = nr >= 13;
    for (let u = 0; u < n; u++) {
      const [x, y] = P[u];
      const c = svgEl("circle", { cx: x, cy: y, r: u === next ? nr + 3 : nr }, svg);
      if (done(u)) setAttrs(c, { fill: "#fff", stroke: COL.grey, "stroke-width": 2 });
      else if (u === next) setAttrs(c, { fill: COL.grn, stroke: COL.ink, "stroke-width": 3 });
      else setAttrs(c, { fill: COL.a, stroke: "#fff", "stroke-width": 1.5 });
      if (showNum) {
        const t = svgEl("text", { x, y: y + 6, "text-anchor": "middle", "font-size": 16, "font-weight": 600 }, svg);
        t.style.fill = done(u) ? COL.grey : "#fff";
        t.textContent = String(el.pos[u] + 1);
      }
    }
    // legend
    const ly = 412;
    const key = (y: number, txt: string, f: (x: number, y: number) => void) => {
      f(GX0 + 4, y - 7);
      svgEl("text", { x: GX0 + 34, y }, svg).textContent = txt;
    };
    const sq = (fill: string, stroke = "") => (x: number, y: number) => {
      const r = svgEl("rect", { x, y: y - 9, width: 18, height: 18, fill }, svg);
      if (stroke) setStroke(r, stroke, 2.5);
    };
    key(ly, "non-zero of A", sq(COL.a));
    key(ly + 28, "fill-in (zero in A)", sq(COL.fill));
    key(ly + 56, "fill-in from the last step", sq(COL.fill, COL.ink));
    key(ly + 84, "next pivot  (number = order)", (x, y) => {
      svgEl("circle", { cx: x + 9, cy: y, r: 9, fill: COL.grn, stroke: COL.ink, "stroke-width": 2.5 }, svg);
    });
    key(ly + 112, "eliminated", (x, y) => { svgEl("rect", { x, y: y - 9, width: 18, height: 18, fill: COL.done, stroke: COL.grey }, svg); });

    // ---------- counters ----------
    const nnzA = n + 2 * mat.edges, lowA = n + mat.edges;
    let fillSoFar = 0;
    for (const f of el.fill.values()) if (f < s) fillSoFar++;
    let lSoFar = 0, opsSoFar = 0;
    for (let k = 0; k < s; k++) { lSoFar += 1 + el.c[k]; opsSoFar += (el.c[k] + 1) ** 2; }
    count.replaceChildren();
    const row = (k: string, v: string, cls = "") => {
      htmlEl("b", "", k, count);
      htmlEl("span", "v " + cls, v, count);
    };
    row("columns eliminated", `${s} / ${n}`);
    row("non-zeros of A", `${fmtInt(nnzA)}  (lower: ${fmtInt(lowA)})`);
    row("non-zeros of L so far", fmtInt(lSoFar));
    row("fill-in so far", fmtInt(fillSoFar), "fi-fill");
    row("operations so far", fmtInt(opsSoFar));
    row("at the end: nnz(L)", `${fmtInt(el.nnzL)}  (fill ${fmtInt(el.nnzL - lowA)})`);

    if (s === 0) stepInfo.textContent = "Click step ▶: eliminate column 1. Its non-zeros below the pivot are the neighbours of the green node.";
    else {
      const k = s - 1, c = el.c[k];
      let created = 0;
      for (const f of el.fill.values()) if (f === k) created++;
      stepInfo.textContent = `Step ${s}: the pivot had c = ${c} remaining neighbour${c === 1 ? "" : "s"}; ` +
        `elimination links them all: ${c * (c - 1) / 2} pair${c * (c - 1) / 2 === 1 ? "" : "s"}, ${created} new (fill-in). ` +
        `Operations (c + 1)² = ${(c + 1) ** 2}.`;
    }

    // ---------- table of orderings ----------
    const big = SIZE[st.t].big, bigLabel = st.t === "grid" ? `${big}×${big}` : `n = ${fmtInt(big)}`;
    tabBox.replaceChildren();
    const tab = htmlEl("table", "fi-tab", "", tabBox);
    const hr = htmlEl("tr", "", "", tab);
    for (const h of ["order", "nnz(L)", "operations", bigLabel + ": nnz(L)"]) htmlEl("th", "", h, hr);
    for (const o of ORDERS) {
      const e = o.key === st.o ? el : elimFor(mat, o.key);
      const tr = htmlEl("tr", o.key === st.o ? "cur" : "", "", tab);
      htmlEl("td", "", o.label, tr);
      htmlEl("td", "", fmtInt(e.nnzL), tr);
      htmlEl("td", "", fmtInt(e.ops), tr);
      htmlEl("td", "", fmtInt(bigNnz(st.t, o.key)), tr);
    }
    const dense = n * (n + 1) / 2;
    htmlEl("div", "qb-note", `Full lower triangle: ${fmtInt(dense)}. ${st.t === "grid" ? `${bigLabel} grid` : bigLabel}, A lower: ${fmtInt(st.t === "grid" ? big * big + 2 * big * (big - 1) : st.t === "band" ? 3 * big - 3 : 2 * big - 1)}.`, tabBox);

    note.textContent = BIO[st.t](st.size) + " Operations per column (c + 1)², c = non-zeros below the pivot; a full matrix needs about n³/3.";
    playB.classList.toggle("active", timer !== null);
  }

  function setStroke(e: Element, c: string, w: number): void { e.setAttribute("stroke", c); e.setAttribute("stroke-width", String(w)); }

  tBtns.forEach((b, i) => b.classList.toggle("active", TYPES[i].key === st.t));
  setSliderRange();
  draw();
}

autoMount(NAME, mount);

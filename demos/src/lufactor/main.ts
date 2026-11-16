// Gaussian elimination as LU, for the lu unit. Two tabs:
//  1. Factorise: elimination on the mixing matrix A one multiplier at a time. Each
//     multiplier l_ik = u_ik / u_kk is written into L at the spot where it made its zero
//     in U; the product L·U is shown next to them and equals A at every step. Exact
//     fractions (lib/frac.ts). A zero pivot stops elimination (with what lu_nopivot would
//     store: NaN, -Inf), and a button swaps the rows and restarts.
//  2. Solve: L and U are kept; each sample (consumer) b costs only forward substitution
//     L y = b and back substitution U x = y, written out line by line. The sources and the
//     samples are drawn in isotope space (n = 3 with a "sum to 1" first row), and an
//     operation counter compares "factor once, then substitute" (F + m S) with
//     "re-eliminate for every sample" (m (F + S)), for this n and for n = 1000.
// Operation counts as in the lecture: a division, a multiplication and a subtraction count
// one each. Factorising: sum_k (n-k)(1 + 2(n-k)) (13 for n = 3); one forward + back
// substitution with a unit-diagonal L: 2n^2 - n (15 for n = 3).
// Auto-mounts on <div data-demo="lufactor">. Optional attributes:
//   data-matrix="mixing|zeropivot"  (the lecture's two matrices; default mixing)
//   data-a="1 1 1 -26 -10 -18 2 6 14"  (n×n, row by row, n = 2..4; overrides data-matrix)
//   data-b="1 -19.6 5.6; 1 -14.8 6.8"  (consumers, ';'-separated)
//   data-names="C3 plants,maize,fish"  data-tab="factor|solve"  data-n="3|1000"  data-seed="3"
import { svgEl, htmlEl, button, presets, scale, rng, autoMount } from "../lib/ui";
import { Frac, ZERO, parseFrac } from "../lib/frac";
import { M, cloneM, identity, matMul } from "../lib/rowops";
import { fr } from "../lib/rowopsui";

const NAME = "lufactor";
// Okabe-Ito, as in make-lu-figs.R
const COL = { blue: "#0072B2", sky: "#56B4E9", orange: "#E69F00", red: "#D55E00", green: "#009E73", grey: "#999999" };
const STYLE = `
[data-demo="${NAME}"].qb-demo { flex-direction: column; gap: .4em; }
[data-demo="${NAME}"] .lu-row { display: flex; flex-wrap: wrap; gap: .4em 1.2em; align-items: center; width: 100%; }
[data-demo="${NAME}"] .lu-group { display: flex; gap: .4em; align-items: center; flex-wrap: wrap; }
[data-demo="${NAME}"] .lu-tabs .qb-btn { font-size: .95em; }
[data-demo="${NAME}"] .lu-count { font-size: .85em; color: var(--qb-muted); }
[data-demo="${NAME}"] .lu-mats { display: flex; align-items: center; gap: .9em; flex-wrap: wrap; }
[data-demo="${NAME}"] .lu-op { font-size: 1.6em; color: var(--qb-muted); }
[data-demo="${NAME}"] .lu-cap { font-size: .8em; color: var(--qb-muted); text-align: center; }
[data-demo="${NAME}"] .lu-cap b { color: var(--qb-fg); font-size: 1.15em; }
[data-demo="${NAME}"] table { margin: 0; width: auto; border: 0; }
[data-demo="${NAME}"] table td { border: 0; line-height: 1.2; vertical-align: middle; }
[data-demo="${NAME}"] table.lu-m { border-collapse: separate; border-spacing: 0; font-family: ui-monospace, monospace; }
[data-demo="${NAME}"] table.lu-m td { padding: .1em .4em; text-align: center; min-width: 1.7em; height: 1.75em; }
[data-demo="${NAME}"] table.lu-m td.rl { font-family: system-ui, sans-serif; font-size: .75em; color: var(--qb-muted); padding-right: .25em; min-width: 0; }
[data-demo="${NAME}"] table.lu-m td.bk { min-width: 0; width: .3em; padding: 0; }
[data-demo="${NAME}"] table.lu-m td.bl { border-left: 2px solid var(--qb-fg); }
[data-demo="${NAME}"] table.lu-m td.br { border-right: 2px solid var(--qb-fg); }
[data-demo="${NAME}"] table.lu-m tr.top td.bk { border-top: 2px solid var(--qb-fg); }
[data-demo="${NAME}"] table.lu-m tr.bot td.bk { border-bottom: 2px solid var(--qb-fg); }
[data-demo="${NAME}"] tr.chg td:not(.rl):not(.bk) { background: #fff1b8; }
[data-demo="${NAME}"] tr.src td:not(.rl):not(.bk) { background: #dcecfa; }
[data-demo="${NAME}"] td.faint { color: #b8b8b8; }
[data-demo="${NAME}"] td.pv { box-shadow: inset 0 0 0 2.5px ${COL.blue}; border-radius: 50%; font-weight: 700; }
[data-demo="${NAME}"] td.nw { box-shadow: inset 0 0 0 2.5px ${COL.green}; color: #00704f; font-weight: 700; border-radius: 4px; }
[data-demo="${NAME}"] td.zp { box-shadow: inset 0 0 0 2.5px ${COL.red}; color: ${COL.red}; font-weight: 700; border-radius: 50%; }
[data-demo="${NAME}"] td.bad { color: ${COL.red}; font-weight: 700; }
[data-demo="${NAME}"] .ro-fr { display: inline-flex; flex-direction: column; vertical-align: middle; text-align: center; font-size: .78em; line-height: 1.05; }
[data-demo="${NAME}"] .ro-fr > span:first-child { border-bottom: 1.5px solid currentColor; padding: 0 .1em; }
[data-demo="${NAME}"] .lu-title { font-size: 1.15em; font-weight: 700; }
[data-demo="${NAME}"] .lu-text p { margin: .15em 0; line-height: 1.45; }
[data-demo="${NAME}"] .lu-eq { line-height: 1.7; font-variant-numeric: tabular-nums; }
[data-demo="${NAME}"] .lu-eq .h { font-family: system-ui, sans-serif; font-weight: 600; }
[data-demo="${NAME}"] .lu-ok { color: #00704f; font-weight: 700; }
[data-demo="${NAME}"] .lu-no { color: ${COL.red}; font-weight: 700; }
[data-demo="${NAME}"] .lu-key { font-size: .75em; color: var(--qb-muted); display: flex; gap: 1em; flex-wrap: wrap; }
[data-demo="${NAME}"] .lu-key span::before { content: ""; display: inline-block; width: .9em; height: .9em; margin-right: .3em; vertical-align: -.1em; border-radius: 3px; }
[data-demo="${NAME}"] .lu-key .k-pv::before { box-shadow: inset 0 0 0 2px ${COL.blue}; border-radius: 50%; }
[data-demo="${NAME}"] .lu-key .k-src::before { background: #dcecfa; }
[data-demo="${NAME}"] .lu-key .k-chg::before { background: #fff1b8; }
[data-demo="${NAME}"] .lu-key .k-nw::before { box-shadow: inset 0 0 0 2px ${COL.green}; }
[data-demo="${NAME}"] .lu-solve { display: flex; flex-wrap: wrap; gap: 14px; width: 100%; align-items: flex-start; }
[data-demo="${NAME}"] .lu-left { flex: 1 1 520px; min-width: 0; display: flex; flex-direction: column; gap: .35em; }
[data-demo="${NAME}"] .lu-right { flex: 1 1 460px; max-width: 560px; min-width: 0; display: flex; flex-direction: column; gap: .3em; }
[data-demo="${NAME}"] .lu-right svg { width: 100%; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .ttl { font-size: 19px; font-weight: 600; }
[data-demo="${NAME}"] svg .sm { font-size: 18px; fill: #555; }
.reveal [data-demo="${NAME}"].qb-demo { font-size: 19px; gap: .3em; }
.reveal [data-demo="${NAME}"] table.lu-m td { height: 1.55em; padding: .04em .35em; }
.reveal [data-demo="${NAME}"] .lu-eq { line-height: 1.55; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

// ---------------------------------------------------------------- the lecture's examples
interface Example { a: number[][]; names: string[]; b: number[][] }
const EXAMPLES: Record<string, Example> = {
  mixing: {
    a: [[1, 1, 1], [-26, -10, -18], [2, 6, 14]],
    names: ["C3 plants", "maize", "fish"],
    b: [[1, -19.6, 5.6], [1, -14.8, 6.8], [1, -17.2, 11.2], [1, -24, 9]],
  },
  zeropivot: {
    a: [[1, 1, 1], [-26, -26, -10], [2, 10, 6]],
    names: ["C3 plants", "fw fish", "maize"],
    b: [[1, -20, 6.5]],
  },
};

// ---------------------------------------------------------------- operation counts
const factorOps = (n: number): number => { let s = 0; for (let j = 1; j < n; j++) s += j * (1 + 2 * j); return s; };
const substOps = (n: number): number => 2 * n * n - n;

// ---------------------------------------------------------------- formatting
/** Exact decimal when the fraction terminates within 6 digits (6.4, 0.375), else a stacked fraction. */
function dec(f: Frac): string {
  let d = f.d, k2 = 0, k5 = 0;
  while (d % 2n === 0n) { d /= 2n; k2++; }
  while (d % 5n === 0n) { d /= 5n; k5++; }
  const k = Math.max(k2, k5);
  if (d !== 1n || k > 6) return fr(f);
  if (k === 0) return fr(f);
  const neg = f.n < 0n, a = neg ? -f.n : f.n;
  const s = (a * 10n ** BigInt(k) / f.d).toString().padStart(k + 1, "0");
  return (neg ? "−" : "") + s.slice(0, -k) + "." + s.slice(-k);
}
const par = (s: string, f: Frac) => (f.sign() < 0 ? `(${s})` : s);
const sub = (s: string) => s.replace(/[0-9]/g, c => "₀₁₂₃₄₅₆₇₈₉"[Number(c)]);
const vecTxt = (v: Frac[], f = dec) => `(${v.map(x => f(x)).join(", ")})`;
function sci(v: number): string {
  if (v < 1e4) return String(Math.round(v));
  const [m, e] = v.toExponential(1).split("e");
  return `${m}·10${e.replace("+", "").replace(/[-0-9]/g, c => (c === "-" ? "⁻" : "⁰¹²³⁴⁵⁶⁷⁸⁹"[Number(c)]))}`;
}

type CellCls = (i: number, j: number) => string;
function matrix(m: (string)[][], cap: string, cell: CellCls = () => "", rowCls: (i: number) => string = () => "", labels = false): string {
  const n = m.length;
  let h = `<div><div class="lu-cap">${cap}</div><table class="lu-m">`;
  for (let i = 0; i < n; i++) {
    h += `<tr class="${i === 0 ? "top" : ""} ${i === n - 1 ? "bot" : ""} ${rowCls(i)}">`;
    if (labels) h += `<td class="rl">R${i + 1}</td>`;
    h += `<td class="bk bl"></td>`;
    for (let j = 0; j < m[i].length; j++) h += `<td class="${cell(i, j)}">${m[i][j]}</td>`;
    h += `<td class="bk br"></td></tr>`;
  }
  return h + `</table></div>`;
}

// ---------------------------------------------------------------- elimination steps
type Step =
  | { kind: "start"; L: M; U: M }
  | { kind: "elim"; k: number; i: number; l: Frac; L: M; U: M; before: M; ops: number }
  | { kind: "zero"; k: number; i: number; L: M; U: M; ops: number; singular: boolean }
  | { kind: "done"; L: M; U: M; ops: number };

function eliminate(A: M): Step[] {
  const n = A.length;
  let L = identity(n), U = cloneM(A), ops = 0;
  const steps: Step[] = [{ kind: "start", L, U }];
  for (let k = 0; k < n - 1; k++) {
    for (let i = k + 1; i < n; i++) {
      if (U[k][k].isZero()) {
        const singular = U.slice(k).every(r => r[k].isZero());
        steps.push({ kind: "zero", k, i, L, U, ops, singular });
        return steps;
      }
      const l = U[i][k].div(U[k][k]);
      const before = U;
      U = cloneM(U); L = cloneM(L);
      L[i][k] = l;
      for (let j = k; j < n; j++) U[i][j] = j === k ? ZERO : before[i][j].sub(l.mul(before[k][j]));
      ops += 1 + 2 * (n - k - 1);
      steps.push({ kind: "elim", k, i, l, L, U, before, ops });
    }
  }
  steps.push({ kind: "done", L, U, ops });
  return steps;
}

// ---------------------------------------------------------------- the widget
function mount(host: HTMLElement): void {
  injectStyle();
  const ds = host.dataset;
  const seed = Number.isFinite(Number(ds.seed)) && ds.seed ? Number(ds.seed) : 3;

  // ---- custom matrix from data-a, else a named example ----
  const parseList = (s: string) => s.trim().split(/[\s,]+/).filter(t => t).map(parseFrac);
  let custom: { A: M; names: string[]; B: Frac[][] } | null = null;
  if (ds.a) {
    const v = parseList(ds.a), n = Math.round(Math.sqrt(v.length));
    if (n >= 2 && n <= 4 && n * n === v.length && v.every(x => x)) {
      const A = Array.from({ length: n }, (_, i) => v.slice(i * n, i * n + n) as Frac[]);
      const B = (ds.b ?? "").split(";").map(parseList).filter(r => r.length === n && r.every(x => x)) as Frac[][];
      const names = (ds.names ?? "").split(",").map(s => s.trim()).filter(s => s);
      custom = { A, names: names.length === n ? names : A.map((_, j) => `source ${j + 1}`), B };
    }
  }
  const toFrac = (rows: number[][]): Frac[][] => rows.map(r => r.map(v => parseFrac(String(v))!));
  const load = (key: string) => {
    if (key === "custom" && custom) return { A: custom.A, names: custom.names, B: custom.B };
    const ex = EXAMPLES[key] ?? EXAMPLES.mixing;
    const B = ds.b && !custom && key === (ds.matrix ?? "mixing")
      ? (ds.b.split(";").map(parseList).filter(r => r.length === ex.a.length && r.every(x => x)) as Frac[][])
      : toFrac(ex.b);
    const names = ds.names && !custom && key === (ds.matrix ?? "mixing") ? ds.names.split(",").map(s => s.trim()) : ex.names;
    return { A: toFrac(ex.a), names: names.length === ex.a.length ? names : ex.names, B: B.length ? B : toFrac(ex.b) };
  };
  const startKey = custom ? "custom" : (ds.matrix === "zeropivot" ? "zeropivot" : "mixing");

  const st = {
    key: startKey,
    A0: [] as M,             // the matrix as given (columns = sources)
    names: [] as string[],
    B: [] as Frac[][],       // preset consumers
    perm: [] as number[],    // row order used for elimination
    steps: [] as Step[],
    cur: 0,
    tab: (ds.tab === "solve" ? "solve" : "factor") as "factor" | "solve",
    nShow: ds.n === "1000" ? 1000 : 3,
    solved: [] as { b: Frac[]; x: Frac[]; label: string }[],
    current: -1,
    r: rng(seed),
  };

  // ---- top: tabs and matrix choice ----
  const top = htmlEl("div", "lu-row", "", host);
  const tabs = htmlEl("div", "lu-group lu-tabs", "", top);
  const tabBtns = presets(tabs, [
    { label: "1. Factorise A = LU", apply: () => showTab("factor") },
    { label: "2. Solve for samples", apply: () => showTab("solve") },
  ]);
  const mg = htmlEl("div", "lu-group", "", top);
  htmlEl("span", "lu-count", "matrix:", mg);
  const mItems = custom
    ? [{ label: "given", key: "custom" }]
    : [{ label: "C3 · maize · fish", key: "mixing" }, { label: "zero pivot (two sources at −26 ‰)", key: "zeropivot" }];
  const mBtns = presets(mg, mItems.map(it => ({ label: it.label, apply: () => setMatrix(it.key) })));

  const facBox = htmlEl("div", "", "", host);
  const solBox = htmlEl("div", "", "", host);
  for (const b of [facBox, solBox]) Object.assign(b.style, { display: "flex", flexDirection: "column", gap: ".35em", width: "100%" });

  // ---- tab 1: factorise ----
  const nav = htmlEl("div", "lu-row", "", facBox);
  const navg = htmlEl("div", "lu-group", "", nav);
  const bBack = button(navg, "◀ back", () => go(st.cur - 1));
  const bNext = button(navg, "next ▶", () => go(st.cur + 1));
  button(navg, "run to end", () => go(st.steps.length - 1));
  button(navg, "restart", () => go(0));
  const count = htmlEl("span", "lu-count", "", navg);
  const facStage = htmlEl("div", "", "", facBox);
  Object.assign(facStage.style, { display: "flex", flexDirection: "column", gap: ".35em" });
  const swapRow = htmlEl("div", "lu-row", "", facBox);

  // ---- tab 2: solve ----
  const sRow = htmlEl("div", "lu-row", "", solBox);
  const sGroup = htmlEl("div", "lu-group", "", sRow);
  const solveArea = htmlEl("div", "lu-solve", "", solBox);
  const left = htmlEl("div", "lu-left", "", solveArea);
  const right = htmlEl("div", "lu-right", "", solveArea);
  const W = 540, HT = 215, HC = 235;
  const svgT = svgEl("svg", { viewBox: `0 0 ${W} ${HT}`, role: "img", "aria-label": "Sources and samples in isotope space" }, right);
  const nRow = htmlEl("div", "lu-group", "", right);
  htmlEl("span", "lu-count", "count operations for a matrix of size", nRow);
  const nBtns = presets(nRow, [3, 1000].map(n => ({ label: `n = ${n}`, apply: () => { st.nShow = n; drawChart(); } })));
  const svgC = svgEl("svg", { viewBox: `0 0 ${W} ${HC}`, role: "img", "aria-label": "Operation count against the number of samples" }, right);

  // ---------------------------------------------------------------- state changes
  function setMatrix(key: string): void {
    const d = load(key);
    st.key = key; st.A0 = d.A; st.names = d.names; st.B = d.B;
    st.perm = d.A.map((_, i) => i);
    mBtns.forEach((b, i) => b.classList.toggle("active", mItems[i].key === key));
    refactor();
    go(st.tab === "solve" ? st.steps.length - 1 : 0);
  }
  function refactor(): void {
    const A = st.perm.map(i => st.A0[i]);
    st.steps = eliminate(A);
    st.solved = []; st.current = -1;
    buildSampleButtons();
    if (st.tab === "solve" && factors()) solveSample(st.B[0], "#1");
    drawSolve();
  }
  function go(k: number): void {
    st.cur = Math.max(0, Math.min(st.steps.length - 1, k));
    drawFactor();
  }
  function showTab(t: "factor" | "solve"): void {
    st.tab = t;
    facBox.style.display = t === "factor" ? "flex" : "none";
    solBox.style.display = t === "solve" ? "flex" : "none";
    tabBtns.forEach((b, i) => b.classList.toggle("active", (i === 0) === (t === "factor")));
    if (t === "solve" && st.solved.length === 0 && factors()) solveSample(st.B[0], "#1");
    drawSolve();
  }
  const last = () => st.steps[st.steps.length - 1];
  const factors = (): { L: M; U: M } | null => {
    const s = last();
    return s.kind === "done" && s.U.every((r, i) => !r[i].isZero()) ? { L: s.L, U: s.U } : null;
  };

  // ---------------------------------------------------------------- tab 1 drawing
  const n = () => st.A0.length;
  const Rn = (i: number) => `R${i + 1}`;
  const lname = (i: number, k: number) => `ℓ${sub(String(i + 1) + String(k + 1))}`;
  const uname = (i: number, k: number) => `u${sub(String(i + 1) + String(k + 1))}`;

  function drawFactor(): void {
    const s = st.steps[st.cur], N = n();
    bBack.disabled = st.cur === 0;
    bNext.disabled = st.cur === st.steps.length - 1;
    count.textContent = `step ${st.cur + 1} of ${st.steps.length}`;
    // Which multipliers are known: everything filled in up to this step.
    const known = new Set<string>();
    for (let t = 1; t <= st.cur; t++) { const q = st.steps[t]; if (q.kind === "elim") known.add(`${q.i},${q.k}`); }
    const Ls = s.L.map((r, i) => r.map((v, j) => fr(v)));
    const lCell: CellCls = (i, j) => {
      if (j >= i) return "";
      if (s.kind === "elim" && i === s.i && j === s.k) return "nw";
      return known.has(`${i},${j}`) ? "" : "faint";
    };
    let Us = s.U.map(r => r.map(v => fr(v)));
    let uCell: CellCls = () => "", uRow: (i: number) => string = () => "";
    let title = "", text = "";
    const A = st.perm.map(i => st.A0[i]);
    const prod = matMul(s.L, s.U);
    const same = prod.every((r, i) => r.every((v, j) => v.eq(A[i][j])));

    if (s.kind === "start") {
      title = "Start: L = I, U = A";
      text = `<p>U starts as a copy of A. L starts as the identity; the grey 0s below its diagonal are the multipliers we do not know yet. `
        + `The right-hand side <b>b</b> is not used anywhere in this tab.</p>`
        + `<p>Click <b>next ▶</b>: each step computes one multiplier ℓ = (entry to remove) ÷ (pivot), writes it into L, and subtracts ℓ × the pivot row in U.</p>`;
      uCell = (i, j) => (i === 0 && j === 0 ? "pv" : "");
    } else if (s.kind === "elim") {
      const { k, i, l, before } = s;
      const lt = fr(l), neg = l.sign() < 0;
      title = `Step: ${lname(i, k)} = ${uname(i, k)} ÷ ${uname(k, k)} = ${par(fr(before[i][k]), before[i][k])} ÷ ${par(fr(before[k][k]), before[k][k])} = ${lt}`;
      const rowB = vecTxt(before[i], fr), rowK = vecTxt(before[k], fr), rowA = vecTxt(s.U[i], fr);
      text = `<p>Pivot ${uname(k, k)} = ${fr(before[k][k])} (circled). The multiplier goes into L at row ${i + 1}, column ${k + 1}: the same spot where U gets its new 0 (green in both).</p>`
        + `<p class="lu-eq">${Rn(i)} ← ${Rn(i)} − ${par(lt, l)}·${Rn(k)}:&nbsp; ${rowB} ${neg ? "+" : "−"} ${fr(l.abs())}·${rowK} = ${rowA}</p>`
        + `<p>Operations: 1 division + ${N - k - 1} multiplication${N - k - 1 === 1 ? "" : "s"} + ${N - k - 1} subtraction${N - k - 1 === 1 ? "" : "s"} = ${1 + 2 * (N - k - 1)}. `
        + `So far: <b>${s.ops}</b>.</p>`;
      uCell = (r, c) => (r === i && c === k ? "nw" : r === k && c === k ? "pv" : "");
      uRow = r => (r === i ? "chg" : r === k ? "src" : "");
    } else if (s.kind === "zero") {
      const { k, i } = s;
      const below = s.U[i][k];
      if (s.singular) {
        title = `<span class="lu-no">Pivot ${uname(k, k)} = 0, and everything below it is 0 too</span>`;
        text = `<p>No row has a non-zero entry in column ${k + 1} from row ${k + 1} down. No reordering of the rows helps: A is <b>singular</b>, `
          + `two sources cannot be told apart.</p>`;
        uCell = (r, c) => (r === k && c === k ? "zp" : "");
      } else {
        title = `<span class="lu-no">Zero pivot: ${lname(i, k)} = ${fr(below)} ÷ 0</span>`;
        // What lu_nopivot stores in doubles: l = x/0 = ±Inf, then row i: 0 - Inf*0 = NaN, u - Inf*v = ±Inf.
        const inf = below.sign() > 0 ? 1 : -1;
        Us = s.U.map(r => r.map(v => fr(v)));
        for (let j = k; j < N; j++) {
          const v = s.U[k][j];
          Us[i][j] = v.isZero() ? "NaN" : (inf * v.sign() > 0 ? "−Inf" : "Inf");
        }
        text = `<p>The pivot ${uname(k, k)} is <b>0</b>, but ${uname(i, k)} = ${fr(below)} still has to be removed. The multiplier would be ${fr(below)} ÷ 0: `
          + `elimination without row swaps <b>stops here</b>. In doubles there is no error message: ${lname(i, k)} = ${below.sign() > 0 ? "" : "−"}Inf, and row ${i + 1} becomes NaN and Inf (red), as in the lecture.</p>`
          + `<p>The system itself is fine: ${Rn(i)} has a non-zero number in column ${k + 1}. Swap the rows and start again.</p>`;
        uCell = (r, c) => (r === k && c === k ? "zp" : r === i && c >= k ? "bad" : "");
        uRow = r => (r === i ? "chg" : "");
      }
    } else {
      title = `Done: A = LU, after ${s.ops} operations`;
      const lst = [] as string[];
      for (let i = 1; i < N; i++) for (let k = 0; k < i; k++) lst.push(`${lname(i, k)} = ${fr(s.L[i][k])}`);
      text = `<p>L is unit lower triangular and holds every multiplier (${lst.join(", ")}); U is upper triangular. `
        + `Nothing was thrown away: L·U gives back A exactly.</p>`
        + `<p><b>b</b> never entered. Store L and U, and every new sample costs only forward and back substitution: tab <b>2. Solve for samples</b>.</p>`;
      uCell = (i, j) => (i === j ? "pv" : "");
    }
    const permNote = st.perm.some((p, i) => p !== i) ? ` (rows reordered: ${st.perm.map(p => Rn(p)).join(", ")} of the original)` : "";
    const prodS = s.kind === "zero" && !s.singular
      ? `<div class="lu-cap" style="max-width:9em">no product: L and U contain Inf and NaN</div>`
      : matrix(prod.map(r => r.map(v => fr(v))), `<b>L · U</b> ${same ? '<span class="lu-ok">= A ✓</span>' : '<span class="lu-no">≠ A</span>'}`);
    facStage.innerHTML =
      `<div class="lu-title">${title}</div>`
      + `<div class="lu-mats">`
      + matrix(A.map(r => r.map(v => fr(v))), `<b>A</b>${permNote}`, () => "", () => "", true)
      + `<span class="lu-op">=</span>`
      + matrix(Ls, "<b>L</b>  (multipliers)", lCell)
      + `<span class="lu-op">·</span>`
      + matrix(Us, "<b>U</b>  (eliminated A)", uCell, uRow, true)
      + `<span class="lu-op" style="margin-left:.6em">check:</span>`
      + prodS
      + `</div>`
      + `<div class="lu-text">${text}</div>`
      + `<div class="lu-key"><span class="k-pv">pivot</span><span class="k-src">pivot row</span><span class="k-chg">row changed</span><span class="k-nw">new multiplier in L = new 0 in U</span></div>`;

    swapRow.innerHTML = "";
    if (s.kind === "zero" && !s.singular) {
      const k = s.k;
      const U = s.U;
      let r = -1;
      for (let i = k + 1; i < N; i++) if (!U[i][k].isZero()) { r = i; break; }
      button(swapRow, `swap ${Rn(k)} ↔ ${Rn(r)} and factorise again`, () => {
        [st.perm[k], st.perm[r]] = [st.perm[r], st.perm[k]];
        refactor(); go(0);
      }, "active");
      htmlEl("span", "lu-count", "(that row order is what pivoting chooses; next unit)", swapRow);
    }
  }

  // ---------------------------------------------------------------- tab 2
  function buildSampleButtons(): void {
    sGroup.innerHTML = "";
    htmlEl("span", "lu-count", "sample:", sGroup);
    st.B.forEach((b, j) => button(sGroup, `consumer #${j + 1}`, () => { solveSample(b, `#${j + 1}`); drawSolve(); }));
    button(sGroup, "solve another sample ▶", () => { newSample(); drawSolve(); }, "active");
    button(sGroup, "reset count", () => { st.solved = []; st.current = -1; st.r = rng(seed); drawSolve(); });
  }
  function solveSample(b: Frac[], label: string): void {
    const F = factors();
    if (!F) return;
    st.solved.push({ b, x: solveLU(F.L, F.U, st.perm.map(i => b[i])).x, label });
    st.current = st.solved.length - 1;
  }
  function newSample(): void {
    // A diet with fractions in steps of 0.05, each at least 0.1: a consumer inside the mixing polygon.
    const N = n();
    let f: number[];
    do {
      f = Array.from({ length: N - 1 }, () => 2 + Math.floor(st.r.unif() * 14));
      f.push(20 - f.reduce((a, v) => a + v, 0));
    } while (f[N - 1] < 2);
    const fx = f.map(v => Frac.of(v, 20));
    const b = st.A0.map(row => row.reduce((s, a, j) => s.add(a.mul(fx[j])), ZERO));
    solveSample(b, `new ${st.solved.filter(s => s.label.startsWith("new")).length + 1}`);
  }
  function solveLU(L: M, U: M, b: Frac[]): { y: Frac[]; x: Frac[] } {
    const N = b.length, y: Frac[] = [], x: Frac[] = Array(N).fill(ZERO);
    for (let i = 0; i < N; i++) { let s = b[i]; for (let j = 0; j < i; j++) s = s.sub(L[i][j].mul(y[j])); y.push(s); }
    for (let i = N - 1; i >= 0; i--) { let s = y[i]; for (let j = i + 1; j < N; j++) s = s.sub(U[i][j].mul(x[j])); x[i] = s.div(U[i][i]); }
    return { y, x };
  }

  function drawSolve(): void {
    nBtns.forEach((b, i) => b.classList.toggle("active", [3, 1000][i] === st.nShow));
    const F = factors(), N = n();
    if (!F) {
      left.innerHTML = `<div class="lu-title lu-no">No factors to reuse</div><div class="lu-text"><p>Elimination stopped at a zero pivot `
        + `(tab 1). Without U there is nothing to substitute with. Go back to tab 1 and swap the rows first.</p></div>`;
      sGroup.querySelectorAll("button").forEach(b => (b.disabled = true));
      drawTriangle(); drawChart();
      return;
    }
    sGroup.querySelectorAll("button").forEach(b => (b.disabled = false));
    const s = st.solved[st.current];
    if (!s) { left.innerHTML = `<div class="lu-text"><p>Pick a sample.</p></div>`; drawTriangle(); drawChart(); return; }
    const pb = st.perm.map(i => s.b[i]);
    const { y, x } = solveLU(F.L, F.U, pb);
    const ys = (i: number) => `y${sub(String(i + 1))}`, xs = (i: number) => `f${sub(String(i + 1))}`;
    const bs = (i: number) => `b${sub(String(i + 1))}`;
    const term = (c: Frac, v: Frac) => `${par(fr(c), c)}·${par(dec(v), v)}`;
    let h = `<div class="lu-title">Consumer ${s.label}: b = ${vecTxt(s.b)}</div>`;
    if (st.perm.some((p, i) => p !== i)) h += `<div class="lu-count">rows reordered like A: P b = ${vecTxt(pb)}</div>`;
    h += `<div class="lu-mats">` + matrix(F.L.map(r => r.map(v => fr(v))), "<b>L</b> (kept)") + matrix(F.U.map(r => r.map(v => fr(v))), "<b>U</b> (kept)") + `</div>`;
    h += `<div class="lu-eq"><div class="h">Forward substitution L y = b, top to bottom:</div>`;
    for (let i = 0; i < N; i++) {
      const ks = [...Array(i).keys()];
      h += `<div>${ys(i)} = ${bs(i)}${ks.map(j => ` − ${lname(i, j)}${ys(j)}`).join("")} = ${dec(pb[i])}${ks.map(j => ` − ${term(F.L[i][j], y[j])}`).join("")}`
        + `${ks.length ? ` = <b>${dec(y[i])}</b>` : ""}</div>`;
    }
    h += `<div class="h">Back substitution U f = y, bottom to top:</div>`;
    for (let i = N - 1; i >= 0; i--) {
      const ks = [...Array(N - i - 1).keys()].map(t => i + 1 + t);
      const num = ks.length ? `(${dec(y[i])}${ks.map(j => ` − ${term(F.U[i][j], x[j])}`).join("")})` : dec(y[i]);
      h += `<div>${xs(i)} = ${num} ÷ ${par(fr(F.U[i][i]), F.U[i][i])} = <b>${dec(x[i])}</b></div>`;
    }
    h += `</div>`;
    const outside = x.some(v => v.sign() < 0);
    const isMix = st.A0[0].every(v => v.isOne());
    if (isMix) {
      h += `<div>Diet: ${x.map((v, j) => `${st.names[j]} <b>${pctTxt(v)}</b>`).join(", ")}.`
        + (outside ? ` <span class="lu-no">A negative fraction: this consumer lies outside the source triangle (a model problem, not a numerical one).</span>` : "")
        + `</div>`;
    }
    h += `<div class="lu-count">This sample: ${substOps(N)} operations (forward + back). Re-eliminating A would have cost ${factorOps(N)} more.</div>`;
    left.innerHTML = h;
    drawTriangle(); drawChart();
  }
  const pctTxt = (v: Frac) => `${(100 * Number(v.n) / Number(v.d)).toFixed(1).replace(/\.0$/, "")} %`.replace("-", "−");

  // ---- isotope space: sources (columns of A, rows 2 and 3) and solved samples ----
  function drawTriangle(): void {
    svgT.replaceChildren();
    const N = n(), A = st.A0;
    if (N !== 3 || !A[0].every(v => v.isOne())) {
      svgEl("text", { x: 10, y: 30, class: "sm" }, svgT).textContent = "(isotope plot only for 3 sources with a sum-to-1 row)";
      return;
    }
    const num = (f: Frac) => Number(f.n) / Number(f.d);
    const src = [0, 1, 2].map(j => [num(A[1][j]), num(A[2][j])]);
    const pts = st.solved.map(s => [num(s.b[1]), num(s.b[2])]);
    const xsA = [...src, ...pts].map(p => p[0]), ysA = [...src, ...pts].map(p => p[1]);
    const x0 = Math.floor(Math.min(...xsA) - 2), x1 = Math.ceil(Math.max(...xsA) + 2);
    const y0 = Math.floor(Math.min(...ysA) - 1.5), y1 = Math.ceil(Math.max(...ysA) + 1.5);
    const L = 56, R = W - 20, T = 40, B = HT - 40;
    const sx = scale(x0, x1, L, R), sy = scale(y0, y1, B, T);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgT).textContent = "Sources (■) and solved samples (●)";
    const tick = (lo: number, hi: number, most: number) => {
      const st2 = [1, 2, 5, 10, 20].find(d => (hi - lo) / d <= most) ?? 50;
      const t = []; for (let v = Math.ceil(lo / st2) * st2; v <= hi; v += st2) t.push(v); return t;
    };
    for (const v of tick(x0, x1, 7)) {
      svgEl("line", { x1: sx(v), x2: sx(v), y1: T, y2: B, class: "grid" }, svgT);
      svgEl("text", { x: sx(v), y: B + 18, "text-anchor": "middle", class: "sm" }, svgT).textContent = String(v).replace("-", "−");
    }
    for (const v of tick(y0, y1, 4)) {
      svgEl("line", { x1: L, x2: R, y1: sy(v), y2: sy(v), class: "grid" }, svgT);
      svgEl("text", { x: L - 6, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svgT).textContent = String(v).replace("-", "−");
    }
    svgEl("text", { x: (L + R) / 2, y: HT - 2, "text-anchor": "middle", class: "sm" }, svgT).textContent = "δ¹³C (‰)";
    const yl = svgEl("text", { x: 16, y: (T + B) / 2, "text-anchor": "middle", class: "sm", transform: `rotate(-90 16 ${(T + B) / 2})` }, svgT);
    yl.textContent = "δ¹⁵N";
    svgEl("polygon", { points: src.map(p => `${sx(p[0])},${sy(p[1])}`).join(" "), fill: COL.sky, "fill-opacity": 0.25, stroke: COL.blue, "stroke-width": 2 }, svgT);
    src.forEach((p, j) => {
      svgEl("rect", { x: sx(p[0]) - 7, y: sy(p[1]) - 7, width: 14, height: 14, fill: COL.blue }, svgT);
      const right = sx(p[0]) < (L + R) / 2;
      const t = svgEl("text", { x: sx(p[0]) + (right ? 11 : -11), y: sy(p[1]) + (p[1] > (y0 + y1) / 2 ? -6 : 18), "text-anchor": right ? "start" : "end", "font-weight": 600 }, svgT);
      t.textContent = st.names[j];
      t.style.fill = COL.blue;
    });
    st.solved.forEach((s, k) => {
      const out = s.x.some(v => v.sign() < 0);
      const cur = k === st.current;
      svgEl("circle", { cx: sx(pts[k][0]), cy: sy(pts[k][1]), r: cur ? 8 : 5.5, fill: out ? COL.red : COL.orange, stroke: cur ? "#222" : "none", "stroke-width": 2 }, svgT);
    });
  }

  // ---- operation counter: factor once vs re-eliminate every sample ----
  function drawChart(): void {
    svgC.replaceChildren();
    const nn = st.nShow, Fo = factorOps(nn), So = substOps(nn), m = st.solved.length;
    const M = Math.max(10, m);
    const once = (k: number) => (k === 0 ? 0 : Fo + k * So), again = (k: number) => k * (Fo + So);
    const L = 84, R = W - 16, T = 92, B = HC - 34;
    const raw = again(M) / 4, mag = 10 ** Math.floor(Math.log10(raw));
    const unit = [1, 2, 2.5, 5, 10].map(d => d * mag).find(d => d >= raw)!;
    const ymax = 4 * unit;
    const sx = scale(0, M, L, R), sy = scale(0, ymax, B, T);
    svgEl("text", { x: 8, y: 22, class: "ttl" }, svgC).textContent = `Operations after m samples, n = ${nn === 3 ? "3 (this A)" : "1000"}`;
    const t1 = svgEl("text", { x: 8, y: 44 }, svgC);
    t1.innerHTML = `<tspan fill="${COL.red}" font-weight="700">re-eliminate each time: ${sci(again(m))}</tspan>`;
    const t2 = svgEl("text", { x: 8, y: 66 }, svgC);
    const ratio = m ? again(m) / once(m) : 1;
    t2.innerHTML = `<tspan fill="${COL.blue}" font-weight="700">factor once: ${sci(once(m))}</tspan>`
      + `<tspan fill="#555">  m = ${m}${m ? `, ${ratio.toFixed(ratio < 10 ? 2 : 1)}× fewer` : ""}</tspan>`;
    for (let k = 0; k <= 4; k++) {
      const v = unit * k;
      svgEl("line", { x1: L, x2: R, y1: sy(v), y2: sy(v), class: "grid" }, svgC);
      svgEl("text", { x: L - 6, y: sy(v) + 6, "text-anchor": "end", class: "sm" }, svgC).textContent = sci(v);
    }
    const step = M <= 10 ? 1 : M <= 20 ? 2 : M <= 50 ? 5 : 10;
    for (let k = 0; k <= M; k += step) svgEl("text", { x: sx(k), y: B + 20, "text-anchor": "middle", class: "sm" }, svgC).textContent = String(k);
    svgEl("line", { x1: L, x2: R, y1: B, y2: B, class: "axis" }, svgC);
    svgEl("text", { x: R, y: HC - 0, "text-anchor": "end", class: "sm" }, svgC).textContent = "samples m";
    const path = (f: (k: number) => number) => [...Array(M + 1).keys()].map(k => `${k ? "L" : "M"}${sx(k).toFixed(1)},${sy(f(k)).toFixed(1)}`).join("");
    svgEl("path", { d: path(again), fill: "none", stroke: COL.red, "stroke-width": 3 }, svgC);
    svgEl("path", { d: path(once), fill: "none", stroke: COL.blue, "stroke-width": 3 }, svgC);
    if (m > 0) {
      svgEl("circle", { cx: sx(m), cy: sy(again(m)), r: 6, fill: COL.red }, svgC);
      svgEl("circle", { cx: sx(m), cy: sy(once(m)), r: 6, fill: COL.blue }, svgC);
    }
  }

  // ---------------------------------------------------------------- start
  setMatrix(startKey);
  showTab(st.tab);
}

autoMount(NAME, mount);


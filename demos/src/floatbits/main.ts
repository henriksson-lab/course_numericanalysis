// IEEE double inspector. Type a number or a small expression (or click a preset) and see
// the 64 bits of the stored double: sign | 11 exponent bits | 52 fraction bits. Clicking a
// bit toggles it. Readouts: the exact decimal value stored (BigInt arithmetic), 17 digits,
// shortest round-trip print, the unbiased exponent, the kind (normal/subnormal/zero/Inf/NaN),
// the gap to the next double (next bit pattern) and the relative gap. A number line shows
// the doubles around x and, as an orange marker, the exact value before rounding (the typed
// decimal, or the exact result of the last operation on the stored operands).
// Auto-mounts on <div data-demo="floatbits">. Optional: data-x="0.1" (initial expression).
import { svgEl, htmlEl, presets, button, releaseFocus, autoMount } from "../lib/ui";
import {
  bitsOf, fromBits, fields, kindOf, bitLength, ratOfBits, decimalOf, fmtDecimal,
  evaluate, nextUp, nextDown, sub, div, toNum, rat, cmp, absR, Rat,
} from "./float";

const NAME = "floatbits";
// Okabe-Ito, as in fig01-double-bits.svg: sign red, exponent sky blue, fraction yellow.
const COL = {
  sign: "#D55E00", exp: "#56B4E9", frac: "#F0E442",
  stored: "#0072B2", exact: "#E69F00", tick: "#777",
};
const PALE: Record<string, string> = { "#D55E00": "#f6d9c7", "#56B4E9": "#ddf0fb", "#F0E442": "#fbf8d2" };
const STYLE = `
[data-demo="${NAME}"] { flex-direction: column; gap: .45em; max-width: 1150px; }
[data-demo="${NAME}"] .fb-top { display: flex; flex-wrap: wrap; gap: .4em .8em; align-items: center; }
[data-demo="${NAME}"] .fb-top input[type=text] {
  font: inherit; font-family: ui-monospace, monospace; width: 13em; padding: .15em .4em;
  border: 1px solid #999; border-radius: 6px; }
[data-demo="${NAME}"] .fb-err { color: #b00; font-size: .9em; }
[data-demo="${NAME}"] svg { width: 100%; height: auto; }
[data-demo="${NAME}"] svg text { font-size: 18px; }
[data-demo="${NAME}"] svg .mono { font-family: ui-monospace, monospace; }
[data-demo="${NAME}"] svg .bit { cursor: pointer; }
[data-demo="${NAME}"] .fb-bottom { display: flex; flex-wrap: wrap; gap: .6em 1.2em; width: 100%; }
[data-demo="${NAME}"] .fb-read { flex: 1 1 520px; min-width: 0; display: grid;
  grid-template-columns: max-content 1fr; gap: .1em .8em; font-size: 16px; align-content: start; }
[data-demo="${NAME}"] .fb-read .k { color: #555; white-space: nowrap; }
[data-demo="${NAME}"] .fb-read .v { font-family: ui-monospace, monospace; overflow-wrap: anywhere; }
[data-demo="${NAME}"] .fb-line { flex: 1 1 480px; min-width: 0; }
[data-demo="${NAME}"] .fb-line .qb-note { font-size: 15px; color: #444; }
.reveal [data-demo="${NAME}"] { font-size: 18px; margin: 0; }
.reveal [data-demo="${NAME}"] .fb-read { font-size: 18px; }
.reveal [data-demo="${NAME}"] .fb-line .qb-note { font-size: 17px; }
.reveal [data-demo="${NAME}"] .qb-btn { font-size: 17px; padding: .15em .55em; }
`;

function injectStyle(): void {
  if (document.getElementById("qb-style-" + NAME)) return;
  const s = document.createElement("style");
  s.id = "qb-style-" + NAME;
  s.textContent = STYLE;
  document.head.appendChild(s);
}

const PRESETS = ["0.1", "0.2", "0.1+0.2", "0.3", "1", "1+ε", "2^53", "2^53+1", "1e308*10", "5e-324", "1/3", "248956422"];
const sup = (k: number) => String(k).replace("-", "−");
const pow2 = (k: number) => "2^" + sup(k);
const sci = (v: number, d = 3) => (Math.abs(v) >= 1e-3 && Math.abs(v) < 1e6 ? v.toPrecision(d) : v.toExponential(d - 1).replace("e+", "e"));
const ratFmt = (r: Rat, maxSig = 60) => fmtDecimal(decimalOf(r, maxSig + 5), maxSig);

function mount(host: HTMLElement): void {
  injectStyle();
  const st: { bits: bigint; ex: Rat | null; what: string; src: string } = { bits: 0n, ex: null, what: "", src: "" };

  // ---- controls ----
  const top = htmlEl("div", "fb-top", "", host);
  htmlEl("span", "qb-big", "x =", top);
  const inp = htmlEl("input", "", "", top);
  inp.type = "text";
  inp.spellcheck = false;
  inp.setAttribute("aria-label", "number or expression");
  releaseFocus(inp);
  button(top, "show", () => setExpr(inp.value));
  const err = htmlEl("span", "fb-err", "", top);
  htmlEl("span", "qb-note", "numbers, + − * / ^, eps; or click a bit to flip it", top);
  inp.addEventListener("keydown", ev => { if (ev.key === "Enter") { setExpr(inp.value); inp.blur(); } });
  const pbtns = presets(host, PRESETS.map(p => ({ label: p, apply: () => setExpr(p) })));

  // ---- the 64 bits ----
  const W = 1150, CW = 17, GAP = 10, X0 = 21, Y0 = 30, CH = 40;
  const bitsSvg = svgEl("svg", { viewBox: `0 0 ${W} 112`, role: "img" }, host);
  const cellX = (j: number) => X0 + j * CW + (j >= 1 ? GAP : 0) + (j >= 12 ? GAP : 0);  // j = 0 is the sign bit (bit 63)
  const cells: { r: SVGRectElement; t: SVGTextElement; col: string }[] = [];
  for (let j = 0; j < 64; j++) {
    const col = j === 0 ? COL.sign : j < 12 ? COL.exp : COL.frac;
    const g = svgEl("g", { class: "bit" }, bitsSvg);
    const r = svgEl("rect", { x: cellX(j), y: Y0, width: CW - 1.5, height: CH, rx: 2, stroke: "#fff", "stroke-width": 0 }, g);
    const t = svgEl("text", { x: cellX(j) + (CW - 1.5) / 2, y: Y0 + CH / 2 + 6, "text-anchor": "middle", class: "mono", "pointer-events": "none" }, g);
    g.addEventListener("click", () => {
      st.bits ^= 1n << BigInt(63 - j);
      st.ex = null; st.what = "bits"; st.src = "";
      pbtns.forEach(b => b.classList.remove("active"));
      draw();
    });
    cells.push({ r, t, col });
  }
  const head = (x: number, txt: string, anchor = "middle") =>
    svgEl("text", { x, y: Y0 - 9, "text-anchor": anchor, "font-weight": 600 }, bitsSvg).textContent = txt;
  head(cellX(0) - 2, "s", "start");
  head((cellX(1) + cellX(11) + CW) / 2, "exponent (11 bits)");
  head((cellX(12) + cellX(63) + CW) / 2, "fraction f (52 bits)");
  const lblS = svgEl("text", { x: cellX(0) + 7, y: Y0 + CH + 26, "text-anchor": "middle", "font-weight": 600 }, bitsSvg);
  const lblE = svgEl("text", { x: (cellX(1) + cellX(11) + CW) / 2, y: Y0 + CH + 26, "text-anchor": "middle", class: "mono" }, bitsSvg);
  const lblF = svgEl("text", { x: (cellX(12) + cellX(63) + CW) / 2, y: Y0 + CH + 26, "text-anchor": "middle", class: "mono" }, bitsSvg);

  // ---- readouts and number line ----
  const bottom = htmlEl("div", "fb-bottom", "", host);
  const read = htmlEl("div", "fb-read", "", bottom);
  const rows: Record<string, HTMLElement> = {};
  const rowK: Record<string, HTMLElement> = {};
  for (const [key, label] of [
    ["exact", "before rounding"], ["stored", "stored exactly"], ["p17", "printed, 17 digits"],
    ["short", "shortest print"], ["form", "(−1)^s × 1.f × 2^e"], ["kind", "kind"],
    ["gap", "gap to next double"], ["rel", "relative gap"],
  ]) {
    rowK[key] = htmlEl("span", "k", label, read);
    rows[key] = htmlEl("span", "v", "", read);
  }
  const lineBox = htmlEl("div", "fb-line", "", bottom);
  const LW = 560, LH = 214;
  const lsvg = svgEl("svg", { viewBox: `0 0 ${LW} ${LH}`, role: "img" }, lineBox);
  const lineNote = htmlEl("div", "qb-note", "", lineBox);

  function setExpr(src: string): void {
    try {
      const v = evaluate(src);
      st.bits = bitsOf(v.v);
      st.ex = v.ex; st.what = v.what; st.src = src.trim();
      err.textContent = "";
    } catch (e) {
      err.textContent = "cannot read: " + (e as Error).message;
      return;
    }
    pbtns.forEach(b => b.classList.toggle("active", b.textContent === src));
    draw();
  }

  function draw(): void {
    const b = st.bits;
    const x = fromBits(b);
    const f = fields(b);
    const kind = kindOf(f);
    if (st.what === "bits" || document.activeElement !== inp) inp.value = st.what === "bits" ? String(x) : st.src;

    // bits
    for (let j = 0; j < 64; j++) {
      const one = ((b >> BigInt(63 - j)) & 1n) === 1n;
      cells[j].r.setAttribute("fill", one ? cells[j].col : PALE[cells[j].col]);
      cells[j].t.textContent = one ? "1" : "0";
      cells[j].t.setAttribute("font-weight", one ? "700" : "400");
    }
    lblS.textContent = f.s ? "−" : "+";
    const e = f.E === 0 ? -1022 : f.E - 1023;
    lblE.textContent = kind === "inf" || kind === "nan" ? `${f.E} = all ones` :
      f.E === 0 ? `0 → e = −1022` : `${f.E} − 1023 = ${sup(e)}`;
    const mant = (f.E === 0 ? 0 : 1) + Number(f.F) / 2 ** 52;   // exact in a double
    lblF.textContent = kind === "inf" ? "f = 0 → Inf" : kind === "nan" ? "f ≠ 0 → NaN" :
      `${f.E === 0 ? "0.f" : "1.f"} = ${mant.toPrecision(20)}${f.E === 0 ? " (no hidden 1)" : ""}`;

    // readouts
    const r = ratOfBits(b);
    const stored = r ? decimalOf(r) : null;
    if (st.ex) {
      rowK.exact.textContent = st.what === "typed" ? "you typed" : "exact result";
      const same = r && cmp(st.ex, r) === 0;
      rows.exact.textContent = ratFmt(st.ex) + (same ? "  (no rounding needed)" : "");
    } else {
      rowK.exact.textContent = "you typed";
      rows.exact.textContent = st.what === "bits" ? "(set by clicking bits)" : st.src || "—";
    }
    rows.stored.textContent = stored
      ? fmtDecimal(stored, 60) + (stored.digits.length > 60 ? `  (${stored.digits.length} significant digits)` : "")
      : (kind === "nan" ? "NaN" : f.s ? "−Inf" : "Inf");
    rows.p17.textContent = x.toPrecision(17);
    rows.short.textContent = String(x);
    rows.form.textContent = kind === "normal" ? `${f.s ? "−" : "+"}${mant.toPrecision(17)} × ${pow2(e)}` :
      kind === "subnormal" || kind === "zero" ? `${f.s ? "−" : "+"}${mant.toPrecision(17)} × ${pow2(-1022)}` : "—";
    rows.kind.textContent = kind === "normal" ? "normal (53 significant bits)" :
      kind === "subnormal" ? `subnormal (only ${bitLength(f.F)} significant bits)` :
      kind === "zero" ? (f.s ? "−0 (zero)" : "zero") : kind === "inf" ? (f.s ? "−Inf (overflow)" : "Inf (overflow)") : "NaN (not a number)";
    if (kind === "inf" || kind === "nan") {
      rows.gap.textContent = "—"; rows.rel.textContent = "—";
    } else {
      const gk = (f.E === 0 ? 1 : f.E) - 1075;                 // next pattern is x ± 2^gk
      const nxt = fromBits(b + 1n);
      rows.gap.textContent = `${pow2(gk)} ≈ ${sci(2 ** gk)}` + (Number.isFinite(nxt) ? "" : "  (next pattern is Inf)");
      if (kind === "zero") rows.rel.textContent = "—";
      else {
        const rel = 2 ** gk / Math.abs(x);
        rows.rel.textContent = `${sci(rel)} = ${(rel / 2 ** -52).toPrecision(3)} ε`;
      }
    }
    drawLine(x, r, kind);
  }

  function drawLine(x: number, r: Rat | null, kind: string): void {
    lsvg.replaceChildren();
    const mid = 120;
    if (!r) {
      const t = svgEl("text", { x: LW / 2, y: 90, "text-anchor": "middle" }, lsvg);
      t.textContent = kind === "nan" ? "NaN is not on the number line" : "Inf: beyond the largest double";
      const t2 = svgEl("text", { x: LW / 2, y: 125, "text-anchor": "middle", class: "mono" }, lsvg);
      t2.textContent = "largest double ≈ 1.7976931348623157e308";
      lineNote.textContent = st.ex && kind === "inf"
        ? `The exact result ${ratFmt(st.ex, 6)} is larger than the largest double, so it overflows to Inf (no warning).`
        : kind === "nan" ? "Exponent all ones with a nonzero fraction: NaN." : "Exponent all ones, fraction 0: Inf.";
      return;
    }
    // doubles k = -3..3 around x
    const vals: { k: number; v: number }[] = [{ k: 0, v: x }];
    let lo = x, hi = x;
    for (let k = 1; k <= 3; k++) {
      lo = nextDown(lo); hi = nextUp(hi);
      if (Number.isFinite(lo)) vals.push({ k: -k, v: lo });
      if (Number.isFinite(hi)) vals.push({ k, v: hi });
    }
    vals.sort((a, b) => a.k - b.k);
    const vmin = vals[0].v, vmax = vals[vals.length - 1].v;
    const span = vmax - vmin;                                  // exact: neighbours are close
    const L = 30, R = LW - 30;
    const px = (off: number) => L + (off / span) * (R - L);    // off = v - vmin
    const pos = (v: number) => px(v - vmin);
    svgEl("line", { x1: L - 20, x2: R + 20, y1: mid, y2: mid, stroke: "#444", "stroke-width": 2 }, lsvg);
    for (const { k, v } of vals) {
      const X = pos(v);
      svgEl("line", { x1: X, x2: X, y1: mid - (k ? 16 : 26), y2: mid + (k ? 16 : 26), stroke: k ? COL.tick : COL.stored, "stroke-width": k ? 2.5 : 5 }, lsvg);
    }
    // labels: x (top row), next up (second row), next down (below the line)
    const label = (v: number, y: number, col: string, bold: boolean) => {
      const txt = Number.isInteger(v) && Math.abs(v) < 1e21 ? v.toFixed(0) : v.toPrecision(20).replace(/e\+?/, "e");
      const hw = txt.length * 5.6;
      const X = Math.min(Math.max(pos(v), hw + 2), LW - hw - 2);
      svgEl("line", { x1: pos(v), x2: X, y1: y < mid ? mid - (bold ? 26 : 16) : mid + 16, y2: y < mid ? y + 6 : y - 18, stroke: col, "stroke-width": 1, "stroke-dasharray": "3 3" }, lsvg);
      const t = svgEl("text", { x: X, y, "text-anchor": "middle", class: "mono", fill: col, "font-weight": bold ? 700 : 400 }, lsvg);
      t.style.fill = col;
      t.textContent = txt;
    };
    label(x, 22, COL.stored, true);
    const up = vals.find(o => o.k === 1), dn = vals.find(o => o.k === -1);
    if (up) label(up.v, 62, "#555", false);
    if (dn) label(dn.v, 182, "#555", false);
    // exact value before rounding
    let note = "Blue: the stored double x. Grey: its neighbours (3 on each side).";
    if (st.ex) {
      const vminR = ratOfBits(bitsOf(vmin))!;
      const frac = toNum(div(sub(st.ex, vminR), ratOfBits(bitsOf(span))!));
      const X = L + frac * (R - L);
      if (X > L - 25 && X < R + 25) {
        svgEl("path", { d: `M${X},${mid + 4} l-9,17 h18 z`, fill: COL.exact }, lsvg);
      }
      const d = sub(st.ex, r);
      const sgn = cmp(d, rat(0n));
      if (sgn === 0) note += ` Orange ▲: ${st.what === "typed" ? "the number you typed" : "the exact result"} is itself a double: no rounding.`;
      else {
        const gap = sgn > 0 ? (nextUp(x) - x) : (x - nextDown(x));
        const inGaps = Math.abs(toNum(div(d, ratOfBits(bitsOf(gap))!)));
        const relErr = toNum(div(absR(d), absR(st.ex)));
        const what = st.what === "typed" ? "the number you typed" : "the exact result of the last operation";
        note += ` Orange ▲: ${what}. It lies ${inGaps.toFixed(3)} gaps ${sgn > 0 ? "above" : "below"} x and was rounded to x` +
          (Math.abs(inGaps - 0.5) < 1e-12 ? " (a tie: rounded to the neighbour whose last bit is 0)" : "") +
          `; relative error ${sci(relErr)}` + (kind === "normal" ? " (≤ ε/2 = 1.11e-16)." : " (subnormal: can be far above ε/2).");
      }
    }
    const gDn = x - nextDown(x), gUp = nextUp(x) - x;
    const gt = svgEl("text", { x: L - 20, y: LH - 5, "text-anchor": "start" }, lsvg);
    gt.textContent = `gap below: ${sci(gDn)}`;
    if (Number.isFinite(gUp)) {
      const gt2 = svgEl("text", { x: R + 20, y: LH - 5, "text-anchor": "end" }, lsvg);
      gt2.textContent = `gap above: ${sci(gUp)}`;
    }
    lineNote.textContent = note;
  }

  setExpr(host.dataset.x ?? "0.1");
}

autoMount(NAME, mount);

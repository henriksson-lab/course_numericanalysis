// Exact arithmetic on IEEE doubles for the floatbits demo: bit patterns (via DataView),
// exact rationals (BigInt), exact decimal expansions, neighbours, and a tiny expression
// parser that also records the exact result of the last operation before rounding.

// ---------------- bits ----------------
const buf = new DataView(new ArrayBuffer(8));

export function bitsOf(x: number): bigint {
  buf.setFloat64(0, x);
  return buf.getBigUint64(0);
}
export function fromBits(b: bigint): number {
  buf.setBigUint64(0, BigInt.asUintN(64, b));
  return buf.getFloat64(0);
}
export interface Fields { s: number; E: number; F: bigint }
export function fields(b: bigint): Fields {
  return { s: Number(b >> 63n), E: Number((b >> 52n) & 0x7ffn), F: b & ((1n << 52n) - 1n) };
}
export type Kind = "normal" | "subnormal" | "zero" | "inf" | "nan";
export function kindOf(f: Fields): Kind {
  if (f.E === 2047) return f.F === 0n ? "inf" : "nan";
  if (f.E === 0) return f.F === 0n ? "zero" : "subnormal";
  return "normal";
}
export function bitLength(n: bigint): number {
  if (n < 0n) n = -n;
  return n === 0n ? 0 : n.toString(2).length;
}

/** Next double away from zero (increment the bit pattern; sign bit untouched). */
export function nextAway(b: bigint): bigint { return b + 1n; }
/** Next double towards +Inf / -Inf. */
export function nextUp(x: number): number {
  if (Number.isNaN(x) || x === Infinity) return x;
  if (x === 0) return 5e-324;
  const b = bitsOf(x);
  return fromBits(x > 0 ? b + 1n : b - 1n);
}
export function nextDown(x: number): number { return -nextUp(-x); }

// ---------------- rationals ----------------
export interface Rat { n: bigint; d: bigint }   // d > 0, reduced
function gcd(a: bigint, b: bigint): bigint {
  if (a < 0n) a = -a;
  while (b) { [a, b] = [b, a % b]; }
  return a;
}
export function rat(n: bigint, d: bigint = 1n): Rat {
  if (d === 0n) throw new Error("division by zero");
  if (d < 0n) { n = -n; d = -d; }
  const g = gcd(n, d);
  return g > 1n ? { n: n / g, d: d / g } : { n, d };
}
export const add = (a: Rat, b: Rat): Rat => rat(a.n * b.d + b.n * a.d, a.d * b.d);
export const sub = (a: Rat, b: Rat): Rat => rat(a.n * b.d - b.n * a.d, a.d * b.d);
export const mul = (a: Rat, b: Rat): Rat => rat(a.n * b.n, a.d * b.d);
export const div = (a: Rat, b: Rat): Rat => rat(a.n * b.d, a.d * b.n);
export const neg = (a: Rat): Rat => ({ n: -a.n, d: a.d });
export const absR = (a: Rat): Rat => ({ n: a.n < 0n ? -a.n : a.n, d: a.d });
export const cmp = (a: Rat, b: Rat): number => { const c = a.n * b.d - b.n * a.d; return c > 0n ? 1 : c < 0n ? -1 : 0; };
export function powR(a: Rat, k: number): Rat {
  const e = BigInt(Math.abs(k));
  const r = rat(a.n ** e, a.d ** e);
  return k >= 0 ? r : div(rat(1n), r);
}
/** Approximate value as a double (fine for ratios of moderate size). */
export function toNum(a: Rat): number {
  let { n, d } = a;
  const sh = Math.max(bitLength(n), bitLength(d)) - 900;
  if (sh > 0) { n >>= BigInt(sh); d >>= BigInt(sh); if (d === 0n) return n === 0n ? 0 : (n > 0n ? Infinity : -Infinity); }
  return Number(n) / Number(d);
}
/** Exact value of a finite double, or null for Inf/NaN. */
export function ratOfBits(b: bigint): Rat | null {
  const { s, E, F } = fields(b);
  if (E === 2047) return null;
  const m = E === 0 ? F : F + (1n << 52n);
  const e = (E === 0 ? 1 : E) - 1075;          // value = m * 2^e
  let r = e >= 0 ? rat(m << BigInt(e)) : rat(m, 1n << BigInt(-e));
  if (s) r = neg(r);
  return r;
}
export const ratOf = (x: number): Rat | null => ratOfBits(bitsOf(x));

/** Exact rational value of a decimal literal like "1.5e-3". */
export function parseDecimal(s: string): Rat {
  const m = /^(\d*)(?:\.(\d*))?(?:[eE]([+-]?\d+))?$/.exec(s);
  if (!m || (m[1] === "" && (m[2] ?? "") === "")) throw new Error("bad number " + s);
  const frac = m[2] ?? "";
  let n = BigInt((m[1] || "0") + frac);
  let e = Number(m[3] ?? 0) - frac.length;
  if (Math.abs(e) > 5000) throw new Error("exponent too large");
  return e >= 0 ? rat(n * 10n ** BigInt(e)) : rat(n, 10n ** BigInt(-e));
}

// ---------------- decimal expansions ----------------
export interface Decimal { neg: boolean; digits: string; E: number; terminates: boolean }
/**
 * Decimal expansion of r: value = ±0.d1d2d3… × 10^(E+1), i.e. d1.d2d3… × 10^E.
 * Exact (all digits) when the expansion terminates; otherwise `maxSig` digits (truncated).
 */
export function decimalOf(r: Rat, maxSig = 80): Decimal {
  const negv = r.n < 0n;
  let n = negv ? -r.n : r.n;
  const d = r.d;
  if (n === 0n) return { neg: false, digits: "0", E: 0, terminates: true };
  let p2 = 0, p5 = 0, t = d;
  while (t % 2n === 0n) { t /= 2n; p2++; }
  while (t % 5n === 0n) { t /= 5n; p5++; }
  if (t === 1n) {
    const k = Math.max(p2, p5);
    const N = n * 10n ** BigInt(k) / d;           // exact integer
    const s = N.toString();
    const E = s.length - 1 - k;
    return { neg: negv, digits: s.replace(/0+$/, "") || "0", E, terminates: true };
  }
  // non-terminating: long division
  const ip = n / d;
  let rem = n % d;
  let s = ip === 0n ? "" : ip.toString();
  let E = s.length - 1;
  let lead = s.length === 0;
  while (s.length < maxSig) {
    rem *= 10n;
    const q = rem / d; rem %= d;
    if (lead) { if (q === 0n) { E--; continue; } lead = false; }
    s += q.toString();
  }
  return { neg: negv, digits: s, E, terminates: false };
}

/** Format a Decimal; at most `maxSig` significant digits are shown, then "…". */
export function fmtDecimal(x: Decimal, maxSig = 60): string {
  let dg = x.digits;
  let tail = "";
  if (dg.length > maxSig) { dg = dg.slice(0, maxSig); tail = "…"; }
  else if (!x.terminates) tail = "…";
  const sgn = x.neg ? "−" : "";
  if (x.E >= -7 && x.E < 21) {
    if (x.E < 0) return sgn + "0." + "0".repeat(-x.E - 1) + dg + tail;
    if (dg.length <= x.E + 1) return sgn + dg + (tail ? tail : "0".repeat(x.E + 1 - dg.length));
    return sgn + dg.slice(0, x.E + 1) + "." + dg.slice(x.E + 1) + tail;
  }
  return sgn + dg[0] + (dg.length > 1 ? "." + dg.slice(1) : "") + tail + "e" + (x.E >= 0 ? "+" : "") + x.E;
}

// ---------------- expressions ----------------
export interface Val { v: number; ex: Rat | null; what: string }

/**
 * Evaluate a small expression in doubles: numbers, eps/ε, Inf, NaN, + - * / ^, parentheses.
 * `ex` is the exact value before the final rounding: the typed literal itself, or the
 * exact result of the last operation applied to the (already stored) operands.
 */
export function evaluate(src: string): Val {
  const s = src.replace(/[−–]/g, "-").replace(/[×·]/g, "*").replace(/\*\*/g, "^").replace(/\s+/g, "");
  let i = 0;
  const peek = () => s[i];
  const fail = (m: string): never => { throw new Error(m); };
  function primary(): Val {
    const c = peek();
    if (c === "(") { i++; const v = sum(); if (peek() !== ")") fail("missing )"); i++; return v; }
    if (c === "-") { i++; const v = power(); return { v: -v.v, ex: v.ex ? neg(v.ex) : null, what: v.what }; }
    if (c === "+") { i++; return power(); }
    const num = /^(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?/.exec(s.slice(i));
    if (num) { i += num[0].length; return { v: Number(num[0]), ex: parseDecimal(num[0]), what: "typed" }; }
    const w = /^(eps|ε|Inf|inf|Infinity|NaN|nan)/.exec(s.slice(i));
    if (w) {
      i += w[0].length;
      if (w[0] === "eps" || w[0] === "ε") return { v: 2 ** -52, ex: rat(1n, 1n << 52n), what: "typed" };
      if (w[0].toLowerCase() === "nan") return { v: NaN, ex: null, what: "typed" };
      return { v: Infinity, ex: null, what: "typed" };
    }
    return fail(i < s.length ? "unexpected '" + c + "'" : "unexpected end");
  }
  function power(): Val {
    const a = primary();
    if (peek() !== "^") return a;
    i++;
    const b = power();                              // right associative
    const v = a.v ** b.v;
    const ra = ratOf(a.v);
    let ex: Rat | null = null;
    if (ra && Number.isInteger(b.v) && Math.abs(b.v) <= 1100 && !(ra.n === 0n && b.v < 0)) ex = powR(ra, b.v);
    return { v, ex, what: "op" };
  }
  function product(): Val {
    let a = power();
    while (peek() === "*" || peek() === "/") {
      const op = s[i++];
      const b = power();
      const ra = ratOf(a.v), rb = ratOf(b.v);
      let ex: Rat | null = null;
      if (ra && rb && !(op === "/" && rb.n === 0n)) ex = op === "*" ? mul(ra, rb) : div(ra, rb);
      a = { v: op === "*" ? a.v * b.v : a.v / b.v, ex, what: "op" };
    }
    return a;
  }
  function sum(): Val {
    let a = product();
    while (peek() === "+" || peek() === "-") {
      const op = s[i++];
      const b = product();
      const ra = ratOf(a.v), rb = ratOf(b.v);
      const ex = ra && rb ? (op === "+" ? add(ra, rb) : sub(ra, rb)) : null;
      a = { v: op === "+" ? a.v + b.v : a.v - b.v, ex, what: "op" };
    }
    return a;
  }
  if (!s) fail("empty");
  const r = sum();
  if (i < s.length) fail("unexpected '" + s[i] + "'");
  return r;
}

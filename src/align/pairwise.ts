import type { SubstitutionMatrix } from "./matrix.ts";

export type AlignMode = "global" | "local";

export interface Scoring {
  matrix: SubstitutionMatrix;
  openGap: number;
  extendGap: number;
  substitute?: string;
}

export interface AlignOptions extends Scoring {
  mode: AlignMode;
  maxCells?: number;
}

export interface AlignmentStats {
  length: number;
  identities: number;
  positives: number;
  mismatches: number;
  gaps: number;
  gapOpens: number;
  identity: number;
}

export interface Alignment {
  mode: AlignMode;
  score: number;
  targetStart: number;
  targetEnd: number;
  queryStart: number;
  queryEnd: number;
  operations: string;
  cigar: string;
  coordinates: [number[], number[]];
  alignedTarget: string;
  alignedQuery: string;
  midline: string;
  stats: AlignmentStats;
}

export interface ScoreResult {
  mode: AlignMode;
  score: number;
  targetEnd: number;
  queryEnd: number;
}

const STATE_M = 0;
const STATE_X = 1;
const STATE_Y = 2;
const STATE_START = 3;

export const DEFAULT_MAX_CELLS = 1 << 28;

export class AlignmentInputError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "AlignmentInputError";
  }
}

export function encodeForMatrix(sequence: string, matrix: SubstitutionMatrix, substitute?: string): Uint8Array {
  const out = new Uint8Array(sequence.length);
  let fallback = -1;
  if (substitute !== undefined) {
    fallback = matrix.index[substitute.charCodeAt(0)] ?? -1;
    if (fallback < 0 || substitute.length !== 1) {
      throw new AlignmentInputError(`substitute symbol '${substitute}' is not in matrix ${matrix.name}`);
    }
  }
  for (let i = 0; i < sequence.length; i++) {
    const code = sequence.charCodeAt(i);
    let index = code < 256 ? matrix.index[code]! : -1;
    if (index < 0) {
      if (fallback < 0) {
        throw new AlignmentInputError(`symbol '${sequence[i]}' at position ${i + 1} is not in matrix ${matrix.name}`);
      }
      index = fallback;
    }
    out[i] = index;
  }
  return out;
}

function validateScoring(options: Scoring): void {
  if (!Number.isFinite(options.openGap) || !Number.isFinite(options.extendGap)) {
    throw new AlignmentInputError("gap scores must be finite numbers");
  }
}

interface DpOutcome {
  score: number;
  endI: number;
  endJ: number;
  endState: number;
  trace: Uint8Array | null;
}

function runDp(t: Uint8Array, q: Uint8Array, options: AlignOptions, withTrace: boolean): DpOutcome {
  const n = t.length;
  const m = q.length;
  const width = m + 1;
  const local = options.mode === "local";
  const open = options.openGap;
  const extend = options.extendGap;
  const size = options.matrix.size;
  const scores = options.matrix.scores;
  const NEG = Number.NEGATIVE_INFINITY;

  let trace: Uint8Array | null = null;
  if (withTrace) {
    const cells = (n + 1) * width;
    const limit = options.maxCells ?? DEFAULT_MAX_CELLS;
    if (cells > limit) {
      throw new RangeError(`alignment needs ${cells} traceback cells, above the limit of ${limit}; use score-only alignment`);
    }
    trace = new Uint8Array(cells);
  }

  let Mp = new Float64Array(width).fill(NEG);
  let Xp = new Float64Array(width).fill(NEG);
  let Yp = new Float64Array(width).fill(NEG);
  let Mc = new Float64Array(width);
  let Xc = new Float64Array(width);
  let Yc = new Float64Array(width);

  let best = 0;
  let bestI = 0;
  let bestJ = 0;

  Mc[0] = local ? NEG : 0;
  Xc[0] = NEG;
  Yc[0] = NEG;
  for (let j = 1; j <= m; j++) {
    Mc[j] = NEG;
    Xc[j] = NEG;
    const fromM = Mc[j - 1]! + open;
    const fromY = Yc[j - 1]! + extend;
    if (fromM >= fromY) {
      Yc[j] = fromM;
      if (trace) trace[j] = STATE_M << 4;
    } else {
      Yc[j] = fromY;
      if (trace) trace[j] = STATE_Y << 4;
    }
  }

  for (let i = 1; i <= n; i++) {
    const tmpM = Mp;
    Mp = Mc;
    Mc = tmpM;
    const tmpX = Xp;
    Xp = Xc;
    Xc = tmpX;
    const tmpY = Yp;
    Yp = Yc;
    Yc = tmpY;

    const rowBase = i * width;
    const scoreRow = t[i - 1]! * size;

    Mc[0] = NEG;
    Yc[0] = NEG;
    {
      const fromM = Mp[0]! + open;
      const fromX = Xp[0]! + extend;
      const fromY = Yp[0]! + open;
      let value = fromM;
      let pointer = STATE_M;
      if (fromX > value) {
        value = fromX;
        pointer = STATE_X;
      }
      if (fromY > value) {
        value = fromY;
        pointer = STATE_Y;
      }
      Xc[0] = value;
      if (trace) trace[rowBase] = pointer << 2;
    }

    for (let j = 1; j <= m; j++) {
      const dM = Mp[j - 1]!;
      const dX = Xp[j - 1]!;
      const dY = Yp[j - 1]!;
      let diag = dM;
      let mPointer = STATE_M;
      if (dX > diag) {
        diag = dX;
        mPointer = STATE_X;
      }
      if (dY > diag) {
        diag = dY;
        mPointer = STATE_Y;
      }
      if (local && diag <= 0) {
        diag = 0;
        mPointer = STATE_START;
      }
      const mValue = diag + scores[scoreRow + q[j - 1]!]!;
      Mc[j] = mValue;

      const uM = Mp[j]! + open;
      const uX = Xp[j]! + extend;
      const uY = Yp[j]! + open;
      let xValue = uM;
      let xPointer = STATE_M;
      if (uX > xValue) {
        xValue = uX;
        xPointer = STATE_X;
      }
      if (uY > xValue) {
        xValue = uY;
        xPointer = STATE_Y;
      }
      Xc[j] = xValue;

      const lM = Mc[j - 1]! + open;
      const lY = Yc[j - 1]! + extend;
      const lX = Xc[j - 1]! + open;
      let yValue = lM;
      let yPointer = STATE_M;
      if (lX > yValue) {
        yValue = lX;
        yPointer = STATE_X;
      }
      if (lY > yValue) {
        yValue = lY;
        yPointer = STATE_Y;
      }
      Yc[j] = yValue;

      if (trace) trace[rowBase + j] = mPointer | (xPointer << 2) | (yPointer << 4);
      if (local && mValue > best) {
        best = mValue;
        bestI = i;
        bestJ = j;
      }
    }
  }

  if (local) {
    return { score: best, endI: bestI, endJ: bestJ, endState: STATE_M, trace };
  }
  let score = Mc[m]!;
  let endState = STATE_M;
  if (Xc[m]! > score) {
    score = Xc[m]!;
    endState = STATE_X;
  }
  if (Yc[m]! > score) {
    score = Yc[m]!;
    endState = STATE_Y;
  }
  if (n === 0 && m === 0) score = 0;
  return { score, endI: n, endJ: m, endState, trace };
}

function traceback(outcome: DpOutcome, width: number, local: boolean): { ops: string; startI: number; startJ: number } {
  const trace = outcome.trace!;
  const ops: string[] = [];
  let i = outcome.endI;
  let j = outcome.endJ;
  let state = outcome.endState;
  if (local && outcome.score <= 0) return { ops: "", startI: 0, startJ: 0 };
  for (;;) {
    if (state === STATE_M) {
      if (!local && i === 0 && j === 0) break;
      const pointer = trace[i * width + j]! & 3;
      ops.push("M");
      i--;
      j--;
      if (pointer === STATE_START) break;
      state = pointer;
    } else if (state === STATE_X) {
      const pointer = (trace[i * width + j]! >> 2) & 3;
      ops.push("D");
      i--;
      state = pointer;
    } else {
      const pointer = (trace[i * width + j]! >> 4) & 3;
      ops.push("I");
      j--;
      state = pointer;
    }
    if (i < 0 || j < 0) throw new Error("traceback left the matrix");
  }
  return { ops: ops.reverse().join(""), startI: i, startJ: j };
}

export function summarize(
  target: string,
  query: string,
  ops: string,
  targetStart: number,
  queryStart: number,
  matrix: SubstitutionMatrix,
): Omit<Alignment, "mode" | "score"> {
  const tParts: string[] = [];
  const qParts: string[] = [];
  const mid: string[] = [];
  const coordT: number[] = [targetStart];
  const coordQ: number[] = [queryStart];
  const cigar: string[] = [];
  let i = targetStart;
  let j = queryStart;
  let identities = 0;
  let positives = 0;
  let mismatches = 0;
  let gaps = 0;
  let gapOpens = 0;
  let runOp = "";
  let runLength = 0;
  for (let k = 0; k < ops.length; k++) {
    const op = ops[k]!;
    if (op !== runOp) {
      if (runLength > 0) {
        cigar.push(`${runLength}${runOp}`);
        coordT.push(i);
        coordQ.push(j);
      }
      if (op !== "M") gapOpens++;
      runOp = op;
      runLength = 0;
    }
    runLength++;
    if (op === "M") {
      const a = target[i]!;
      const b = query[j]!;
      tParts.push(a);
      qParts.push(b);
      if (a.toUpperCase() === b.toUpperCase()) {
        identities++;
        mid.push("|");
      } else {
        mismatches++;
        mid.push(".");
      }
      const ia = matrix.index[a.charCodeAt(0)]!;
      const ib = matrix.index[b.charCodeAt(0)]!;
      if (ia >= 0 && ib >= 0 && matrix.scores[ia * matrix.size + ib]! > 0) positives++;
      i++;
      j++;
    } else if (op === "D") {
      tParts.push(target[i]!);
      qParts.push("-");
      mid.push("-");
      gaps++;
      i++;
    } else {
      tParts.push("-");
      qParts.push(query[j]!);
      mid.push("-");
      gaps++;
      j++;
    }
  }
  if (runLength > 0) {
    cigar.push(`${runLength}${runOp}`);
    coordT.push(i);
    coordQ.push(j);
  }
  if (coordT.length === 1) {
    coordT.push(i);
    coordQ.push(j);
  }
  const length = ops.length;
  return {
    targetStart,
    targetEnd: i,
    queryStart,
    queryEnd: j,
    operations: ops,
    cigar: cigar.join(""),
    coordinates: [coordT, coordQ],
    alignedTarget: tParts.join(""),
    alignedQuery: qParts.join(""),
    midline: mid.join(""),
    stats: {
      length,
      identities,
      positives,
      mismatches,
      gaps,
      gapOpens,
      identity: length === 0 ? 0 : identities / length,
    },
  };
}

const LIN_DIAG = 1;
const LIN_LEFT = 2;
const LIN_UP = 4;
const LIN_START = 8;
const LIN_END = 16;

function linearTraceback(
  t: Uint8Array,
  q: Uint8Array,
  options: AlignOptions,
): { score: number; ops: string; startI: number; startJ: number } {
  const n = t.length;
  const m = q.length;
  const width = m + 1;
  const local = options.mode === "local";
  const gap = options.extendGap;
  const size = options.matrix.size;
  const scores = options.matrix.scores;
  const cells = (n + 1) * width;
  const limit = options.maxCells ?? DEFAULT_MAX_CELLS;
  if (cells > limit) {
    throw new RangeError(`alignment needs ${cells} traceback cells, above the limit of ${limit}; use score-only alignment`);
  }
  const trace = new Uint8Array(cells);
  let prev = new Float64Array(width);
  let cur = new Float64Array(width);
  const candidates: number[] = [];
  const candidateScores: number[] = [];
  let maximum = 0;

  for (let j = 1; j <= m; j++) {
    prev[j] = local ? 0 : j * gap;
    trace[j] = local ? LIN_START : LIN_LEFT;
  }
  for (let i = 1; i <= n; i++) {
    const rowBase = i * width;
    const scoreRow = t[i - 1]! * size;
    cur[0] = local ? 0 : i * gap;
    trace[rowBase] = local ? LIN_START : LIN_UP;
    for (let j = 1; j <= m; j++) {
      let score = prev[j - 1]! + scores[scoreRow + q[j - 1]!]!;
      let bits = LIN_DIAG;
      const left = cur[j - 1]! + gap;
      if (left > score) {
        score = left;
        bits = LIN_LEFT;
      } else if (left === score) {
        bits |= LIN_LEFT;
      }
      const up = prev[j]! + gap;
      if (up > score) {
        score = up;
        bits = LIN_UP;
      } else if (up === score) {
        bits |= LIN_UP;
      }
      if (local) {
        if (score <= 0) {
          score = 0;
          bits = LIN_START;
        } else if (bits & LIN_DIAG && score >= maximum) {
          if (score > maximum) maximum = score;
          candidates.push(rowBase + j);
          candidateScores.push(score);
        }
      }
      cur[j] = score;
      trace[rowBase + j] = bits;
    }
    const swap = prev;
    prev = cur;
    cur = swap;
  }

  const ops: string[] = [];
  if (!local) {
    let i = n;
    let j = m;
    while (i > 0 || j > 0) {
      const bits = trace[i * width + j]!;
      if (bits & LIN_LEFT) {
        ops.push("I");
        j--;
      } else if (bits & LIN_UP) {
        ops.push("D");
        i--;
      } else {
        ops.push("M");
        i--;
        j--;
      }
    }
    return { score: prev[m]!, ops: ops.reverse().join(""), startI: 0, startJ: 0 };
  }

  if (maximum <= 0) return { score: 0, ops: "", startI: 0, startJ: 0 };
  for (let k = 0; k < candidates.length; k++) {
    if (candidateScores[k] === maximum) trace[candidates[k]!]! |= LIN_END;
  }
  const reachable = new Uint8Array(cells);
  reachable.fill(1, 0, width);
  for (let i = 1; i <= n; i++) {
    const rowBase = i * width;
    reachable[rowBase] = 1;
    for (let j = 1; j <= m; j++) {
      const k = rowBase + j;
      let bits = trace[k]!;
      if (!reachable[k - width - 1]) bits &= ~LIN_DIAG;
      if (!reachable[k - 1]) bits &= ~LIN_LEFT;
      if (!reachable[k - width]) bits &= ~LIN_UP;
      if (bits & (LIN_START | LIN_LEFT | LIN_UP | LIN_DIAG)) {
        reachable[k] = bits & LIN_END ? 0 : 1;
      } else {
        reachable[k] = 0;
        bits = 0;
      }
      trace[k] = bits;
    }
  }
  let end = -1;
  for (let k = width; k < cells; k++) {
    if (trace[k]! & LIN_END && trace[k]! & LIN_DIAG) {
      end = k;
      break;
    }
  }
  if (end < 0) throw new Error("no valid local alignment endpoint");
  let i = Math.floor(end / width);
  let j = end - i * width;
  ops.push("M");
  i--;
  j--;
  for (;;) {
    const bits = trace[i * width + j]!;
    if (bits & LIN_LEFT) {
      ops.push("I");
      j--;
    } else if (bits & LIN_UP) {
      ops.push("D");
      i--;
    } else if (bits & LIN_DIAG) {
      ops.push("M");
      i--;
      j--;
    } else {
      break;
    }
  }
  return { score: maximum, ops: ops.reverse().join(""), startI: i, startJ: j };
}

export function align(target: string, query: string, options: AlignOptions): Alignment {
  validateScoring(options);
  const t = encodeForMatrix(target, options.matrix, options.substitute);
  const q = encodeForMatrix(query, options.matrix, options.substitute);
  const local = options.mode === "local";
  if (options.openGap === options.extendGap) {
    const { score, ops, startI, startJ } = linearTraceback(t, q, options);
    const body = summarize(target, query, ops, startI, startJ, options.matrix);
    return { mode: options.mode, score, ...body };
  }
  const outcome = runDp(t, q, options, true);
  const { ops, startI, startJ } = traceback(outcome, q.length + 1, local);
  const body = summarize(target, query, ops, startI, startJ, options.matrix);
  return { mode: options.mode, score: outcome.score, ...body };
}

export function alignScore(target: string, query: string, options: AlignOptions): ScoreResult {
  validateScoring(options);
  const t = encodeForMatrix(target, options.matrix, options.substitute);
  const q = encodeForMatrix(query, options.matrix, options.substitute);
  const outcome = runDp(t, q, options, false);
  return { mode: options.mode, score: outcome.score, targetEnd: outcome.endI, queryEnd: outcome.endJ };
}

export function alignEncoded(t: Uint8Array, q: Uint8Array, options: AlignOptions): ScoreResult {
  validateScoring(options);
  const outcome = runDp(t, q, options, false);
  return { mode: options.mode, score: outcome.score, targetEnd: outcome.endI, queryEnd: outcome.endJ };
}

export function globalAlign(target: string, query: string, scoring: Scoring): Alignment {
  return align(target, query, { ...scoring, mode: "global" });
}

export function localAlign(target: string, query: string, scoring: Scoring): Alignment {
  return align(target, query, { ...scoring, mode: "local" });
}

export function scoreAlignedPair(alignedTarget: string, alignedQuery: string, scoring: Scoring): number {
  if (alignedTarget.length !== alignedQuery.length) {
    throw new AlignmentInputError("aligned sequences must have the same length");
  }
  const matrix = scoring.matrix;
  let total = 0;
  let previous = "M";
  for (let k = 0; k < alignedTarget.length; k++) {
    const a = alignedTarget[k]!;
    const b = alignedQuery[k]!;
    if (a === "-" && b === "-") throw new AlignmentInputError(`column ${k + 1} is gap against gap`);
    if (a === "-") {
      total += previous === "I" ? scoring.extendGap : scoring.openGap;
      previous = "I";
    } else if (b === "-") {
      total += previous === "D" ? scoring.extendGap : scoring.openGap;
      previous = "D";
    } else {
      const ia = encodeForMatrix(a, matrix, scoring.substitute)[0]!;
      const ib = encodeForMatrix(b, matrix, scoring.substitute)[0]!;
      total += matrix.scores[ia * matrix.size + ib]!;
      previous = "M";
    }
  }
  return total;
}

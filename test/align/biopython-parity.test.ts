import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { describe, it } from "node:test";
import { getMatrix, matchMismatchMatrix, type SubstitutionMatrix } from "../../src/align/matrix.ts";
import { type AlignMode, align, alignScore, scoreAlignedPair } from "../../src/align/pairwise.ts";

interface OracleCase {
  id: string;
  mode: AlignMode;
  matrix?: string;
  match?: number;
  mismatch?: number;
  open: number;
  extend: number;
  target: string;
  query: string;
  score: number;
  count: number | null;
  complete: boolean;
  alignments: [number[], number[]][];
}

const fixture = JSON.parse(
  readFileSync(path.join(import.meta.dirname, "..", "fixtures", "biopython", "align-cases.json"), "utf8"),
) as { biopython: string; cases: OracleCase[] };

function matrixFor(c: OracleCase): SubstitutionMatrix {
  if (c.matrix) return getMatrix(c.matrix);
  return matchMismatchMatrix("ACGT", c.match!, c.mismatch!);
}

const key = (coords: [number[], number[]]) => JSON.stringify(coords);

describe(`Biopython ${fixture.biopython} parity (${fixture.cases.length} cases)`, () => {
  it("scores match exactly", () => {
    const failures: string[] = [];
    for (const c of fixture.cases) {
      const options = { mode: c.mode, matrix: matrixFor(c), openGap: c.open, extendGap: c.extend };
      const full = align(c.target, c.query, options);
      const scoreOnly = alignScore(c.target, c.query, options);
      if (full.score !== c.score || scoreOnly.score !== c.score) {
        failures.push(`${c.id}: expected ${c.score}, full ${full.score}, score-only ${scoreOnly.score}`);
      }
    }
    assert.deepEqual(failures, []);
  });

  it("returned alignment re-scores to the optimal score", () => {
    const failures: string[] = [];
    for (const c of fixture.cases) {
      const options = { mode: c.mode, matrix: matrixFor(c), openGap: c.open, extendGap: c.extend };
      const result = align(c.target, c.query, options);
      const rescored = scoreAlignedPair(result.alignedTarget, result.alignedQuery, options);
      if (rescored !== c.score) failures.push(`${c.id}: rescored ${rescored}, expected ${c.score}`);
    }
    assert.deepEqual(failures, []);
  });

  it("returned alignment is one of Biopython's optimal alignments", () => {
    const failures: string[] = [];
    let checked = 0;
    for (const c of fixture.cases) {
      if (!c.complete) continue;
      checked++;
      const options = { mode: c.mode, matrix: matrixFor(c), openGap: c.open, extendGap: c.extend };
      const result = align(c.target, c.query, options);
      const allowed = new Set(c.alignments.map(key));
      if (!allowed.has(key(result.coordinates))) {
        failures.push(`${c.id}: ${key(result.coordinates)} not among ${c.alignments.length} optimal alignments`);
      }
    }
    assert.ok(checked > fixture.cases.length * 0.8, `only ${checked} cases had complete optimal sets`);
    assert.deepEqual(failures.slice(0, 10), []);
  });

  it("returned alignment matches Biopython's first alignment", () => {
    const failures: string[] = [];
    for (const c of fixture.cases) {
      const options = { mode: c.mode, matrix: matrixFor(c), openGap: c.open, extendGap: c.extend };
      const result = align(c.target, c.query, options);
      if (key(result.coordinates) !== key(c.alignments[0]!)) {
        failures.push(`${c.id}: ours ${key(result.coordinates)} first ${key(c.alignments[0]!)}`);
      }
    }
    assert.deepEqual(failures.slice(0, 10), [], `${failures.length} cases differ from Biopython's first alignment`);
  });
});

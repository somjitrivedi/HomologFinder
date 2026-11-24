import assert from "node:assert/strict";
import { describe, it } from "node:test";
import { formatAlignment } from "../../src/align/format.ts";
import { getMatrix, matchMismatchMatrix } from "../../src/align/matrix.ts";
import {
  AlignmentInputError,
  type AlignOptions,
  align,
  alignEncoded,
  alignScore,
  encodeForMatrix,
  globalAlign,
  localAlign,
  scoreAlignedPair,
} from "../../src/align/pairwise.ts";
import { STANDARD_AMINO_ACIDS, mulberry32, randomSequence } from "../../src/util/random.ts";

const blosum = { matrix: getMatrix("BLOSUM62"), openGap: -12, extendGap: -1 };
const dna = { matrix: matchMismatchMatrix("ACGT", 1, -1), openGap: -2, extendGap: -1 };

describe("global alignment (Needleman-Wunsch with affine gaps)", () => {
  it("aligns identical sequences without gaps", () => {
    const a = globalAlign("HEAGAWGHEE", "HEAGAWGHEE", blosum);
    assert.equal(a.cigar, "10M");
    assert.equal(a.stats.identities, 10);
    assert.equal(a.stats.identity, 1);
    assert.equal(a.score, scoreAlignedPair("HEAGAWGHEE", "HEAGAWGHEE", blosum));
  });

  it("charges open for the first gap position and extend for the rest", () => {
    const a = globalAlign("AAAAGGGGAAAA", "AAAAAAAA", dna);
    assert.equal(a.score, 8 + -2 + -1 * 3);
    assert.equal(a.cigar, "4M4D4M");
    assert.equal(a.stats.gaps, 4);
    assert.equal(a.stats.gapOpens, 1);
  });

  it("handles empty sequences", () => {
    const both = globalAlign("", "", blosum);
    assert.equal(both.score, 0);
    assert.equal(both.cigar, "");
    assert.deepEqual(both.coordinates, [
      [0, 0],
      [0, 0],
    ]);
    const left = globalAlign("", "ACD", blosum);
    assert.equal(left.score, -12 - 1 - 1);
    assert.equal(left.alignedTarget, "---");
    const right = globalAlign("ACD", "", blosum);
    assert.equal(right.cigar, "3D");
  });

  it("uses Biopython-style coordinates", () => {
    const a = globalAlign("AAAAGGGGAAAA", "AAAAAAAA", dna);
    assert.deepEqual(a.coordinates, [
      [0, 4, 8, 12],
      [0, 4, 4, 8],
    ]);
  });
});

describe("local alignment (Smith-Waterman with affine gaps)", () => {
  it("finds the embedded match", () => {
    const a = localAlign("TTTTTHEAGAWGHEETTTTT", "PPPHEAGAWGHEEPPP", blosum);
    assert.equal(a.alignedTarget, "HEAGAWGHEE");
    assert.equal(a.targetStart, 5);
    assert.equal(a.targetEnd, 15);
    assert.equal(a.queryStart, 3);
    assert.equal(a.queryEnd, 13);
  });

  it("returns an empty alignment with score 0 when nothing scores positively", () => {
    const a = localAlign("AAAA", "TTTT", dna);
    assert.equal(a.score, 0);
    assert.equal(a.cigar, "");
    assert.equal(a.stats.length, 0);
  });

  it("never starts or ends with a gap", () => {
    const random = mulberry32(5);
    for (let k = 0; k < 200; k++) {
      const a = localAlign(
        randomSequence(random, STANDARD_AMINO_ACIDS, 40),
        randomSequence(random, STANDARD_AMINO_ACIDS, 40),
        { ...blosum, openGap: -4, extendGap: -1 },
      );
      if (a.score === 0) continue;
      assert.equal(a.operations[0], "M");
      assert.equal(a.operations.at(-1), "M");
    }
  });
});

describe("consistency", () => {
  const random = mulberry32(17);
  const configs: Omit<AlignOptions, "mode">[] = [
    blosum,
    { matrix: getMatrix("PAM250"), openGap: -10, extendGap: -0.5 },
    { matrix: getMatrix("BLOSUM45"), openGap: -3, extendGap: -3 },
  ];

  it("score-only, encoded and traceback paths agree and alignments re-score exactly", () => {
    for (const config of configs) {
      for (const mode of ["global", "local"] as const) {
        for (let k = 0; k < 40; k++) {
          const t = randomSequence(random, STANDARD_AMINO_ACIDS, 1 + Math.floor(random() * 80));
          const q = randomSequence(random, STANDARD_AMINO_ACIDS, 1 + Math.floor(random() * 80));
          const options = { ...config, mode };
          const full = align(t, q, options);
          assert.equal(alignScore(t, q, options).score, full.score);
          const encoded = alignEncoded(encodeForMatrix(t, config.matrix), encodeForMatrix(q, config.matrix), options);
          assert.equal(encoded.score, full.score);
          assert.equal(scoreAlignedPair(full.alignedTarget, full.alignedQuery, options), full.score);
          assert.equal(full.alignedTarget.replace(/-/g, ""), t.slice(full.targetStart, full.targetEnd));
          assert.equal(full.alignedQuery.replace(/-/g, ""), q.slice(full.queryStart, full.queryEnd));
        }
      }
    }
  });

  it("local score is never below global score and never negative", () => {
    for (let k = 0; k < 100; k++) {
      const t = randomSequence(random, STANDARD_AMINO_ACIDS, 30);
      const q = randomSequence(random, STANDARD_AMINO_ACIDS, 30);
      const local = alignScore(t, q, { ...blosum, mode: "local" }).score;
      assert.ok(local >= 0);
      assert.ok(local >= alignScore(t, q, { ...blosum, mode: "global" }).score);
    }
  });

  it("is symmetric in score when the matrix is symmetric", () => {
    for (let k = 0; k < 50; k++) {
      const t = randomSequence(random, STANDARD_AMINO_ACIDS, 25);
      const q = randomSequence(random, STANDARD_AMINO_ACIDS, 35);
      for (const mode of ["global", "local"] as const) {
        assert.equal(alignScore(t, q, { ...blosum, mode }).score, alignScore(q, t, { ...blosum, mode }).score);
      }
    }
  });
});

describe("input handling", () => {
  it("rejects symbols missing from the matrix unless a substitute is given", () => {
    assert.throws(() => globalAlign("MAGU", "MAG", blosum), (error: unknown) => {
      assert.ok(error instanceof AlignmentInputError);
      assert.match(error.message, /'U' at position 4/);
      return true;
    });
    const a = globalAlign("MAGU", "MAGX", { ...blosum, substitute: "X" });
    assert.equal(a.alignedTarget, "MAGU");
    assert.throws(() => globalAlign("A", "A", { ...blosum, substitute: "U" }), /substitute/);
  });

  it("is case-insensitive", () => {
    assert.equal(globalAlign("heagawghee", "PAWHEAE", blosum).score, globalAlign("HEAGAWGHEE", "PAWHEAE", blosum).score);
  });

  it("rejects non-finite gap scores", () => {
    assert.throws(() => globalAlign("A", "A", { ...blosum, openGap: Number.NaN }), /finite/);
  });

  it("enforces the traceback cell limit but allows score-only", () => {
    const t = "A".repeat(300);
    const options = { ...blosum, mode: "global" as const, maxCells: 10_000 };
    assert.throws(() => align(t, t, options), RangeError);
    assert.throws(() => align(t, t, { ...options, openGap: -1, extendGap: -1 }), RangeError);
    assert.equal(alignScore(t, t, options).score, 300 * 4);
  });

  it("rejects malformed aligned pairs", () => {
    assert.throws(() => scoreAlignedPair("A-", "A", blosum), /same length/);
    assert.throws(() => scoreAlignedPair("A-", "A-", blosum), /gap against gap/);
  });
});

describe("formatAlignment", () => {
  it("prints blocks with coordinates and a midline", () => {
    const a = localAlign("TTTTTHEAGAWGHEETTTTT", "PPPHEAGAWGHEEPPP", blosum);
    const text = formatAlignment(a, { width: 6, targetName: "t", queryName: "q" });
    assert.match(text, /score 62/);
    assert.match(text, /t\s+6 HEAGAW 11/);
    assert.match(text, /\|\|\|\|\|\|/);
    assert.match(text, /q\s+10 GHEE 13/);
  });
});

import assert from "node:assert/strict";
import { writeFileSync } from "node:fs";
import path from "node:path";
import { describe, it } from "node:test";
import {
  SubstitutionMatrix,
  builtinMatrixNames,
  getMatrix,
  loadMatrix,
  matchMismatchMatrix,
  parseNcbiMatrix,
} from "../../src/align/matrix.ts";
import { tempDir } from "../helpers.ts";

describe("substitution matrices", () => {
  it("ships BLOSUM, PAM and NUC.4.4", () => {
    assert.deepEqual(builtinMatrixNames().sort(), [
      "BLOSUM45",
      "BLOSUM50",
      "BLOSUM62",
      "BLOSUM80",
      "BLOSUM90",
      "NUC.4.4",
      "PAM250",
      "PAM30",
      "PAM70",
    ]);
  });

  it("has the well-known BLOSUM62 values", () => {
    const m = getMatrix("blosum62");
    assert.equal(m.score("W", "W"), 11);
    assert.equal(m.score("A", "A"), 4);
    assert.equal(m.score("A", "R"), -1);
    assert.equal(m.score("w", "c"), -2);
    assert.equal(m.score("*", "*"), 1);
    assert.equal(m.alphabet, "ARNDCQEGHILKMFPSTWYVBZX*");
  });

  it("has the well-known PAM250 and NUC.4.4 values", () => {
    assert.equal(getMatrix("PAM250").score("W", "W"), 17);
    assert.equal(getMatrix("NUC.4.4").score("A", "A"), 5);
    assert.equal(getMatrix("NUC.4.4").score("A", "C"), -4);
  });

  it("every built-in matrix is symmetric and integral", () => {
    for (const name of builtinMatrixNames()) {
      const m = getMatrix(name);
      assert.ok(m.isSymmetric(), name);
      assert.ok(m.integral, name);
    }
  });

  it("caches matrices", () => {
    assert.equal(getMatrix("PAM30"), getMatrix("pam30"));
  });

  it("rejects unknown names and symbols", () => {
    assert.throws(() => getMatrix("BLOSUM99"), /unknown matrix/);
    assert.throws(() => getMatrix("BLOSUM62").score("U", "A"), /not in matrix/);
    assert.equal(getMatrix("BLOSUM62").has("U"), false);
  });

  it("builds match/mismatch matrices", () => {
    const m = matchMismatchMatrix("ACGT", 2, -3);
    assert.equal(m.score("A", "A"), 2);
    assert.equal(m.score("a", "g"), -3);
  });

  it("parses NCBI text matrices and JSON files", () => {
    const text = "# demo\n   A  B\nA  1 -1\nB -1  2\n";
    const m = parseNcbiMatrix(text, "demo");
    assert.equal(m.score("B", "B"), 2);
    const dir = tempDir();
    const txt = path.join(dir, "demo.txt");
    const json = path.join(dir, "demo.json");
    writeFileSync(txt, text);
    writeFileSync(json, JSON.stringify(m.toData()));
    assert.equal(loadMatrix(txt).score("A", "B"), -1);
    assert.equal(loadMatrix(json).score("A", "A"), 1);
    assert.equal(loadMatrix("BLOSUM62"), getMatrix("BLOSUM62"));
  });

  it("validates shape and symbols", () => {
    assert.throws(() => new SubstitutionMatrix({ name: "bad", alphabet: "AB", scores: [[1]] }), /2x2/);
    assert.throws(() => new SubstitutionMatrix({ name: "dup", alphabet: "AA", scores: [[1, 1], [1, 1]] }), /duplicate/);
    assert.throws(() => parseNcbiMatrix("  A B\nB 1 2\nA 1 2\n", "swapped"), /differ/);
  });
});

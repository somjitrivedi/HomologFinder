import assert from "node:assert/strict";
import { describe, it } from "node:test";
import {
  DNA,
  INVALID,
  PROTEIN,
  SKIP,
  alphabetById,
  alphabetByKind,
  decodeToString,
  detectAlphabet,
} from "../../src/alphabet/alphabet.ts";

describe("alphabets", () => {
  it("protein codes follow NCBI matrix order for the first 24 symbols", () => {
    assert.equal(PROTEIN.symbols.slice(0, 24), "ARNDCQEGHILKMFPSTWYVBZX*");
    for (let code = 0; code < PROTEIN.symbols.length; code++) {
      assert.equal(PROTEIN.encode[PROTEIN.symbols.charCodeAt(code)], code);
    }
  });

  it("covers every letter A-Z for protein", () => {
    for (let c = 65; c <= 90; c++) assert.ok(PROTEIN.encode[c]! >= 0, String.fromCharCode(c));
  });

  it("is case-insensitive and skips whitespace", () => {
    assert.equal(PROTEIN.encode["m".charCodeAt(0)], PROTEIN.encode["M".charCodeAt(0)]);
    assert.equal(DNA.encode["a".charCodeAt(0)], 0);
    for (const ws of [" ", "\t", "\r", "\n"]) assert.equal(DNA.encode[ws.charCodeAt(0)], SKIP);
    assert.equal(PROTEIN.encode["1".charCodeAt(0)], INVALID);
    assert.equal(DNA.encode["E".charCodeAt(0)], INVALID);
  });

  it("DNA puts ACGT first so they fit in 2 bits", () => {
    assert.equal(DNA.symbols.slice(0, 4), "ACGT");
  });

  it("round-trips codes to text", () => {
    const codes = Uint8Array.from("MKV*".split("").map((c) => PROTEIN.encode[c.charCodeAt(0)]!));
    assert.equal(decodeToString(PROTEIN, codes), "MKV*");
  });

  it("looks alphabets up by id and kind", () => {
    assert.equal(alphabetById(1), PROTEIN);
    assert.equal(alphabetById(2), DNA);
    assert.equal(alphabetByKind("dna"), DNA);
    assert.throws(() => alphabetById(9));
  });
});

describe("detectAlphabet", () => {
  const detect = (text: string) => detectAlphabet(Buffer.from(text)).kind;

  it("detects nucleotides", () => {
    assert.equal(detect(">x\nACGTACGTNNACGT\n>y\nacgtu\n"), "dna");
  });

  it("detects protein", () => {
    assert.equal(detect(">sp|P1|X\nMKTAYIAKQRQISFVKSHFSRQ\n"), "protein");
  });

  it("ignores header text", () => {
    assert.equal(detect(">PROTEIN EFLQW header words\nACGTACGT\n"), "dna");
  });

  it("treats a protein made of ACGT-like letters with others as protein", () => {
    assert.equal(detect(">x\nACGTACGTACGTE\n"), "protein");
  });

  it("defaults to protein for empty input", () => {
    assert.equal(detect(""), "protein");
  });
});

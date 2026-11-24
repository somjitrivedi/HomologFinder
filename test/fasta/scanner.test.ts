import assert from "node:assert/strict";
import { readFileSync, writeFileSync } from "node:fs";
import path from "node:path";
import { describe, it } from "node:test";
import { gzipSync } from "node:zlib";
import { parseFasta, parseFastaText, readFastaFile, formatFasta } from "../../src/fasta/records.ts";
import { FastaFormatError, FastaScanner, headerId, splitHeader } from "../../src/fasta/scanner.ts";
import { FIXTURES, chunked, tempDir, type ExpectedRecord } from "../helpers.ts";

const FASTA_DIR = path.join(FIXTURES, "fasta");
const SAMPLES = ["mixed", "swissprot-sample", "dna-sample"];

function expected(name: string): ExpectedRecord[] {
  return JSON.parse(readFileSync(path.join(FASTA_DIR, `${name}.expected.json`), "utf8")) as ExpectedRecord[];
}

async function collect(source: AsyncIterable<Uint8Array>) {
  const out = [];
  for await (const record of parseFasta(source)) out.push(record);
  return out;
}

describe("FASTA parsing matches Biopython SeqIO", () => {
  for (const name of SAMPLES) {
    it(`${name}.fasta`, async () => {
      const records = [];
      for await (const record of readFastaFile(path.join(FASTA_DIR, `${name}.fasta`))) {
        records.push({ id: record.id, description: record.description, sequence: record.sequence });
      }
      assert.deepEqual(records, expected(name));
    });
  }

  it("gives identical records for every chunk size, including 1 byte", async () => {
    const data = readFileSync(path.join(FASTA_DIR, "mixed.fasta"));
    const reference = parseFastaText(data);
    for (const size of [1, 2, 3, 7, 64, 1000]) {
      assert.deepEqual(await collect(chunked(data, size)), reference, `chunk size ${size}`);
    }
  });

  it("reads gzip-compressed files transparently", async () => {
    const dir = tempDir();
    const gz = path.join(dir, "sample.fasta.gz");
    writeFileSync(gz, gzipSync(readFileSync(path.join(FASTA_DIR, "swissprot-sample.fasta"))));
    let count = 0;
    for await (const record of readFastaFile(gz)) {
      assert.equal(record.sequence, expected("swissprot-sample")[count]!.sequence);
      count++;
    }
    assert.equal(count, 220);
  });
});

describe("FastaScanner edge cases", () => {
  it("handles CRLF line endings", () => {
    const records = parseFastaText(">a first\r\nACGT\r\nAC\r\n>b\r\nGG\r\n");
    assert.deepEqual(
      records.map((r) => [r.id, r.description, r.sequence]),
      [
        ["a", "first", "ACGTAC"],
        ["b", "", "GG"],
      ],
    );
  });

  it("keeps empty records and a final record without newline", () => {
    const records = parseFastaText(">empty\n>x\nMK");
    assert.deepEqual(
      records.map((r) => [r.id, r.sequence]),
      [
        ["empty", ""],
        ["x", "MK"],
      ],
    );
  });

  it("skips blank lines and leading whitespace before the first header", () => {
    const records = parseFastaText("\n  \n>a\n\nAC\n\nGT\n");
    assert.deepEqual(records.map((r) => r.sequence), ["ACGT"]);
  });

  it("ignores ';' comment lines", () => {
    const records = parseFastaText(";comment\n>a\n;inner\nAC\n");
    assert.deepEqual(records.map((r) => r.sequence), ["AC"]);
  });

  it("rejects sequence data before the first header with a line number", () => {
    assert.throws(() => parseFastaText("\nACGT\n>a\nAC\n"), (error: unknown) => {
      assert.ok(error instanceof FastaFormatError);
      assert.equal(error.line, 2);
      return true;
    });
  });

  it("handles a header split across chunks", () => {
    const headers: string[] = [];
    const scanner = new FastaScanner({
      header: (line) => headers.push(Buffer.from(line).toString()),
      residues: () => {},
      recordEnd: () => {},
    });
    for (const piece of [">sp|P1|", "NAME desc", "ription\r", "\nAC\n>tail"]) scanner.push(Buffer.from(piece));
    scanner.end();
    assert.deepEqual(headers, ["sp|P1|NAME description", "tail"]);
  });

  it("counts lines across chunks", () => {
    const scanner = new FastaScanner({ header: () => {}, residues: () => {}, recordEnd: () => {} });
    scanner.push(Buffer.from(">a\nAC"));
    scanner.push(Buffer.from("GT\nAC\n"));
    assert.equal(scanner.line, 4);
  });

  it("refuses input after end()", () => {
    const scanner = new FastaScanner({ header: () => {}, residues: () => {}, recordEnd: () => {} });
    scanner.end();
    assert.throws(() => scanner.push(Buffer.from(">a")));
  });
});

describe("header helpers", () => {
  it("splits id and description", () => {
    assert.deepEqual(splitHeader("sp|P12345|ABC_HUMAN  Some protein  "), {
      id: "sp|P12345|ABC_HUMAN",
      description: "Some protein  ",
    });
    assert.deepEqual(splitHeader(""), { id: "", description: "" });
  });

  it("extracts the id bytes", () => {
    assert.equal(Buffer.from(headerId(Buffer.from("  abc\tdef"))).toString(), "abc");
  });

  it("formats FASTA with wrapping", () => {
    assert.equal(formatFasta("x y", "ABCDEFG", 3), ">x y\nABC\nDEF\nG\n");
    assert.equal(formatFasta("x", "ABC", 0), ">x\nABC\n");
  });
});

import assert from "node:assert/strict";
import { existsSync, readFileSync, statSync, writeFileSync } from "node:fs";
import path from "node:path";
import { describe, it } from "node:test";
import { DNA, PROTEIN } from "../../src/alphabet/alphabet.ts";
import { parseFastaText } from "../../src/fasta/records.ts";
import { StoreFormatError } from "../../src/store/format.ts";
import { packFasta, packFastaSource } from "../../src/store/pack.ts";
import { SequenceStore } from "../../src/store/store.ts";
import { InvalidResidueError, StoreWriter } from "../../src/store/writer.ts";
import { mulberry32 } from "../../src/util/random.ts";
import { FIXTURES, chunked, tempDir, type ExpectedRecord } from "../helpers.ts";

const FASTA_DIR = path.join(FIXTURES, "fasta");

function expected(name: string): ExpectedRecord[] {
  return JSON.parse(readFileSync(path.join(FASTA_DIR, `${name}.expected.json`), "utf8")) as ExpectedRecord[];
}

function assertStoreMatches(store: SequenceStore, records: ExpectedRecord[]): void {
  assert.equal(store.count, records.length);
  records.forEach((record, i) => {
    assert.equal(store.id(i), record.id, `id #${i}`);
    assert.equal(store.description(i), record.description, `description #${i}`);
    assert.equal(store.sequence(i), record.sequence.toUpperCase(), `sequence #${i}`);
    assert.equal(store.length(i), record.sequence.length);
  });
}

describe("protein store", () => {
  const dir = tempDir();
  const db = path.join(dir, "sprot.hfdb");

  it("packs the Swiss-Prot sample and auto-detects protein", async () => {
    const result = await packFasta(path.join(FASTA_DIR, "swissprot-sample.fasta"), db);
    assert.equal(result.alphabet, "protein");
    assert.equal(result.sequenceCount, 220);
    const records = expected("swissprot-sample");
    assert.equal(result.totalResidues, records.reduce((n, r) => n + r.sequence.length, 0));
    assert.equal(result.maxLength, Math.max(...records.map((r) => r.sequence.length)));
    assert.equal(statSync(db).size, result.fileSize);
    assert.ok(!existsSync(`${db}.tmp`) && !existsSync(`${db}.headers.tmp`));
  });

  for (const residues of ["memory", "disk"] as const) {
    for (const headers of ["memory", "disk"] as const) {
      it(`reads back every record (residues ${residues}, headers ${headers})`, () => {
        const store = SequenceStore.open(db, { residues, headers });
        try {
          assertStoreMatches(store, expected("swissprot-sample"));
        } finally {
          store.close();
        }
      });
    }
  }

  it("uses one byte per residue", () => {
    const store = SequenceStore.open(db);
    try {
      assert.equal(store.header.residuesLength, store.totalResidues);
      assert.equal(store.alphabet, PROTEIN);
    } finally {
      store.close();
    }
  });

  it("finds sequences by id", () => {
    const store = SequenceStore.open(db);
    try {
      const records = expected("swissprot-sample");
      for (const index of [0, 57, 219]) assert.deepEqual(store.findById(records[index]!.id), [index]);
      assert.deepEqual(store.findById("sp|NOPE|NOPE_HUMAN"), []);
    } finally {
      store.close();
    }
  });

  it("can back residues with a SharedArrayBuffer", () => {
    const store = SequenceStore.open(db, { shared: true });
    try {
      assert.ok(store.residueBuffer!.buffer instanceof SharedArrayBuffer);
      assert.equal(store.sequence(3), expected("swissprot-sample")[3]!.sequence);
    } finally {
      store.close();
    }
  });

  it("copies codes into a caller buffer", () => {
    const store = SequenceStore.open(db, { residues: "disk" });
    try {
      const scratch = new Uint8Array(store.header.maxLength);
      const codes = store.codes(5, scratch);
      assert.equal(codes.buffer, scratch.buffer);
      assert.equal(codes.length, store.length(5));
    } finally {
      store.close();
    }
  });

  it("rejects out-of-range indices and use after close", () => {
    const store = SequenceStore.open(db, { residues: "disk" });
    assert.throws(() => store.length(-1), RangeError);
    assert.throws(() => store.length(220), RangeError);
    assert.throws(() => store.length(1.5), RangeError);
    store.close();
    assert.throws(() => store.sequence(0), /closed/);
  });
});

describe("mixed FASTA edge cases", () => {
  it("stores empty records, lowercase, U/O/J and stop symbols", async () => {
    const dir = tempDir();
    const db = path.join(dir, "mixed.hfdb");
    await packFasta(path.join(FASTA_DIR, "mixed.fasta"), db, { alphabet: "protein" });
    const store = SequenceStore.open(db);
    try {
      assertStoreMatches(store, expected("mixed"));
      const empty = store.findById("empty_record")[0]!;
      assert.equal(store.length(empty), 0);
      assert.equal(store.sequence(empty), "");
    } finally {
      store.close();
    }
  });
});

describe("DNA store", () => {
  const dir = tempDir();
  const db = path.join(dir, "dna.hfdb");

  it("auto-detects DNA and packs 4 bases per byte", async () => {
    const result = await packFasta(path.join(FASTA_DIR, "dna-sample.fasta"), db);
    assert.equal(result.alphabet, "dna");
    const store = SequenceStore.open(db);
    try {
      assert.equal(store.alphabet, DNA);
      assert.equal(store.header.residuesLength, Math.ceil(store.totalResidues / 4));
      assert.ok(result.exceptionRuns > 0);
    } finally {
      store.close();
    }
  });

  for (const residues of ["memory", "disk"] as const) {
    it(`round-trips ambiguity codes and N runs (${residues})`, () => {
      const store = SequenceStore.open(db, { residues });
      try {
        assertStoreMatches(store, expected("dna-sample"));
      } finally {
        store.close();
      }
    });
  }

  it("merges consecutive identical ambiguity symbols into runs per record", async () => {
    const local = path.join(dir, "runs.hfdb");
    await packFastaSource(chunked(Buffer.from(">a\nACNNNNGT\nNN\n>b\nNNRRA\n"), 3), local, { alphabet: "dna" });
    const store = SequenceStore.open(local);
    try {
      assert.deepEqual(
        store.exceptions(0).map((e) => [e.position, e.length, DNA.symbols[e.code]]),
        [
          [2, 4, "N"],
          [8, 2, "N"],
        ],
      );
      assert.deepEqual(
        store.exceptions(1).map((e) => [e.position, e.length, DNA.symbols[e.code]]),
        [
          [10, 2, "N"],
          [12, 2, "R"],
        ],
      );
      assert.equal(store.sequence(0), "ACNNNNGTNN");
      assert.equal(store.sequence(1), "NNRRA");
    } finally {
      store.close();
    }
  });

  it("round-trips random sequences at every bit offset", async () => {
    const random = mulberry32(99);
    const symbols = "ACGTACGTACGTACGTN";
    const records: string[] = [];
    let text = "";
    for (let r = 0; r < 200; r++) {
      const length = Math.floor(random() * 40);
      let s = "";
      for (let i = 0; i < length; i++) s += symbols[Math.floor(random() * symbols.length)];
      records.push(s);
      text += `>r${r}\n${s}\n`;
    }
    const local = path.join(dir, "random.hfdb");
    await packFastaSource(chunked(Buffer.from(text), 5), local, { alphabet: "dna" });
    for (const residues of ["memory", "disk"] as const) {
      const store = SequenceStore.open(local, { residues });
      try {
        records.forEach((s, i) => assert.equal(store.sequence(i), s, `record ${i}`));
      } finally {
        store.close();
      }
    }
  });
});

describe("invalid residues", () => {
  it("fails with the record id by default and leaves no files behind", async () => {
    const dir = tempDir();
    const db = path.join(dir, "bad.hfdb");
    await assert.rejects(
      packFastaSource(chunked(Buffer.from(">ok\nMK\n>bad_one x\nMK1V\n"), 4), db, { alphabet: "protein" }),
      (error: unknown) => {
        assert.ok(error instanceof InvalidResidueError);
        assert.equal(error.recordId, "bad_one");
        assert.equal(error.character, "1");
        return true;
      },
    );
    assert.ok(!existsSync(db) && !existsSync(`${db}.tmp`) && !existsSync(`${db}.headers.tmp`));
  });

  it("replaces with X (protein) or N (DNA) when asked", async () => {
    const dir = tempDir();
    const protein = path.join(dir, "p.hfdb");
    const dna = path.join(dir, "d.hfdb");
    const p = await packFastaSource(chunked(Buffer.from(">a\nMK1.V\n"), 64), protein, {
      alphabet: "protein",
      onInvalid: "replace",
    });
    const d = await packFastaSource(chunked(Buffer.from(">a\nAC#GT\n"), 64), dna, { alphabet: "dna", onInvalid: "replace" });
    assert.equal(p.replacedResidues, 2);
    assert.equal(d.replacedResidues, 1);
    const ps = SequenceStore.open(protein);
    const ds = SequenceStore.open(dna);
    try {
      assert.equal(ps.sequence(0), "MKXXV");
      assert.equal(ds.sequence(0), "ACNGT");
    } finally {
      ps.close();
      ds.close();
    }
  });

  it("rejects malformed FASTA and cleans up", async () => {
    const dir = tempDir();
    const db = path.join(dir, "malformed.hfdb");
    await assert.rejects(packFastaSource(chunked(Buffer.from("MKV\n>a\nMK\n"), 64), db), /before the first/);
    assert.ok(!existsSync(`${db}.tmp`));
  });
});

describe("store format validation", () => {
  it("rejects files that are not stores", () => {
    const dir = tempDir();
    const junk = path.join(dir, "junk.hfdb");
    writeFileSync(junk, Buffer.alloc(200, 7));
    assert.throws(() => SequenceStore.open(junk), StoreFormatError);
    writeFileSync(junk, Buffer.alloc(10));
    assert.throws(() => SequenceStore.open(junk), StoreFormatError);
  });

  it("rejects truncated stores", async () => {
    const dir = tempDir();
    const db = path.join(dir, "t.hfdb");
    await packFasta(path.join(FASTA_DIR, "mixed.fasta"), db);
    const bytes = readFileSync(db);
    writeFileSync(db, bytes.subarray(0, bytes.length - 20));
    assert.throws(() => SequenceStore.open(db), /truncated/);
  });

  it("writes an empty store", () => {
    const dir = tempDir();
    const db = path.join(dir, "empty.hfdb");
    const stats = new StoreWriter(db, PROTEIN).finish();
    assert.equal(stats.sequenceCount, 0);
    const store = SequenceStore.open(db);
    try {
      assert.equal(store.count, 0);
      assert.deepEqual(store.findById("x"), []);
    } finally {
      store.close();
    }
  });

  it("parses the same records that the store returns", async () => {
    const dir = tempDir();
    const db = path.join(dir, "same.hfdb");
    const text = readFileSync(path.join(FASTA_DIR, "mixed.fasta"));
    await packFastaSource(chunked(text, 1), db, { alphabet: "protein" });
    const store = SequenceStore.open(db);
    try {
      parseFastaText(text).forEach((record, i) => {
        assert.equal(store.headerLine(i), record.header);
        assert.equal(store.sequence(i), record.sequence.toUpperCase());
      });
    } finally {
      store.close();
    }
  });
});

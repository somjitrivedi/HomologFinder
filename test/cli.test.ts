import assert from "node:assert/strict";
import { spawnSync } from "node:child_process";
import path from "node:path";
import { describe, it } from "node:test";
import { FIXTURES, tempDir } from "./helpers.ts";

const CLI = path.join(import.meta.dirname, "..", "src", "cli.ts");

function cli(...args: string[]) {
  const result = spawnSync(process.execPath, [CLI, ...args], { encoding: "utf8" });
  return { code: result.status, out: result.stdout, err: result.stderr };
}

describe("cli", () => {
  const dir = tempDir();
  const db = path.join(dir, "sample.hfdb");

  it("packs, describes and fetches sequences", () => {
    const pack = cli("pack", path.join(FIXTURES, "fasta", "swissprot-sample.fasta"), "-o", db, "-q");
    assert.equal(pack.code, 0, pack.err);
    assert.match(pack.out, /alphabet\s+protein/);

    const info = cli("info", db);
    assert.equal(info.code, 0, info.err);
    assert.match(info.out, /alphabet\s+protein/);

    const get = cli("get", db, "#0", "--width", "0");
    assert.equal(get.code, 0, get.err);
    const [header, sequence] = get.out.trim().split("\n");
    assert.ok(header!.startsWith(">sp|"));
    const id = header!.slice(1).split(" ")[0]!;
    assert.equal(cli("get", db, id, "--width", "0").out, get.out);
    assert.ok(sequence!.length > 0);
  });

  it("aligns stored sequences and raw strings", () => {
    const stored = cli("align", "--db", db, "#0", "#1", "--json");
    assert.equal(stored.code, 0, stored.err);
    const parsed = JSON.parse(stored.out) as { score: number; matrix: string };
    assert.equal(parsed.matrix, "BLOSUM62");
    assert.ok(parsed.score >= 0);

    const raw = cli("align", "HEAGAWGHEE", "PAWHEAE", "--mode", "global", "--open", "-10", "--extend", "-1");
    assert.equal(raw.code, 0, raw.err);
    assert.match(raw.out, /mode global/);
    assert.match(raw.out, /matrix BLOSUM62/);

    const nucleotide = cli("align", "ACGTACGT", "ACGAACGT");
    assert.match(nucleotide.out, /matrix NUC\.4\.4/);
  });

  it("lists matrices and reports usage errors with exit code 2", () => {
    assert.match(cli("matrices").out, /BLOSUM62/);
    const bad = cli("align", "ONLYONE");
    assert.equal(bad.code, 2);
    assert.match(bad.err, /needs a target and a query/);
    assert.equal(cli("frobnicate").code, 2);
    assert.equal(cli("get", db, "no_such_id").code, 2);
  });

  it("reports runtime errors with exit code 1", () => {
    const missing = cli("info", path.join(dir, "missing.hfdb"));
    assert.equal(missing.code, 1);
    assert.match(missing.err, /ENOENT/);
  });
});

import { statSync } from "node:fs";
import { parseArgs } from "node:util";
import { formatAlignment } from "./align/format.ts";
import { builtinMatrixNames, getMatrix, loadMatrix, matchMismatchMatrix } from "./align/matrix.ts";
import { type AlignMode, align } from "./align/pairwise.ts";
import type { AlphabetKind } from "./alphabet/alphabet.ts";
import { formatFasta } from "./fasta/records.ts";
import { packFasta } from "./store/pack.ts";
import { SequenceStore } from "./store/store.ts";
import type { InvalidResiduePolicy } from "./store/writer.ts";

const USAGE = `homologfinder <command> [options]

commands
  pack <input.fasta[.gz]> -o <db.hfdb> [--alphabet auto|protein|dna] [--on-invalid error|replace]
  info <db.hfdb>
  get <db.hfdb> <id | #index>... [--width 60]
  align <target> <query> [--db db.hfdb] [--mode local|global] [--matrix BLOSUM62]
        [--open -12] [--extend -1] [--match 1 --mismatch -2] [--substitute X] [--json]
  matrices`;

class UsageError extends Error {}

function mb(bytes: number): string {
  return `${(bytes / 1e6).toFixed(1)} MB`;
}

const NUMERIC_FLAGS = new Set(["--open", "--extend", "--match", "--mismatch", "--width"]);

function joinNegativeNumbers(args: string[]): string[] {
  const out: string[] = [];
  for (let i = 0; i < args.length; i++) {
    const arg = args[i]!;
    const next = args[i + 1];
    if (NUMERIC_FLAGS.has(arg) && next !== undefined && /^-\d*\.?\d+$/.test(next)) {
      out.push(`${arg}=${next}`);
      i++;
    } else {
      out.push(arg);
    }
  }
  return out;
}

function numberOption(value: string | undefined, name: string, fallback: number): number {
  if (value === undefined) return fallback;
  const parsed = Number(value);
  if (!Number.isFinite(parsed)) throw new UsageError(`--${name} must be a number`);
  return parsed;
}

async function commandPack(args: string[]): Promise<void> {
  const { values, positionals } = parseArgs({
    args,
    allowPositionals: true,
    options: {
      output: { type: "string", short: "o" },
      alphabet: { type: "string", default: "auto" },
      "on-invalid": { type: "string", default: "error" },
      quiet: { type: "boolean", short: "q", default: false },
    },
  });
  const input = positionals[0];
  if (!input || !values.output) throw new UsageError("pack needs an input FASTA and -o <output>");
  const alphabet = values.alphabet as AlphabetKind | "auto";
  if (!["auto", "protein", "dna"].includes(alphabet)) throw new UsageError("--alphabet must be auto, protein or dna");
  const onInvalid = values["on-invalid"] as InvalidResiduePolicy;
  if (!["error", "replace"].includes(onInvalid)) throw new UsageError("--on-invalid must be error or replace");
  const inputSize = statSync(input).size;
  let lastReport = 0;
  const result = await packFasta(input, values.output, {
    alphabet,
    onInvalid,
    onProgress: values.quiet
      ? undefined
      : (bytes, sequences) => {
          const now = Date.now();
          if (now - lastReport < 500) return;
          lastReport = now;
          process.stderr.write(`\r${sequences.toLocaleString()} sequences, ${mb(bytes)} read   `);
        },
  });
  if (!values.quiet) process.stderr.write("\r".padEnd(60) + "\r");
  const detection = result.detection
    ? ` (detected, ${(result.detection.nucleotideFraction * 100).toFixed(1)}% nucleotide letters)`
    : "";
  console.log(`packed ${input} -> ${values.output}`);
  console.log(`  alphabet       ${result.alphabet}${detection}`);
  console.log(`  sequences      ${result.sequenceCount.toLocaleString()}`);
  console.log(`  residues       ${result.totalResidues.toLocaleString()} (longest ${result.maxLength.toLocaleString()})`);
  if (result.replacedResidues) console.log(`  replaced       ${result.replacedResidues.toLocaleString()} invalid residues`);
  if (result.exceptionRuns) console.log(`  ambiguity runs ${result.exceptionRuns.toLocaleString()}`);
  console.log(`  input          ${mb(inputSize)} on disk, ${mb(result.inputBytes)} uncompressed`);
  console.log(`  store          ${mb(result.fileSize)}`);
  console.log(`  time           ${(result.elapsedMs / 1000).toFixed(2)} s`);
}

function commandInfo(args: string[]): void {
  const path = args[0];
  if (!path) throw new UsageError("info needs a store path");
  const store = SequenceStore.open(path, { residues: "disk" });
  try {
    const h = store.header;
    console.log(`store          ${path}`);
    console.log(`format         v${h.version}`);
    console.log(`alphabet       ${store.alphabet.kind}`);
    console.log(`sequences      ${h.sequenceCount.toLocaleString()}`);
    console.log(`residues       ${h.totalResidues.toLocaleString()}`);
    console.log(`longest        ${h.maxLength.toLocaleString()}`);
    console.log(`mean length    ${h.sequenceCount ? (h.totalResidues / h.sequenceCount).toFixed(1) : 0}`);
    console.log(`residue bytes  ${mb(h.residuesLength)}`);
    console.log(`header bytes   ${mb(h.headersLength)}`);
    console.log(`ambiguity runs ${h.exceptionCount.toLocaleString()}`);
    console.log(`file size      ${mb(statSync(path).size)}`);
  } finally {
    store.close();
  }
}

function resolveIndex(store: SequenceStore, key: string): number {
  if (/^#\d+$/.test(key)) return Number(key.slice(1));
  const hits = store.findById(key);
  if (hits.length === 0) throw new UsageError(`no sequence with id '${key}'`);
  return hits[0]!;
}

function commandGet(args: string[]): void {
  const { values, positionals } = parseArgs({
    args,
    allowPositionals: true,
    options: { width: { type: "string", default: "60" } },
  });
  const [path, ...keys] = positionals;
  if (!path || keys.length === 0) throw new UsageError("get needs a store path and at least one id or #index");
  const width = numberOption(values.width, "width", 60);
  const store = SequenceStore.open(path, { residues: "disk" });
  try {
    for (const key of keys) {
      const index = resolveIndex(store, key);
      process.stdout.write(formatFasta(store.headerLine(index), store.sequence(index), width));
    }
  } finally {
    store.close();
  }
}

function commandAlign(args: string[]): void {
  const { values, positionals } = parseArgs({
    args,
    allowPositionals: true,
    options: {
      db: { type: "string" },
      mode: { type: "string", default: "local" },
      matrix: { type: "string" },
      open: { type: "string" },
      extend: { type: "string" },
      match: { type: "string" },
      mismatch: { type: "string" },
      substitute: { type: "string" },
      width: { type: "string", default: "60" },
      json: { type: "boolean", default: false },
    },
  });
  if (positionals.length !== 2) throw new UsageError("align needs a target and a query");
  const mode = values.mode as AlignMode;
  if (mode !== "local" && mode !== "global") throw new UsageError("--mode must be local or global");

  let target = positionals[0]!;
  let query = positionals[1]!;
  let targetName = "target";
  let queryName = "query";
  let kind: AlphabetKind = /^[ACGTUN]+$/i.test(target + query) ? "dna" : "protein";
  if (values.db) {
    const store = SequenceStore.open(values.db, { residues: "disk" });
    try {
      const ti = resolveIndex(store, target);
      const qi = resolveIndex(store, query);
      targetName = store.id(ti);
      queryName = store.id(qi);
      target = store.sequence(ti);
      query = store.sequence(qi);
      kind = store.alphabet.kind;
    } finally {
      store.close();
    }
  }

  const useMatchMismatch = values.match !== undefined || values.mismatch !== undefined;
  const matrix = useMatchMismatch
    ? matchMismatchMatrix("ACGTUN", numberOption(values.match, "match", 1), numberOption(values.mismatch, "mismatch", -2))
    : values.matrix
      ? loadMatrix(values.matrix)
      : getMatrix(kind === "dna" ? "NUC.4.4" : "BLOSUM62");
  const dnaDefaults = kind === "dna";
  const openGap = numberOption(values.open, "open", dnaDefaults ? -5 : -12);
  const extendGap = numberOption(values.extend, "extend", dnaDefaults ? -2 : -1);
  const substitute = values.substitute ?? (matrix.has("X") && kind === "protein" ? "X" : matrix.has("N") ? "N" : undefined);

  const result = align(target, query, { mode, matrix, openGap, extendGap, substitute });
  if (values.json) {
    console.log(JSON.stringify({ matrix: matrix.name, openGap, extendGap, ...result }, null, 2));
    return;
  }
  console.log(`matrix ${matrix.name}  open ${openGap}  extend ${extendGap}`);
  console.log(formatAlignment(result, { width: numberOption(values.width, "width", 60), targetName, queryName }));
}

async function main(argv: string[]): Promise<number> {
  const [command, ...raw] = argv;
  const rest = joinNegativeNumbers(raw);
  try {
    switch (command) {
      case "pack":
        await commandPack(rest);
        return 0;
      case "info":
        commandInfo(rest);
        return 0;
      case "get":
        commandGet(rest);
        return 0;
      case "align":
        commandAlign(rest);
        return 0;
      case "matrices":
        console.log(builtinMatrixNames().join("\n"));
        return 0;
      case undefined:
      case "help":
      case "--help":
      case "-h":
        console.log(USAGE);
        return command === undefined ? 1 : 0;
      default:
        throw new UsageError(`unknown command '${command}'`);
    }
  } catch (error) {
    if (error instanceof UsageError) {
      console.error(`error: ${error.message}\n\n${USAGE}`);
      return 2;
    }
    console.error(`error: ${(error as Error).message}`);
    return 1;
  }
}

process.exitCode = await main(process.argv.slice(2));

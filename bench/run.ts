import { spawnSync } from "node:child_process";
import { appendFileSync, existsSync, mkdirSync } from "node:fs";
import { cpus } from "node:os";
import path from "node:path";
import type { BenchReport } from "./report.ts";

const root = path.resolve(import.meta.dirname, "..");
const resultsDir = path.join(root, "bench", "results");
const CASES: Record<string, string[][]> = {
  pack: [[]],
  open: [["", "memory"], ["", "disk"]],
  align: [[]],
};

const requested = process.argv.slice(2).filter((arg) => !arg.startsWith("--"));
const names = requested.length > 0 ? requested : Object.keys(CASES);
const unknown = names.filter((name) => !CASES[name]);
if (unknown.length > 0) {
  console.error(`unknown bench case(s): ${unknown.join(", ")}; available: ${Object.keys(CASES).join(", ")}`);
  process.exit(2);
}

const needsData = names.some((name) => name === "pack" || name === "open");
if (needsData && !existsSync(path.join(root, "data", "uniprot_sprot.fasta.gz"))) {
  console.error("data/uniprot_sprot.fasta.gz missing, run `npm run data:swissprot` first");
  process.exit(2);
}

mkdirSync(resultsDir, { recursive: true });
const environment = {
  date: new Date().toISOString(),
  node: process.version,
  platform: `${process.platform}-${process.arch}`,
  cpu: cpus()[0]?.model.trim() ?? "unknown",
  cores: cpus().length,
};

let failed = false;
for (const name of names) {
  if (name === "open" && !existsSync(path.join(root, "data", "swissprot.hfdb"))) {
    console.error("data/swissprot.hfdb missing, run the pack case first");
    failed = true;
    continue;
  }
  for (const variant of CASES[name]!) {
    const args = variant.map((arg) => arg || path.join(root, "data", "swissprot.hfdb"));
    const child = spawnSync(process.execPath, [path.join(root, "bench", "cases", `${name}.ts`), ...args], {
      encoding: "utf8",
      maxBuffer: 1 << 26,
    });
    const line = child.stdout.split("\n").find((l) => l.startsWith("BENCH_RESULT "));
    if (child.status !== 0 || !line) {
      console.error(`bench ${name} failed (exit ${child.status})\n${child.stderr}`);
      failed = true;
      continue;
    }
    const report = JSON.parse(line.slice("BENCH_RESULT ".length)) as BenchReport;
    appendFileSync(path.join(resultsDir, "history.jsonl"), `${JSON.stringify({ ...environment, ...report })}\n`);
    console.log(`\n${report.case}  (peak RSS ${report.peakRssMb.toFixed(0)} MB)`);
    for (const metric of report.metrics) {
      const value = metric.value >= 1000 ? metric.value.toFixed(0) : metric.value.toPrecision(4);
      console.log(`  ${metric.name.padEnd(22)} ${value.padStart(12)} ${metric.unit}`);
    }
  }
}
console.log(`\n${environment.cpu}, ${environment.cores} cores, Node ${environment.node}`);
process.exit(failed ? 1 : 0);

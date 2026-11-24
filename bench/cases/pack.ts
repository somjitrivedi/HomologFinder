import { statSync } from "node:fs";
import path from "node:path";
import { packFasta } from "../../src/store/pack.ts";
import { emit, peakRssMb } from "../report.ts";

const root = path.resolve(import.meta.dirname, "..", "..");
const input = process.argv[2] ?? path.join(root, "data", "uniprot_sprot.fasta.gz");
const output = process.argv[3] ?? path.join(root, "data", "swissprot.hfdb");

const result = await packFasta(input, output);
emit({
  case: "pack",
  metrics: [
    { name: "wall time", value: result.elapsedMs / 1000, unit: "s" },
    { name: "sequences", value: result.sequenceCount, unit: "seq" },
    { name: "residues", value: result.totalResidues, unit: "res" },
    { name: "throughput", value: result.inputBytes / 1e6 / (result.elapsedMs / 1000), unit: "MB/s uncompressed" },
    { name: "store size", value: result.fileSize / 1e6, unit: "MB" },
  ],
  peakRssMb: peakRssMb(),
  details: { input: path.relative(root, input), inputMb: statSync(input).size / 1e6, alphabet: result.alphabet },
});

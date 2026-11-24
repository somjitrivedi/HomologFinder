import path from "node:path";
import { type Residency, SequenceStore } from "../../src/store/store.ts";
import { mulberry32 } from "../../src/util/random.ts";
import { emit, peakRssMb, time } from "../report.ts";

const root = path.resolve(import.meta.dirname, "..", "..");
const storePath = process.argv[2] ?? path.join(root, "data", "swissprot.hfdb");
const residency = (process.argv[3] ?? "memory") as Residency;
const samples = 100_000;

const opened = time(() => SequenceStore.open(storePath, { residues: residency }));
const store = opened.result;
const random = mulberry32(7);
const indices = Array.from({ length: samples }, () => Math.floor(random() * store.count));

let residues = 0;
const access = time(() => {
  for (const index of indices) residues += store.codes(index).length;
});

const ids = indices.slice(0, 10_000).map((index) => store.id(index));
const lookup = time(() => {
  let found = 0;
  for (const id of ids) found += store.findById(id).length;
  return found;
});

const scan = time(() => {
  let checksum = 0;
  for (let i = 0; i < store.count; i++) {
    const codes = store.codes(i);
    checksum = (checksum + codes.length + (codes[0] ?? 0)) >>> 0;
  }
  return checksum;
});

emit({
  case: `open-${residency}`,
  metrics: [
    { name: "open time", value: opened.ms, unit: "ms" },
    { name: "random access", value: (access.ms * 1000) / samples, unit: "us/seq" },
    { name: "id lookup", value: (lookup.ms * 1000) / ids.length, unit: "us/id" },
    { name: "full scan", value: scan.ms, unit: "ms" },
  ],
  peakRssMb: peakRssMb(),
  details: { sequences: store.count, sampledResidues: residues, idsFound: lookup.result },
});
store.close();

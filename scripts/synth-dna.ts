import { createWriteStream } from "node:fs";
import { once } from "node:events";
import path from "node:path";
import { parseArgs } from "node:util";
import { mulberry32 } from "../src/util/random.ts";

const { values } = parseArgs({
  options: {
    count: { type: "string", default: "100000" },
    min: { type: "string", default: "200" },
    max: { type: "string", default: "4000" },
    seed: { type: "string", default: "42" },
    width: { type: "string", default: "70" },
    out: { type: "string", default: path.resolve(import.meta.dirname, "..", "data", "synth-dna.fasta") },
  },
});

const count = Number(values.count);
const min = Number(values.min);
const max = Number(values.max);
const width = Number(values.width);
const random = mulberry32(Number(values.seed));
const BASES = "ACGT";
const IUPAC = "RYKMSWBDHV";

const out = createWriteStream(values.out);
let residues = 0;
for (let r = 0; r < count; r++) {
  const length = min + Math.floor(random() * (max - min + 1));
  const chars: string[] = new Array(length);
  for (let i = 0; i < length; i++) chars[i] = BASES[Math.floor(random() * 4)]!;
  if (random() < 0.3) {
    const start = Math.floor(random() * length);
    const run = 1 + Math.floor(random() * 50);
    for (let i = start; i < Math.min(length, start + run); i++) chars[i] = "N";
  }
  if (random() < 0.2) chars[Math.floor(random() * length)] = IUPAC[Math.floor(random() * IUPAC.length)]!;
  if (random() < 0.1) {
    const start = Math.floor(random() * length);
    for (let i = start; i < Math.min(length, start + 30); i++) chars[i] = chars[i]!.toLowerCase();
  }
  const sequence = chars.join("");
  residues += length;
  let block = `>synth_${r + 1} length=${length} seed=${values.seed}\n`;
  for (let i = 0; i < length; i += width) block += `${sequence.slice(i, i + width)}\n`;
  if (!out.write(block)) await once(out, "drain");
}
out.end();
await once(out, "finish");
console.log(`wrote ${count} sequences, ${residues} bases to ${values.out}`);

import { existsSync, writeFileSync } from "node:fs";
import path from "node:path";
import { formatFasta, readFastaFile } from "../src/fasta/records.ts";

const root = path.resolve(import.meta.dirname, "..");
const input = process.argv[2] ?? path.join(root, "data", "uniprot_sprot.fasta.gz");
const output = process.argv[3] ?? path.join(root, "test", "fixtures", "fasta", "swissprot-sample.fasta");
const headCount = 200;
const unusualCount = 20;

if (!existsSync(input)) {
  console.error(`${input} missing, run \`npm run data:swissprot\` first`);
  process.exit(2);
}

const chunks: string[] = [];
let taken = 0;
let unusual = 0;
for await (const record of readFastaFile(input)) {
  if (taken < headCount) {
    chunks.push(formatFasta(record.header, record.sequence));
    taken++;
  } else if (unusual < unusualCount && /[UOBZXJ]/.test(record.sequence) && record.sequence.length < 600) {
    chunks.push(formatFasta(record.header, record.sequence));
    unusual++;
  }
  if (taken >= headCount && unusual >= unusualCount) break;
}
writeFileSync(output, chunks.join(""));
console.log(`wrote ${taken + unusual} records (${unusual} with non-standard residues) to ${output}`);

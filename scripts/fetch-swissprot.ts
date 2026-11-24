import { createWriteStream, existsSync, mkdirSync, renameSync, rmSync, statSync } from "node:fs";
import type { IncomingMessage } from "node:http";
import { get } from "node:https";
import { Transform } from "node:stream";
import { pipeline } from "node:stream/promises";
import path from "node:path";

const MIRRORS = [
  "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
  "https://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz",
];

const dataDir = path.resolve(import.meta.dirname, "..", "data");
const target = path.join(dataDir, "uniprot_sprot.fasta.gz");
const force = process.argv.includes("--force");

function request(url: string, redirects = 5): Promise<IncomingMessage> {
  return new Promise((resolve, reject) => {
    const req = get(url, { timeout: 60_000 }, (res) => {
      const status = res.statusCode ?? 0;
      if (status >= 300 && status < 400 && res.headers.location && redirects > 0) {
        res.resume();
        resolve(request(new URL(res.headers.location, url).toString(), redirects - 1));
        return;
      }
      if (status !== 200) {
        res.resume();
        reject(new Error(`HTTP ${status} from ${url}`));
        return;
      }
      resolve(res);
    });
    req.on("timeout", () => req.destroy(new Error(`timeout connecting to ${url}`)));
    req.on("error", reject);
  });
}

async function download(url: string, destination: string): Promise<number> {
  const res = await request(url);
  const total = Number(res.headers["content-length"] ?? 0);
  let received = 0;
  let lastReport = 0;
  const progress = new Transform({
    transform(chunk: Buffer, _encoding, callback) {
      received += chunk.length;
      const now = Date.now();
      if (now - lastReport > 1000) {
        lastReport = now;
        const pct = total ? ` ${((received / total) * 100).toFixed(1)}%` : "";
        process.stdout.write(`\r${(received / 1e6).toFixed(1)} MB${pct}   `);
      }
      callback(null, chunk);
    },
  });
  await pipeline(res, progress, createWriteStream(destination));
  if (total && received !== total) throw new Error(`truncated download: ${received} of ${total} bytes`);
  return received;
}

if (existsSync(target) && !force) {
  console.log(`already present: ${target} (${(statSync(target).size / 1e6).toFixed(1)} MB), use --force to refetch`);
  process.exit(0);
}

mkdirSync(dataDir, { recursive: true });
const partial = `${target}.part`;
let lastError: unknown;
for (const url of MIRRORS) {
  try {
    console.log(`fetching ${url}`);
    const bytes = await download(url, partial);
    renameSync(partial, target);
    process.stdout.write(`\rsaved ${target} (${(bytes / 1e6).toFixed(1)} MB)\n`);
    process.exit(0);
  } catch (error) {
    lastError = error;
    rmSync(partial, { force: true });
    console.error(`\n${(error as Error).message}`);
  }
}
console.error(`all mirrors failed: ${(lastError as Error)?.message}`);
process.exit(1);

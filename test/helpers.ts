import { mkdtempSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import path from "node:path";
import { after } from "node:test";

export const FIXTURES = path.join(import.meta.dirname, "fixtures");

export function tempDir(): string {
  const dir = mkdtempSync(path.join(tmpdir(), "homologfinder-"));
  after(() => rmSync(dir, { recursive: true, force: true }));
  return dir;
}

export async function* chunked(data: Uint8Array, size: number): AsyncGenerator<Uint8Array> {
  for (let i = 0; i < data.length; i += size) yield data.subarray(i, Math.min(data.length, i + size));
}

export interface ExpectedRecord {
  id: string;
  description: string;
  sequence: string;
}

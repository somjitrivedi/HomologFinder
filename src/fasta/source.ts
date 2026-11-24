import { closeSync, createReadStream, openSync, readSync } from "node:fs";
import { createGunzip } from "node:zlib";
import type { Readable } from "node:stream";

export const DEFAULT_CHUNK_SIZE = 1 << 20;

export function isGzipFile(path: string): boolean {
  const fd = openSync(path, "r");
  try {
    const magic = Buffer.alloc(2);
    const read = readSync(fd, magic, 0, 2, 0);
    return read === 2 && magic[0] === 0x1f && magic[1] === 0x8b;
  } finally {
    closeSync(fd);
  }
}

export function openFastaStream(path: string, chunkSize = DEFAULT_CHUNK_SIZE): Readable {
  const raw = createReadStream(path, { highWaterMark: chunkSize });
  if (!isGzipFile(path)) return raw;
  const gunzip = createGunzip({ chunkSize });
  raw.on("error", (error) => gunzip.destroy(error));
  return raw.pipe(gunzip);
}

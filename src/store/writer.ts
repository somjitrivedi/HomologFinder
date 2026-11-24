import { closeSync, openSync, readSync, renameSync, rmSync, writeSync } from "node:fs";
import { type Alphabet, SKIP } from "../alphabet/alphabet.ts";
import { headerId } from "../fasta/scanner.ts";
import { Float64List, Uint32List, Uint8List } from "../util/growable.ts";
import {
  EXCEPTION_ENTRY_SIZE,
  FORMAT_VERSION,
  HEADER_SIZE,
  ID_ENTRY_SIZE,
  align8,
  encodeHeader,
  fnv1a32,
  type StoreHeader,
} from "./format.ts";

export type InvalidResiduePolicy = "error" | "replace";

export class InvalidResidueError extends Error {
  readonly recordIndex: number;
  readonly recordId: string;
  readonly character: string;
  constructor(recordIndex: number, recordId: string, character: string, alphabet: string) {
    super(`invalid ${alphabet} residue '${character}' in record #${recordIndex + 1} (${recordId})`);
    this.name = "InvalidResidueError";
    this.recordIndex = recordIndex;
    this.recordId = recordId;
    this.character = character;
  }
}

export interface StoreWriterOptions {
  onInvalid?: InvalidResiduePolicy;
  bufferSize?: number;
}

export interface StoreStats {
  sequenceCount: number;
  totalResidues: number;
  maxLength: number;
  replacedResidues: number;
  exceptionRuns: number;
  fileSize: number;
}

const U64_CHUNK = 1 << 18;
const DNA_BASES = 4;

class BufferedFile {
  readonly fd: number;
  position: number;
  private readonly buffer: Uint8Array;
  private used = 0;

  constructor(fd: number, position: number, size: number) {
    this.fd = fd;
    this.position = position;
    this.buffer = new Uint8Array(size);
  }

  tell(): number {
    return this.position + this.used;
  }

  byte(value: number): void {
    if (this.used === this.buffer.length) this.flush();
    this.buffer[this.used++] = value;
  }

  bytes(data: Uint8Array): void {
    if (data.length > this.buffer.length - this.used) {
      this.flush();
      if (data.length > this.buffer.length) {
        this.writeDirect(data);
        return;
      }
    }
    this.buffer.set(data, this.used);
    this.used += data.length;
  }

  flush(): void {
    if (this.used === 0) return;
    this.writeDirect(this.buffer.subarray(0, this.used));
    this.used = 0;
  }

  private writeDirect(data: Uint8Array): void {
    let written = 0;
    while (written < data.length) {
      written += writeSync(this.fd, data, written, data.length - written, this.position + written);
    }
    this.position += data.length;
  }
}

export class StoreWriter {
  readonly path: string;
  readonly alphabet: Alphabet;
  private readonly tempPath: string;
  private readonly headersTempPath: string;
  private readonly fd: number;
  private readonly headersFd: number;
  private readonly residues: BufferedFile;
  private readonly headers: BufferedFile;
  private readonly onInvalid: InvalidResiduePolicy;
  private readonly isDna: boolean;
  private readonly replacement: number;

  private readonly sequenceOffsets = new Float64List(1 << 16);
  private readonly headerOffsets = new Float64List(1 << 16);
  private readonly exceptionIndex = new Float64List(1 << 16);
  private readonly idHashes = new Uint32List(1 << 16);
  private readonly exceptionPositions = new Float64List();
  private readonly exceptionLengths = new Uint32List();
  private readonly exceptionCodes = new Uint8List();

  private totalResidues = 0;
  private headerBytes = 0;
  private maxLength = 0;
  private recordStart = 0;
  private recordOpen = false;
  private currentId = "";
  private replaced = 0;
  private packed = 0;
  private packedSlots = 0;
  private runOpen = false;
  private closed = false;

  constructor(path: string, alphabet: Alphabet, options: StoreWriterOptions = {}) {
    this.path = path;
    this.alphabet = alphabet;
    this.tempPath = `${path}.tmp`;
    this.headersTempPath = `${path}.headers.tmp`;
    this.onInvalid = options.onInvalid ?? "error";
    this.isDna = alphabet.kind === "dna";
    this.replacement = alphabet.encode[this.isDna ? 78 : 88]!;
    const bufferSize = options.bufferSize ?? 1 << 22;
    this.fd = openSync(this.tempPath, "w");
    this.headersFd = openSync(this.headersTempPath, "w+");
    this.residues = new BufferedFile(this.fd, HEADER_SIZE, bufferSize);
    this.headers = new BufferedFile(this.headersFd, 0, 1 << 20);
  }

  get count(): number {
    return this.sequenceOffsets.length;
  }

  beginRecord(header: Uint8Array): void {
    if (this.recordOpen) this.endRecord();
    this.recordOpen = true;
    this.recordStart = this.totalResidues;
    this.sequenceOffsets.push(this.totalResidues);
    this.headerOffsets.push(this.headerBytes);
    this.exceptionIndex.push(this.exceptionPositions.length);
    const id = headerId(header);
    this.idHashes.push(fnv1a32(id));
    this.currentId = "";
    if (this.onInvalid === "error") this.currentId = Buffer.from(id).toString("utf8");
    this.headers.bytes(header);
    this.headerBytes += header.length;
    this.runOpen = false;
  }

  addResidues(chunk: Uint8Array, start: number, end: number): void {
    if (!this.recordOpen) throw new Error("addResidues called outside a record");
    if (this.isDna) this.addDna(chunk, start, end);
    else this.addProtein(chunk, start, end);
  }

  endRecord(): void {
    if (!this.recordOpen) return;
    this.recordOpen = false;
    this.runOpen = false;
    const length = this.totalResidues - this.recordStart;
    if (length > this.maxLength) this.maxLength = length;
  }

  private resolveInvalid(byte: number): number {
    if (this.onInvalid === "error") {
      throw new InvalidResidueError(this.count - 1, this.currentId, String.fromCharCode(byte), this.alphabet.kind);
    }
    this.replaced++;
    return this.replacement;
  }

  private addProtein(chunk: Uint8Array, start: number, end: number): void {
    const encode = this.alphabet.encode;
    const out = this.residues;
    let added = 0;
    for (let i = start; i < end; i++) {
      let code = encode[chunk[i]!]!;
      if (code < 0) {
        if (code === SKIP) continue;
        code = this.resolveInvalid(chunk[i]!);
      }
      out.byte(code);
      added++;
    }
    this.totalResidues += added;
  }

  private addDna(chunk: Uint8Array, start: number, end: number): void {
    const encode = this.alphabet.encode;
    for (let i = start; i < end; i++) {
      let code = encode[chunk[i]!]!;
      if (code < 0) {
        if (code === SKIP) continue;
        code = this.resolveInvalid(chunk[i]!);
      }
      let base = code;
      if (code >= DNA_BASES) {
        base = 0;
        const last = this.exceptionPositions.length - 1;
        if (
          this.runOpen &&
          this.exceptionCodes.get(last) === code &&
          this.exceptionPositions.get(last) + this.exceptionLengths.get(last) === this.totalResidues
        ) {
          this.exceptionLengths.set(last, this.exceptionLengths.get(last) + 1);
        } else {
          this.exceptionPositions.push(this.totalResidues);
          this.exceptionLengths.push(1);
          this.exceptionCodes.push(code);
          this.runOpen = true;
        }
      } else {
        this.runOpen = false;
      }
      this.packed |= base << (this.packedSlots * 2);
      this.packedSlots++;
      if (this.packedSlots === 4) {
        this.residues.byte(this.packed);
        this.packed = 0;
        this.packedSlots = 0;
      }
      this.totalResidues++;
    }
  }

  finish(): StoreStats {
    if (this.closed) throw new Error("writer already closed");
    this.endRecord();
    try {
      if (this.isDna && this.packedSlots > 0) {
        this.residues.byte(this.packed);
        this.packed = 0;
        this.packedSlots = 0;
      }
      this.residues.flush();
      this.headers.flush();
      const residuesLength = this.residues.tell() - HEADER_SIZE;
      const count = this.count;

      this.sequenceOffsets.push(this.totalResidues);
      this.headerOffsets.push(this.headerBytes);
      this.exceptionIndex.push(this.exceptionPositions.length);

      const out = new BufferedFile(this.fd, align8(this.residues.tell()), 1 << 22);
      const sequenceOffsetsOffset = out.tell();
      this.writeU64(out, this.sequenceOffsets.view());
      const headerOffsetsOffset = out.tell();
      this.writeU64(out, this.headerOffsets.view());
      const exceptionIndexOffset = out.tell();
      this.writeU64(out, this.exceptionIndex.view());
      const exceptionsOffset = out.tell();
      this.writeExceptions(out);
      const idIndexOffset = out.tell();
      this.writeIdIndex(out, count);
      const headersOffset = out.tell();
      this.copyHeaders(out);
      out.flush();

      const header: StoreHeader = {
        version: FORMAT_VERSION,
        alphabetId: this.alphabet.id,
        sequenceCount: count,
        totalResidues: this.totalResidues,
        maxLength: this.maxLength,
        residuesOffset: HEADER_SIZE,
        residuesLength,
        sequenceOffsetsOffset,
        headerOffsetsOffset,
        exceptionIndexOffset,
        exceptionsOffset,
        exceptionCount: this.exceptionPositions.length,
        idIndexOffset,
        headersOffset,
        headersLength: this.headerBytes,
      };
      writeSync(this.fd, encodeHeader(header), 0, HEADER_SIZE, 0);
      const fileSize = out.tell();
      this.closeHandles();
      rmSync(this.headersTempPath, { force: true });
      renameSync(this.tempPath, this.path);
      return {
        sequenceCount: count,
        totalResidues: this.totalResidues,
        maxLength: this.maxLength,
        replacedResidues: this.replaced,
        exceptionRuns: this.exceptionPositions.length,
        fileSize,
      };
    } catch (error) {
      this.abort();
      throw error;
    }
  }

  abort(): void {
    this.closeHandles();
    rmSync(this.tempPath, { force: true });
    rmSync(this.headersTempPath, { force: true });
  }

  private closeHandles(): void {
    if (this.closed) return;
    this.closed = true;
    closeSync(this.fd);
    closeSync(this.headersFd);
  }

  private writeU64(out: BufferedFile, values: Float64Array): void {
    const scratch = new BigUint64Array(Math.min(U64_CHUNK, Math.max(1, values.length)));
    const bytes = new Uint8Array(scratch.buffer);
    for (let start = 0; start < values.length; start += U64_CHUNK) {
      const n = Math.min(U64_CHUNK, values.length - start);
      for (let i = 0; i < n; i++) scratch[i] = BigInt(values[start + i]!);
      out.bytes(bytes.subarray(0, n * 8));
    }
  }

  private writeExceptions(out: BufferedFile): void {
    const entry = Buffer.alloc(EXCEPTION_ENTRY_SIZE);
    const n = this.exceptionPositions.length;
    for (let i = 0; i < n; i++) {
      entry.writeBigUInt64LE(BigInt(this.exceptionPositions.get(i)), 0);
      entry.writeUInt32LE(this.exceptionLengths.get(i), 8);
      entry.writeUInt8(this.exceptionCodes.get(i), 12);
      out.bytes(entry);
    }
  }

  private writeIdIndex(out: BufferedFile, count: number): void {
    const keys = new BigUint64Array(count);
    const hashes = this.idHashes.view();
    for (let i = 0; i < count; i++) keys[i] = (BigInt(hashes[i]!) << 32n) | BigInt(i);
    keys.sort();
    const pairs = new Uint32Array(count * 2);
    for (let i = 0; i < count; i++) {
      const key = keys[i]!;
      pairs[i * 2] = Number(key >> 32n);
      pairs[i * 2 + 1] = Number(key & 0xffffffffn);
    }
    out.bytes(new Uint8Array(pairs.buffer, 0, count * ID_ENTRY_SIZE));
  }

  private copyHeaders(out: BufferedFile): void {
    const chunk = Buffer.allocUnsafe(1 << 20);
    let position = 0;
    while (position < this.headerBytes) {
      const read = readSync(this.headersFd, chunk, 0, Math.min(chunk.length, this.headerBytes - position), position);
      if (read <= 0) throw new Error("failed to read back header data");
      out.bytes(chunk.subarray(0, read));
      position += read;
    }
  }
}

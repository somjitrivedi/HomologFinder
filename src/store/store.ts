import { closeSync, fstatSync, openSync, readSync } from "node:fs";
import { type Alphabet, alphabetById, decodeToString } from "../alphabet/alphabet.ts";
import { splitHeader } from "../fasta/scanner.ts";
import {
  EXCEPTION_ENTRY_SIZE,
  HEADER_SIZE,
  StoreFormatError,
  decodeHeader,
  fnv1a32,
  type StoreHeader,
} from "./format.ts";

export type Residency = "memory" | "disk";

export interface OpenOptions {
  residues?: Residency;
  headers?: Residency;
  shared?: boolean;
}

const READ_CHUNK = 1 << 26;

function readExact(fd: number, target: Uint8Array, position: number): void {
  let done = 0;
  while (done < target.length) {
    const n = readSync(fd, target, done, Math.min(READ_CHUNK, target.length - done), position + done);
    if (n <= 0) throw new StoreFormatError(`unexpected end of file at byte ${position + done}`);
    done += n;
  }
}

function readU64Array(fd: number, position: number, count: number): Float64Array {
  const raw = new BigUint64Array(count);
  readExact(fd, new Uint8Array(raw.buffer), position);
  const out = new Float64Array(count);
  for (let i = 0; i < count; i++) out[i] = Number(raw[i]!);
  return out;
}

export interface ExceptionRun {
  position: number;
  length: number;
  code: number;
}

export class SequenceStore {
  readonly path: string;
  readonly header: StoreHeader;
  readonly alphabet: Alphabet;
  readonly residueResidency: Residency;
  readonly headerResidency: Residency;
  private fd: number;
  private readonly sequenceOffsets: Float64Array;
  private readonly headerOffsets: Float64Array;
  private readonly exceptionIndex: Float64Array;
  private readonly exceptionPositions: Float64Array;
  private readonly exceptionLengths: Uint32Array;
  private readonly exceptionCodes: Uint8Array;
  private readonly idIndex: Uint32Array;
  private readonly residueData: Uint8Array | null;
  private readonly headerData: Uint8Array | null;

  private constructor(path: string, fd: number, options: OpenOptions) {
    this.path = path;
    this.fd = fd;
    const headerBuffer = Buffer.alloc(HEADER_SIZE);
    const fileSize = fstatSync(fd).size;
    if (fileSize < HEADER_SIZE) throw new StoreFormatError("file too small to be a HomologFinder store");
    readExact(fd, headerBuffer, 0);
    const header = decodeHeader(headerBuffer);
    if (header.headersOffset + header.headersLength > fileSize) throw new StoreFormatError("store file is truncated");
    this.header = header;
    this.alphabet = alphabetById(header.alphabetId);
    const count = header.sequenceCount;

    this.sequenceOffsets = readU64Array(fd, header.sequenceOffsetsOffset, count + 1);
    this.headerOffsets = readU64Array(fd, header.headerOffsetsOffset, count + 1);
    this.exceptionIndex = readU64Array(fd, header.exceptionIndexOffset, count + 1);

    const runs = header.exceptionCount;
    this.exceptionPositions = new Float64Array(runs);
    this.exceptionLengths = new Uint32Array(runs);
    this.exceptionCodes = new Uint8Array(runs);
    if (runs > 0) {
      const raw = Buffer.alloc(runs * EXCEPTION_ENTRY_SIZE);
      readExact(fd, raw, header.exceptionsOffset);
      for (let i = 0; i < runs; i++) {
        const base = i * EXCEPTION_ENTRY_SIZE;
        this.exceptionPositions[i] = Number(raw.readBigUInt64LE(base));
        this.exceptionLengths[i] = raw.readUInt32LE(base + 8);
        this.exceptionCodes[i] = raw[base + 12]!;
      }
    }

    this.idIndex = new Uint32Array(count * 2);
    readExact(fd, new Uint8Array(this.idIndex.buffer), header.idIndexOffset);

    this.residueResidency = options.residues ?? "memory";
    this.headerResidency = options.headers ?? "disk";
    if (this.residueResidency === "memory") {
      const backing = options.shared
        ? new SharedArrayBuffer(header.residuesLength)
        : new ArrayBuffer(header.residuesLength);
      this.residueData = new Uint8Array(backing);
      readExact(fd, this.residueData, header.residuesOffset);
    } else {
      this.residueData = null;
    }
    if (this.headerResidency === "memory") {
      this.headerData = new Uint8Array(header.headersLength);
      readExact(fd, this.headerData, header.headersOffset);
    } else {
      this.headerData = null;
    }
  }

  static open(path: string, options: OpenOptions = {}): SequenceStore {
    const fd = openSync(path, "r");
    try {
      return new SequenceStore(path, fd, options);
    } catch (error) {
      closeSync(fd);
      throw error;
    }
  }

  get count(): number {
    return this.header.sequenceCount;
  }

  get totalResidues(): number {
    return this.header.totalResidues;
  }

  get residueBuffer(): Uint8Array | null {
    return this.residueData;
  }

  close(): void {
    if (this.fd >= 0) {
      closeSync(this.fd);
      this.fd = -1;
    }
  }

  private checkIndex(index: number): void {
    if (!Number.isInteger(index) || index < 0 || index >= this.count) {
      throw new RangeError(`sequence index ${index} out of range 0..${this.count - 1}`);
    }
  }

  private ensureOpen(): number {
    if (this.fd < 0) throw new Error("store is closed");
    return this.fd;
  }

  length(index: number): number {
    this.checkIndex(index);
    return this.sequenceOffsets[index + 1]! - this.sequenceOffsets[index]!;
  }

  offset(index: number): number {
    this.checkIndex(index);
    return this.sequenceOffsets[index]!;
  }

  headerBytes(index: number): Uint8Array {
    this.checkIndex(index);
    const start = this.headerOffsets[index]!;
    const end = this.headerOffsets[index + 1]!;
    if (this.headerData) return this.headerData.subarray(start, end);
    const out = new Uint8Array(end - start);
    readExact(this.ensureOpen(), out, this.header.headersOffset + start);
    return out;
  }

  headerLine(index: number): string {
    return Buffer.from(this.headerBytes(index)).toString("utf8");
  }

  id(index: number): string {
    return splitHeader(this.headerLine(index)).id;
  }

  description(index: number): string {
    return splitHeader(this.headerLine(index)).description;
  }

  codes(index: number, target?: Uint8Array): Uint8Array {
    const length = this.length(index);
    const start = this.sequenceOffsets[index]!;
    if (this.alphabet.kind === "protein") {
      if (this.residueData) {
        const view = this.residueData.subarray(start, start + length);
        if (!target) return view;
        target.set(view);
        return target.subarray(0, length);
      }
      const out = target ? target.subarray(0, length) : new Uint8Array(length);
      readExact(this.ensureOpen(), out, this.header.residuesOffset + start);
      return out;
    }
    return this.unpackDna(index, start, length, target);
  }

  sequence(index: number): string {
    return decodeToString(this.alphabet, this.codes(index));
  }

  exceptions(index: number): ExceptionRun[] {
    this.checkIndex(index);
    const runs: ExceptionRun[] = [];
    for (let e = this.exceptionIndex[index]!; e < this.exceptionIndex[index + 1]!; e++) {
      runs.push({ position: this.exceptionPositions[e]!, length: this.exceptionLengths[e]!, code: this.exceptionCodes[e]! });
    }
    return runs;
  }

  findById(id: string): number[] {
    const hash = fnv1a32(Buffer.from(id, "utf8"));
    const index = this.idIndex;
    let lo = 0;
    let hi = this.count;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      if (index[mid * 2]! < hash) lo = mid + 1;
      else hi = mid;
    }
    const hits: number[] = [];
    for (let i = lo; i < this.count && index[i * 2] === hash; i++) {
      const candidate = index[i * 2 + 1]!;
      if (this.id(candidate) === id) hits.push(candidate);
    }
    return hits;
  }

  private unpackDna(index: number, start: number, length: number, target?: Uint8Array): Uint8Array {
    const out = target ? target.subarray(0, length) : new Uint8Array(length);
    if (length === 0) return out;
    const firstByte = Math.floor(start / 4);
    const lastByte = Math.floor((start + length - 1) / 4);
    let packed: Uint8Array;
    if (this.residueData) {
      packed = this.residueData.subarray(firstByte, lastByte + 1);
    } else {
      packed = new Uint8Array(lastByte - firstByte + 1);
      readExact(this.ensureOpen(), packed, this.header.residuesOffset + firstByte);
    }
    let slot = start - firstByte * 4;
    let byteIndex = 0;
    let current = packed[0]! >> (slot * 2);
    for (let i = 0; i < length; i++) {
      out[i] = current & 3;
      current >>= 2;
      if (++slot === 4) {
        slot = 0;
        current = packed[++byteIndex]!;
      }
    }
    for (let e = this.exceptionIndex[index]!; e < this.exceptionIndex[index + 1]!; e++) {
      const from = this.exceptionPositions[e]! - start;
      out.fill(this.exceptionCodes[e]!, from, from + this.exceptionLengths[e]!);
    }
    return out;
  }
}

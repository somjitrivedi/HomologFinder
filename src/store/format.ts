export const MAGIC = 0x42444648;
export const FORMAT_VERSION = 1;
export const HEADER_SIZE = 128;
export const EXCEPTION_ENTRY_SIZE = 16;
export const ID_ENTRY_SIZE = 8;

export interface StoreHeader {
  version: number;
  alphabetId: number;
  sequenceCount: number;
  totalResidues: number;
  maxLength: number;
  residuesOffset: number;
  residuesLength: number;
  sequenceOffsetsOffset: number;
  headerOffsetsOffset: number;
  exceptionIndexOffset: number;
  exceptionsOffset: number;
  exceptionCount: number;
  idIndexOffset: number;
  headersOffset: number;
  headersLength: number;
}

const U64_FIELDS = [
  ["sequenceCount", 8],
  ["totalResidues", 16],
  ["maxLength", 24],
  ["residuesOffset", 32],
  ["residuesLength", 40],
  ["sequenceOffsetsOffset", 48],
  ["headerOffsetsOffset", 56],
  ["exceptionIndexOffset", 64],
  ["exceptionsOffset", 72],
  ["exceptionCount", 80],
  ["idIndexOffset", 88],
  ["headersOffset", 96],
  ["headersLength", 104],
] as const;

export function encodeHeader(header: StoreHeader): Buffer {
  const buffer = Buffer.alloc(HEADER_SIZE);
  buffer.writeUInt32LE(MAGIC, 0);
  buffer.writeUInt16LE(header.version, 4);
  buffer.writeUInt8(header.alphabetId, 6);
  for (const [field, offset] of U64_FIELDS) buffer.writeBigUInt64LE(BigInt(header[field]), offset);
  return buffer;
}

export class StoreFormatError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "StoreFormatError";
  }
}

export function decodeHeader(buffer: Buffer): StoreHeader {
  if (buffer.length < HEADER_SIZE) throw new StoreFormatError("file too small to be a HomologFinder store");
  if (buffer.readUInt32LE(0) !== MAGIC) throw new StoreFormatError("not a HomologFinder store (bad magic)");
  const version = buffer.readUInt16LE(4);
  if (version !== FORMAT_VERSION) throw new StoreFormatError(`unsupported store version ${version}`);
  const header = { version, alphabetId: buffer.readUInt8(6) } as StoreHeader;
  for (const [field, offset] of U64_FIELDS) header[field] = Number(buffer.readBigUInt64LE(offset));
  return header;
}

export function fnv1a32(bytes: Uint8Array): number {
  let hash = 0x811c9dc5;
  for (let i = 0; i < bytes.length; i++) {
    hash ^= bytes[i]!;
    hash = Math.imul(hash, 0x01000193);
  }
  return hash >>> 0;
}

export function align8(value: number): number {
  return Math.ceil(value / 8) * 8;
}

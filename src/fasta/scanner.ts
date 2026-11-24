export class FastaFormatError extends Error {
  readonly line: number;
  constructor(message: string, line: number) {
    super(`${message} (line ${line})`);
    this.name = "FastaFormatError";
    this.line = line;
  }
}

export interface FastaHandlers {
  header(line: Uint8Array): void;
  residues(chunk: Uint8Array, start: number, end: number): void;
  recordEnd(): void;
}

const NL = 10;
const CR = 13;
const GT = 62;
const SEMI = 59;

function isBlank(chunk: Uint8Array, start: number, end: number): boolean {
  for (let i = start; i < end; i++) {
    const b = chunk[i]!;
    if (b !== 32 && b !== 9 && b !== CR && b !== 11 && b !== 12) return false;
  }
  return true;
}

export class FastaScanner {
  private readonly handlers: FastaHandlers;
  private atLineStart = true;
  private inHeader = false;
  private inComment = false;
  private haveRecord = false;
  private headerParts: Uint8Array[] = [];
  private headerLength = 0;
  private lineNumber = 1;
  private finished = false;

  constructor(handlers: FastaHandlers) {
    this.handlers = handlers;
  }

  get line(): number {
    return this.lineNumber;
  }

  push(chunk: Uint8Array): void {
    if (this.finished) throw new Error("scanner already finished");
    const length = chunk.length;
    let pos = 0;
    while (pos < length) {
      if (this.atLineStart) {
        const first = chunk[pos]!;
        this.atLineStart = false;
        if (first === GT) {
          if (this.haveRecord) this.handlers.recordEnd();
          this.haveRecord = true;
          this.inHeader = true;
          this.headerParts = [];
          this.headerLength = 0;
          pos++;
          continue;
        }
        if (first === SEMI) {
          this.inComment = true;
        }
      }

      const nl = chunk.indexOf(NL, pos);
      const end = nl === -1 ? length : nl;

      if (this.inHeader) {
        if (end > pos) {
          this.headerParts.push(chunk.slice(pos, end));
          this.headerLength += end - pos;
        }
        if (nl !== -1) this.emitHeader();
      } else if (!this.inComment && end > pos) {
        if (this.haveRecord) {
          this.handlers.residues(chunk, pos, end);
        } else if (!isBlank(chunk, pos, end)) {
          throw new FastaFormatError("sequence data before the first '>' header", this.lineNumber);
        }
      }

      if (nl === -1) {
        pos = length;
      } else {
        pos = nl + 1;
        this.atLineStart = true;
        this.inComment = false;
        this.lineNumber++;
      }
    }
  }

  end(): void {
    if (this.finished) return;
    this.finished = true;
    if (this.inHeader) this.emitHeader();
    if (this.haveRecord) this.handlers.recordEnd();
  }

  private emitHeader(): void {
    let line: Uint8Array;
    if (this.headerParts.length === 1) {
      line = this.headerParts[0]!;
    } else {
      line = new Uint8Array(this.headerLength);
      let offset = 0;
      for (const part of this.headerParts) {
        line.set(part, offset);
        offset += part.length;
      }
    }
    let stop = line.length;
    while (stop > 0 && (line[stop - 1] === CR || line[stop - 1] === 32 || line[stop - 1] === 9)) stop--;
    this.handlers.header(stop === line.length ? line : line.subarray(0, stop));
    this.inHeader = false;
    this.headerParts = [];
    this.headerLength = 0;
  }
}

export function splitHeader(header: string): { id: string; description: string } {
  const match = /^\s*(\S*)\s*(.*)$/s.exec(header);
  return { id: match?.[1] ?? "", description: match?.[2] ?? "" };
}

export function headerId(header: Uint8Array): Uint8Array {
  let start = 0;
  while (start < header.length && (header[start] === 32 || header[start] === 9)) start++;
  let end = start;
  while (end < header.length && header[end] !== 32 && header[end] !== 9) end++;
  return header.subarray(start, end);
}

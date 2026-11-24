import { FastaScanner, splitHeader } from "./scanner.ts";
import { openFastaStream } from "./source.ts";

export interface FastaRecord {
  header: string;
  id: string;
  description: string;
  sequence: string;
}

class RecordCollector {
  readonly ready: FastaRecord[] = [];
  private header = "";
  private parts: Buffer[] = [];

  readonly scanner = new FastaScanner({
    header: (line) => {
      this.header = Buffer.from(line).toString("utf8");
      this.parts = [];
    },
    residues: (chunk, start, end) => {
      this.parts.push(Buffer.from(chunk.subarray(start, end)));
    },
    recordEnd: () => {
      const sequence = Buffer.concat(this.parts).toString("latin1").replace(/\s+/g, "");
      const { id, description } = splitHeader(this.header);
      this.ready.push({ header: this.header, id, description, sequence });
      this.parts = [];
    },
  });
}

export async function* parseFasta(source: AsyncIterable<Uint8Array | string>): AsyncGenerator<FastaRecord> {
  const collector = new RecordCollector();
  for await (const chunk of source) {
    collector.scanner.push(typeof chunk === "string" ? Buffer.from(chunk, "utf8") : chunk);
    while (collector.ready.length > 0) yield collector.ready.shift()!;
  }
  collector.scanner.end();
  while (collector.ready.length > 0) yield collector.ready.shift()!;
}

export function parseFastaText(text: string | Uint8Array): FastaRecord[] {
  const collector = new RecordCollector();
  collector.scanner.push(typeof text === "string" ? Buffer.from(text, "utf8") : text);
  collector.scanner.end();
  return collector.ready;
}

export function readFastaFile(path: string): AsyncGenerator<FastaRecord> {
  return parseFasta(openFastaStream(path));
}

export function formatFasta(header: string, sequence: string, width = 60): string {
  const lines = [`>${header}`];
  if (width <= 0) {
    lines.push(sequence);
  } else {
    for (let i = 0; i < sequence.length; i += width) lines.push(sequence.slice(i, i + width));
  }
  return `${lines.join("\n")}\n`;
}

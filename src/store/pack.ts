import { type AlphabetKind, type DetectionResult, alphabetByKind, detectAlphabet } from "../alphabet/alphabet.ts";
import { FastaScanner } from "../fasta/scanner.ts";
import { openFastaStream } from "../fasta/source.ts";
import { type InvalidResiduePolicy, type StoreStats, StoreWriter } from "./writer.ts";

export interface PackOptions {
  alphabet?: AlphabetKind | "auto";
  onInvalid?: InvalidResiduePolicy;
  detectionBytes?: number;
  onProgress?: (inputBytes: number, sequences: number) => void;
}

export interface PackResult extends StoreStats {
  alphabet: AlphabetKind;
  detection: DetectionResult | null;
  inputBytes: number;
  elapsedMs: number;
}

export async function packFastaSource(
  source: AsyncIterable<Uint8Array>,
  output: string,
  options: PackOptions = {},
): Promise<PackResult> {
  const started = performance.now();
  const iterator = source[Symbol.asyncIterator]();
  const pending: Uint8Array[] = [];
  let detection: DetectionResult | null = null;
  let kind: AlphabetKind;
  let exhausted = false;

  if (!options.alphabet || options.alphabet === "auto") {
    const want = options.detectionBytes ?? 1 << 20;
    let buffered = 0;
    while (buffered < want) {
      const next = await iterator.next();
      if (next.done) {
        exhausted = true;
        break;
      }
      pending.push(next.value);
      buffered += next.value.length;
    }
    detection = detectAlphabet(Buffer.concat(pending));
    kind = detection.kind;
  } else {
    kind = options.alphabet;
  }

  const writer = new StoreWriter(output, alphabetByKind(kind), { onInvalid: options.onInvalid ?? "error" });
  const scanner = new FastaScanner({
    header: (line) => writer.beginRecord(line),
    residues: (chunk, start, end) => writer.addResidues(chunk, start, end),
    recordEnd: () => writer.endRecord(),
  });

  let inputBytes = 0;
  const feed = (chunk: Uint8Array) => {
    inputBytes += chunk.length;
    scanner.push(chunk);
    options.onProgress?.(inputBytes, writer.count);
  };

  try {
    for (const chunk of pending) feed(chunk);
    pending.length = 0;
    if (!exhausted) {
      for (;;) {
        const next = await iterator.next();
        if (next.done) break;
        feed(next.value);
      }
    }
    scanner.end();
  } catch (error) {
    writer.abort();
    await iterator.return?.();
    throw error;
  }

  const stats = writer.finish();
  return { ...stats, alphabet: kind, detection, inputBytes, elapsedMs: performance.now() - started };
}

export function packFasta(input: string, output: string, options: PackOptions = {}): Promise<PackResult> {
  return packFastaSource(openFastaStream(input), output, options);
}

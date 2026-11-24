export type AlphabetKind = "protein" | "dna";

export const SKIP = -2;
export const INVALID = -1;

export interface Alphabet {
  readonly kind: AlphabetKind;
  readonly id: number;
  readonly symbols: string;
  readonly encode: Int16Array;
  readonly decode: Uint8Array;
}

export const PROTEIN_SYMBOLS = "ARNDCQEGHILKMFPSTWYVBZX*UOJ-";
export const DNA_SYMBOLS = "ACGTURYKMSWBDHVN-";

function buildEncodeTable(symbols: string): Int16Array {
  const table = new Int16Array(256).fill(INVALID);
  for (let code = 0; code < symbols.length; code++) {
    const upper = symbols.charCodeAt(code);
    table[upper] = code;
    const lower = symbols[code]!.toLowerCase().charCodeAt(0);
    if (lower !== upper) table[lower] = code;
  }
  for (const ws of [9, 10, 11, 12, 13, 32]) table[ws] = SKIP;
  return table;
}

function buildDecodeTable(symbols: string): Uint8Array {
  const table = new Uint8Array(symbols.length);
  for (let code = 0; code < symbols.length; code++) table[code] = symbols.charCodeAt(code);
  return table;
}

export const PROTEIN: Alphabet = {
  kind: "protein",
  id: 1,
  symbols: PROTEIN_SYMBOLS,
  encode: buildEncodeTable(PROTEIN_SYMBOLS),
  decode: buildDecodeTable(PROTEIN_SYMBOLS),
};

export const DNA: Alphabet = {
  kind: "dna",
  id: 2,
  symbols: DNA_SYMBOLS,
  encode: buildEncodeTable(DNA_SYMBOLS),
  decode: buildDecodeTable(DNA_SYMBOLS),
};

export function alphabetById(id: number): Alphabet {
  if (id === PROTEIN.id) return PROTEIN;
  if (id === DNA.id) return DNA;
  throw new Error(`unknown alphabet id ${id}`);
}

export function alphabetByKind(kind: AlphabetKind): Alphabet {
  return kind === "dna" ? DNA : PROTEIN;
}

export function decodeToString(alphabet: Alphabet, codes: Uint8Array): string {
  const bytes = new Uint8Array(codes.length);
  const decode = alphabet.decode;
  for (let i = 0; i < codes.length; i++) bytes[i] = decode[codes[i]!]!;
  return Buffer.from(bytes.buffer, bytes.byteOffset, bytes.byteLength).toString("latin1");
}

export interface DetectionResult {
  kind: AlphabetKind;
  residuesSampled: number;
  nucleotideFraction: number;
}

const NUCLEOTIDE_CORE = new Set([65, 67, 71, 84, 85, 78, 97, 99, 103, 116, 117, 110]);

export function detectAlphabet(sample: Uint8Array): DetectionResult {
  let atLineStart = true;
  let inHeader = false;
  let letters = 0;
  let nucleotide = 0;
  let nonDna = false;
  for (let i = 0; i < sample.length; i++) {
    const b = sample[i]!;
    if (b === 10) {
      atLineStart = true;
      inHeader = false;
      continue;
    }
    if (atLineStart) {
      atLineStart = false;
      if (b === 62 || b === 59) {
        inHeader = true;
        continue;
      }
    }
    if (inHeader) continue;
    const code = DNA.encode[b]!;
    if (code === SKIP) continue;
    letters++;
    if (NUCLEOTIDE_CORE.has(b)) nucleotide++;
    if (code === INVALID) nonDna = true;
  }
  const fraction = letters === 0 ? 0 : nucleotide / letters;
  const kind: AlphabetKind = letters > 0 && !nonDna && fraction >= 0.9 ? "dna" : "protein";
  return { kind, residuesSampled: letters, nucleotideFraction: fraction };
}

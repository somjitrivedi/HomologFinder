import { readFileSync } from "node:fs";
import BLOSUM45 from "./matrices/BLOSUM45.json" with { type: "json" };
import BLOSUM50 from "./matrices/BLOSUM50.json" with { type: "json" };
import BLOSUM62 from "./matrices/BLOSUM62.json" with { type: "json" };
import BLOSUM80 from "./matrices/BLOSUM80.json" with { type: "json" };
import BLOSUM90 from "./matrices/BLOSUM90.json" with { type: "json" };
import NUC44 from "./matrices/NUC.4.4.json" with { type: "json" };
import PAM30 from "./matrices/PAM30.json" with { type: "json" };
import PAM70 from "./matrices/PAM70.json" with { type: "json" };
import PAM250 from "./matrices/PAM250.json" with { type: "json" };

export interface MatrixData {
  name: string;
  alphabet: string;
  scores: number[][];
}

export class SubstitutionMatrix {
  readonly name: string;
  readonly alphabet: string;
  readonly size: number;
  readonly scores: Float64Array;
  readonly index: Int16Array;
  readonly integral: boolean;

  constructor(data: MatrixData) {
    const size = data.alphabet.length;
    if (data.scores.length !== size || data.scores.some((row) => row.length !== size)) {
      throw new Error(`matrix ${data.name} is not ${size}x${size}`);
    }
    if (new Set(data.alphabet.toUpperCase()).size !== size) {
      throw new Error(`matrix ${data.name} has duplicate symbols`);
    }
    this.name = data.name;
    this.alphabet = data.alphabet.toUpperCase();
    this.size = size;
    this.scores = new Float64Array(size * size);
    let integral = true;
    for (let i = 0; i < size; i++) {
      for (let j = 0; j < size; j++) {
        const value = data.scores[i]![j]!;
        if (!Number.isFinite(value)) throw new Error(`matrix ${data.name} has a non-finite score`);
        if (!Number.isInteger(value)) integral = false;
        this.scores[i * size + j] = value;
      }
    }
    this.integral = integral;
    this.index = new Int16Array(256).fill(-1);
    for (let i = 0; i < size; i++) {
      const upper = this.alphabet.charCodeAt(i);
      this.index[upper] = i;
      const lower = this.alphabet[i]!.toLowerCase().charCodeAt(0);
      this.index[lower] = i;
    }
  }

  has(symbol: string): boolean {
    return symbol.length === 1 && this.index[symbol.charCodeAt(0)]! >= 0;
  }

  score(a: string, b: string): number {
    const i = this.index[a.charCodeAt(0)]!;
    const j = this.index[b.charCodeAt(0)]!;
    if (i < 0 || j < 0 || a.length !== 1 || b.length !== 1) {
      throw new Error(`symbol pair '${a}','${b}' is not in matrix ${this.name}`);
    }
    return this.scores[i * this.size + j]!;
  }

  isSymmetric(): boolean {
    for (let i = 0; i < this.size; i++) {
      for (let j = i + 1; j < this.size; j++) {
        if (this.scores[i * this.size + j] !== this.scores[j * this.size + i]) return false;
      }
    }
    return true;
  }

  toData(): MatrixData {
    const scores: number[][] = [];
    for (let i = 0; i < this.size; i++) {
      scores.push(Array.from(this.scores.subarray(i * this.size, (i + 1) * this.size)));
    }
    return { name: this.name, alphabet: this.alphabet, scores };
  }
}

const BUILTIN: Record<string, MatrixData> = {
  BLOSUM45,
  BLOSUM50,
  BLOSUM62,
  BLOSUM80,
  BLOSUM90,
  PAM30,
  PAM70,
  PAM250,
  "NUC.4.4": NUC44,
};

const cache = new Map<string, SubstitutionMatrix>();

export function builtinMatrixNames(): string[] {
  return Object.keys(BUILTIN);
}

export function getMatrix(name: string): SubstitutionMatrix {
  const key = Object.keys(BUILTIN).find((candidate) => candidate.toUpperCase() === name.toUpperCase());
  if (!key) throw new Error(`unknown matrix '${name}', available: ${builtinMatrixNames().join(", ")}`);
  let matrix = cache.get(key);
  if (!matrix) {
    matrix = new SubstitutionMatrix(BUILTIN[key]!);
    cache.set(key, matrix);
  }
  return matrix;
}

export function matchMismatchMatrix(alphabet: string, match: number, mismatch: number): SubstitutionMatrix {
  const scores = Array.from(alphabet, (_, i) => Array.from(alphabet, (_, j) => (i === j ? match : mismatch)));
  return new SubstitutionMatrix({ name: `match${match}/mismatch${mismatch}`, alphabet, scores });
}

export function parseNcbiMatrix(text: string, name: string): SubstitutionMatrix {
  const lines = text
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter((line) => line.length > 0 && !line.startsWith("#"));
  if (lines.length === 0) throw new Error(`matrix ${name} is empty`);
  const columns = lines[0]!.split(/\s+/);
  const scores: number[][] = [];
  const rows: string[] = [];
  for (const line of lines.slice(1)) {
    const fields = line.split(/\s+/);
    rows.push(fields[0]!);
    scores.push(fields.slice(1).map(Number));
  }
  if (rows.join("") !== columns.join("")) throw new Error(`matrix ${name} row and column symbols differ`);
  return new SubstitutionMatrix({ name, alphabet: columns.join(""), scores });
}

export function loadMatrix(nameOrPath: string): SubstitutionMatrix {
  if (builtinMatrixNames().some((n) => n.toUpperCase() === nameOrPath.toUpperCase())) return getMatrix(nameOrPath);
  const text = readFileSync(nameOrPath, "utf8");
  if (nameOrPath.toLowerCase().endsWith(".json")) return new SubstitutionMatrix(JSON.parse(text) as MatrixData);
  return parseNcbiMatrix(text, nameOrPath);
}

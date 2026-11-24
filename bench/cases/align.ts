import { getMatrix } from "../../src/align/matrix.ts";
import { type AlignMode, align, alignScore } from "../../src/align/pairwise.ts";
import { STANDARD_AMINO_ACIDS, mulberry32, randomSequence } from "../../src/util/random.ts";
import { emit, peakRssMb, time } from "../report.ts";

const length = Number(process.argv[2] ?? 300);
const pairs = Number(process.argv[3] ?? 200);
const random = mulberry32(11);
const matrix = getMatrix("BLOSUM62");
const inputs = Array.from({ length: pairs }, () => [
  randomSequence(random, STANDARD_AMINO_ACIDS, length),
  randomSequence(random, STANDARD_AMINO_ACIDS, length),
]);
const cells = pairs * length * length;

function measure(mode: AlignMode, withTrace: boolean): number {
  const options = { mode, matrix, openGap: -12, extendGap: -1 };
  for (const [t, q] of inputs.slice(0, 10)) (withTrace ? align : alignScore)(t!, q!, options);
  const { ms } = time(() => {
    for (const [t, q] of inputs) (withTrace ? align : alignScore)(t!, q!, options);
  });
  return cells / (ms / 1000) / 1e6;
}

emit({
  case: "align",
  metrics: [
    { name: "local score-only", value: measure("local", false), unit: "MCUPS" },
    { name: "local + traceback", value: measure("local", true), unit: "MCUPS" },
    { name: "global score-only", value: measure("global", false), unit: "MCUPS" },
    { name: "global + traceback", value: measure("global", true), unit: "MCUPS" },
  ],
  peakRssMb: peakRssMb(),
  details: { length, pairs, matrix: matrix.name, openGap: -12, extendGap: -1 },
});

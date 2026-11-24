import type { Alignment } from "./pairwise.ts";

export function formatAlignment(
  alignment: Alignment,
  options: { width?: number; targetName?: string; queryName?: string } = {},
): string {
  const width = options.width ?? 60;
  const targetName = options.targetName ?? "target";
  const queryName = options.queryName ?? "query";
  const label = Math.max(targetName.length, queryName.length, 6);
  const s = alignment.stats;
  const lines: string[] = [
    `mode ${alignment.mode}  score ${alignment.score}  length ${s.length}`,
    `identities ${s.identities}/${s.length} (${(s.identity * 100).toFixed(1)}%)  positives ${s.positives}/${s.length}  gaps ${s.gaps}/${s.length} in ${s.gapOpens} opens`,
    `target ${alignment.targetStart + 1}-${alignment.targetEnd}  query ${alignment.queryStart + 1}-${alignment.queryEnd}  cigar ${alignment.cigar || "-"}`,
    "",
  ];
  let t = alignment.targetStart;
  let q = alignment.queryStart;
  const digits = String(Math.max(alignment.targetEnd, alignment.queryEnd)).length;
  for (let k = 0; k < alignment.alignedTarget.length; k += width) {
    const tSeg = alignment.alignedTarget.slice(k, k + width);
    const qSeg = alignment.alignedQuery.slice(k, k + width);
    const mSeg = alignment.midline.slice(k, k + width);
    const tCount = tSeg.replace(/-/g, "").length;
    const qCount = qSeg.replace(/-/g, "").length;
    lines.push(`${targetName.padEnd(label)} ${String(t + (tCount ? 1 : 0)).padStart(digits)} ${tSeg} ${t + tCount}`);
    lines.push(`${"".padEnd(label)} ${"".padStart(digits)} ${mSeg}`);
    lines.push(`${queryName.padEnd(label)} ${String(q + (qCount ? 1 : 0)).padStart(digits)} ${qSeg} ${q + qCount}`);
    lines.push("");
    t += tCount;
    q += qCount;
  }
  return lines.join("\n");
}

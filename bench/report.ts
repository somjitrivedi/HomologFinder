export interface BenchMetric {
  name: string;
  value: number;
  unit: string;
}

export interface BenchReport {
  case: string;
  metrics: BenchMetric[];
  peakRssMb: number;
  details?: Record<string, unknown>;
}

export function peakRssMb(): number {
  return process.resourceUsage().maxRSS / 1024;
}

export function emit(report: BenchReport): void {
  process.stdout.write(`BENCH_RESULT ${JSON.stringify(report)}\n`);
}

export function time<T>(fn: () => T): { result: T; ms: number } {
  const start = performance.now();
  const result = fn();
  return { result, ms: performance.now() - start };
}

import { spawnSync } from "node:child_process";
import { existsSync } from "node:fs";
import path from "node:path";

const root = path.resolve(import.meta.dirname, "..");
const oracleDir = path.join(root, "tools", "oracle");
const venv = path.join(oracleDir, ".venv");
const venvPython =
  process.platform === "win32" ? path.join(venv, "Scripts", "python.exe") : path.join(venv, "bin", "python");
const systemPython = process.env.PYTHON ?? (process.platform === "win32" ? "python" : "python3");

function run(command: string, args: string[]): void {
  console.log(`> ${path.basename(command)} ${args.map((a) => path.relative(root, a) || a).join(" ")}`);
  const result = spawnSync(command, args, { cwd: root, stdio: "inherit" });
  if (result.status !== 0) {
    console.error(`command failed with exit code ${result.status ?? result.error?.message}`);
    process.exit(1);
  }
}

if (!existsSync(venvPython)) {
  run(systemPython, ["-m", "venv", venv]);
  run(venvPython, ["-m", "pip", "install", "--quiet", "-r", path.join(oracleDir, "requirements.txt")]);
}

const steps = process.argv.slice(2);
const all = steps.length === 0;
if (all || steps.includes("matrices")) {
  run(venvPython, ["-I", path.join(oracleDir, "export_matrices.py"), path.join(root, "src", "align", "matrices")]);
}
if (all || steps.includes("fasta")) {
  run(venvPython, ["-I", path.join(oracleDir, "gen_fasta_fixtures.py"), path.join(root, "test", "fixtures", "fasta")]);
}
if (all || steps.includes("align")) {
  run(venvPython, [
    "-I",
    path.join(oracleDir, "gen_align_fixtures.py"),
    path.join(root, "test", "fixtures", "biopython", "align-cases.json"),
  ]);
}

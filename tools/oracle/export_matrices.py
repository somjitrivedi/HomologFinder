import json
import sys
from pathlib import Path

from Bio.Align import substitution_matrices

NAMES = [
    "BLOSUM45",
    "BLOSUM50",
    "BLOSUM62",
    "BLOSUM80",
    "BLOSUM90",
    "PAM30",
    "PAM70",
    "PAM250",
    "NUC.4.4",
]


def main(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in NAMES:
        matrix = substitution_matrices.load(name)
        alphabet = "".join(matrix.alphabet)
        rows = []
        for a in alphabet:
            row = []
            for b in alphabet:
                value = float(matrix[a][b])
                row.append(int(value) if value.is_integer() else value)
            rows.append(row)
        payload = {"name": name, "alphabet": alphabet, "scores": rows}
        target = out_dir / f"{name}.json"
        target.write_text(json.dumps(payload, separators=(",", ":")) + "\n", encoding="utf-8", newline="\n")
        print(f"{name}: {len(alphabet)} symbols -> {target}")


if __name__ == "__main__":
    main(Path(sys.argv[1]))

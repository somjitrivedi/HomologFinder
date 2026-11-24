import gzip
import json
import sys
from pathlib import Path

from Bio import SeqIO


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def main(fixture_dir: Path) -> None:
    sources = sorted(p for p in fixture_dir.iterdir() if p.name.endswith((".fasta", ".fasta.gz")))
    for source in sources:
        records = []
        with open_text(source) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                parts = record.description.split(None, 1)
                records.append({
                    "id": record.id,
                    "description": parts[1] if len(parts) > 1 else "",
                    "sequence": str(record.seq),
                })
        stem = source.name.removesuffix(".gz").removesuffix(".fasta")
        target = fixture_dir / f"{stem}.expected.json"
        target.write_text(json.dumps(records, separators=(",", ":")) + "\n", encoding="utf-8", newline="\n")
        print(f"{source.name}: {len(records)} records -> {target.name}")


if __name__ == "__main__":
    main(Path(sys.argv[1]))

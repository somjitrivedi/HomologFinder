import json
import random
import sys
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices

CAP = 400
PROTEIN = "ARNDCQEGHILKMFPSTWYV"
PROTEIN_EXTRA = "BZX"
DNA = "ACGT"

PROTEIN_MATRICES = ["BLOSUM62", "BLOSUM45", "BLOSUM50", "BLOSUM80", "BLOSUM90", "PAM30", "PAM70", "PAM250"]
PROTEIN_GAPS = [(-12, -1), (-11, -1), (-10, -0.5), (-8, -2), (-5, -5), (-3, -1), (-1, -1), (-0.5, -0.5)]
DNA_SCORING = [
    {"match": 1, "mismatch": -1, "open": -2, "extend": -1},
    {"match": 2, "mismatch": -3, "open": -5, "extend": -2},
    {"match": 1, "mismatch": -2, "open": -2.5, "extend": -0.5},
    {"match": 5, "mismatch": -4, "open": -10, "extend": -1},
]


def random_seq(rng, alphabet, length):
    return "".join(rng.choice(alphabet) for _ in range(length))


def mutate(rng, seq, alphabet, rate):
    out = []
    for ch in seq:
        r = rng.random()
        if r < rate * 0.6:
            out.append(rng.choice(alphabet))
        elif r < rate * 0.8:
            continue
        elif r < rate:
            out.append(ch)
            out.append(random_seq(rng, alphabet, rng.randint(1, 4)))
        else:
            out.append(ch)
    return "".join(out) or rng.choice(alphabet)


def make_pair(rng, alphabet, max_len):
    kind = rng.random()
    if kind < 0.35:
        return random_seq(rng, alphabet, rng.randint(1, max_len)), random_seq(rng, alphabet, rng.randint(1, max_len))
    base = random_seq(rng, alphabet, rng.randint(1, max_len))
    other = mutate(rng, base, alphabet, rng.choice([0.05, 0.15, 0.3, 0.5]))
    if kind > 0.8:
        flank = random_seq(rng, alphabet, rng.randint(1, 15))
        other = flank + other + random_seq(rng, alphabet, rng.randint(0, 15))
    return base, other


def run_case(aligner, target, query):
    score = aligner.score(target, query)
    alignments = aligner.align(target, query)
    try:
        count = len(alignments)
    except OverflowError:
        count = None
    complete = count is not None and count <= CAP
    coords = []
    if complete:
        for alignment in alignments:
            coords.append(alignment.coordinates.tolist())
    else:
        coords.append(next(iter(alignments)).coordinates.tolist())
    return {
        "score": float(score),
        "count": count,
        "complete": complete,
        "alignments": coords,
    }


def main(out_path: Path, seed: int) -> None:
    rng = random.Random(seed)
    cases = []
    for matrix_name in PROTEIN_MATRICES:
        matrix = substitution_matrices.load(matrix_name)
        for open_gap, extend_gap in PROTEIN_GAPS:
            for mode in ("global", "local"):
                aligner = Align.PairwiseAligner(mode=mode, substitution_matrix=matrix,
                                                open_gap_score=open_gap, extend_gap_score=extend_gap)
                for _ in range(6):
                    alphabet = PROTEIN + (PROTEIN_EXTRA if rng.random() < 0.2 else "")
                    target, query = make_pair(rng, alphabet, 60)
                    result = run_case(aligner, target, query)
                    if mode == "local" and result["score"] <= 0:
                        continue
                    cases.append({
                        "id": f"{matrix_name}/{open_gap}/{extend_gap}/{mode}/{len(cases)}",
                        "mode": mode,
                        "matrix": matrix_name,
                        "open": open_gap,
                        "extend": extend_gap,
                        "target": target,
                        "query": query,
                        **result,
                    })
    for scoring in DNA_SCORING:
        for mode in ("global", "local"):
            aligner = Align.PairwiseAligner(mode=mode, match_score=scoring["match"],
                                            mismatch_score=scoring["mismatch"],
                                            open_gap_score=scoring["open"], extend_gap_score=scoring["extend"])
            for _ in range(25):
                target, query = make_pair(rng, DNA, 50)
                result = run_case(aligner, target, query)
                if mode == "local" and result["score"] <= 0:
                    continue
                cases.append({
                    "id": f"dna/{scoring['match']}/{scoring['mismatch']}/{scoring['open']}/{scoring['extend']}/{mode}/{len(cases)}",
                    "mode": mode,
                    "match": scoring["match"],
                    "mismatch": scoring["mismatch"],
                    "open": scoring["open"],
                    "extend": scoring["extend"],
                    "target": target,
                    "query": query,
                    **result,
                })
    matrix = substitution_matrices.load("NUC.4.4")
    for mode in ("global", "local"):
        aligner = Align.PairwiseAligner(mode=mode, substitution_matrix=matrix, open_gap_score=-10, extend_gap_score=-1)
        for _ in range(25):
            target, query = make_pair(rng, "ACGTN" if rng.random() < 0.2 else DNA, 60)
            result = run_case(aligner, target, query)
            if mode == "local" and result["score"] <= 0:
                continue
            cases.append({
                "id": f"NUC.4.4/-10/-1/{mode}/{len(cases)}",
                "mode": mode,
                "matrix": "NUC.4.4",
                "open": -10,
                "extend": -1,
                "target": target,
                "query": query,
                **result,
            })
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"biopython": __import__("Bio").__version__, "seed": seed, "cases": cases}
    out_path.write_text(json.dumps(payload, separators=(",", ":")) + "\n", encoding="utf-8", newline="\n")
    complete = sum(1 for c in cases if c["complete"])
    print(f"{len(cases)} alignment cases ({complete} with complete optimal sets) -> {out_path}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 20261124)

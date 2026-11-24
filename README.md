# HomologFinder

A sequence alignment and homology search engine in TypeScript (a mini-BLAST). Give it a DNA or protein sequence and it finds similar sequences in a large database using k-mer indexing and Smith-Waterman alignment.

Status: **Phases 0-2 complete** (setup and benchmark harness, sequence storage, alignment algorithms). The full roadmap is at the bottom.

## Requirements

- Node.js 22.18 or newer. TypeScript runs directly through Node's built-in type stripping, so no build step is needed during development.
- Python 3.10+ only to regenerate the Biopython test fixtures. It is not needed to run or test the engine.
- No native dependencies.

```sh
npm install
npm run data:swissprot      # download UniProt Swiss-Prot (~94 MB) into data/
npm test                    # 83 tests: parser, store, matrices, alignment, Biopython parity, CLI
npm run typecheck
npm run bench               # pack / open / align benchmarks, appended to bench/results/history.jsonl
npm run build               # optional: compile to dist/ (bin: homologfinder)
```

## CLI

```sh
node src/cli.ts pack data/uniprot_sprot.fasta.gz -o data/swissprot.hfdb
node src/cli.ts info data/swissprot.hfdb
node src/cli.ts get  data/swissprot.hfdb "sp|P69905|HBA_HUMAN" "#42"
node src/cli.ts align --db data/swissprot.hfdb "sp|P69905|HBA_HUMAN" "sp|P68871|HBB_HUMAN"
node src/cli.ts align HEAGAWGHEE PAWHEAE --mode global --matrix BLOSUM62 --open -10 --extend -1
node src/cli.ts align ACGTTGCA ACGTGCA --match 1 --mismatch -2 --open -3 --extend -1 --json
node src/cli.ts matrices
```

`pack` reads plain or gzip FASTA (detected by magic bytes) and works out whether the input is protein or DNA unless `--alphabet` is given. Invalid residues stop the pack and report the record id; `--on-invalid replace` swaps them for `X` (protein) or `N` (DNA) and counts them instead.

## Phase 1: sequence storage

### Streaming parser

`FastaScanner` is a byte-level state machine that takes chunks of any size, down to 1 byte. Records are never assembled as strings on the packing path. It handles CRLF line endings, blank lines, `;` comment lines, empty records, a missing final newline, and headers split across chunks. Its output matches Biopython `SeqIO` on every test fixture.

### Encodings

| Alphabet | Symbols | Storage |
|---|---|---|
| protein | `ARNDCQEGHILKMFPSTWYVBZX*UOJ-` (first 24 in NCBI matrix order) | 1 byte per residue |
| DNA | `ACGT` packed, plus `URYKMSWBDHVN-` as exceptions | 2 bits per base, 4 bases per byte, plus a run-length table of ambiguity codes |

Lowercase is folded to uppercase. Ambiguity runs (for example a stretch of `N`) are stored once as `(position, length, symbol)`. The packed bases under a run are zero, and the run is laid back over them when the sequence is read.

### `.hfdb` file format (v1, little-endian)

| Section | Contents |
|---|---|
| header (128 B) | magic `HFDB`, version, alphabet, sequence count, total residues, longest sequence, offset and length of every section |
| residues | packed residue bytes, all sequences back to back |
| sequence offsets | `u64[count+1]`: residue offset of each sequence |
| header offsets | `u64[count+1]`: byte offset of each header line |
| exception index | `u64[count+1]`: first ambiguity run of each sequence |
| exceptions | 16 B per run: `u64 position, u32 length, u8 symbol` |
| id index | `(u32 FNV-1a hash of id, u32 sequence index)`, sorted for binary search |
| headers | raw header lines, UTF-8 |

The writer streams residues straight to disk and keeps headers in a temporary side file. It writes the tables at the end, then renames `x.hfdb.tmp` to `x.hfdb`, so a failed pack never leaves a partial store behind. Memory use stays flat: only the offset tables grow.

`SequenceStore.open(path, { residues, headers, shared })`:
- `residues: "memory"` loads the residue section into one buffer. With `shared: true` that buffer is a `SharedArrayBuffer`, ready for worker threads in Phase 7.
- `residues: "disk"` serves sequences with positional reads, which lets the OS page cache act like a memory map.
- Headers are read from disk by default.

## Phase 2: alignment algorithms

`align(target, query, { mode, matrix, openGap, extendGap })` runs Needleman-Wunsch (global) or Smith-Waterman (local) with affine gaps, using Gotoh's three-state recurrence. `alignScore` and `alignEncoded` compute the same score in linear memory, without a traceback.

- **Gap convention (same as Biopython):** a gap of length L scores `openGap + (L-1) * extendGap`. BLAST's "existence 11, extension 1" equals `--open -12 --extend -1`, which is the CLI default for protein.
- **Matrices:** BLOSUM45/50/62/80/90, PAM30/70/250, NUC.4.4, exported from Biopython 1.88. Custom matrices can be loaded from NCBI text or JSON files, and `matchMismatchMatrix` builds simple match/mismatch scoring.
- **Output:** score, start and end positions, the aligned strings, a midline, a CIGAR string, Biopython-style `coordinates`, and identity, positive, gap and gap-open counts.
- **Exact Biopython parity:** when several alignments tie for the best score, HomologFinder applies Biopython's tie-breaking rules, so it returns the same first alignment Biopython does:
  - Affine gaps: Gotoh traceback with the preference order M, then Ix, then Iy.
  - Linear gaps (`open == extend`), where Biopython switches to plain NW/SW: horizontal, then vertical, then diagonal, plus Biopython's rule against zero-score extensions past a local endpoint.
- **Memory guard:** a traceback alignment over more than `maxCells` cells (default 2^28) throws a `RangeError`. Score-only alignment has no such limit.

### Verification

`test/fixtures/biopython/align-cases.json` holds 1,017 cases generated by `tools/oracle/gen_align_fixtures.py`. They cover 8 protein matrices × 8 gap schemes × global and local modes, plus DNA match/mismatch and NUC.4.4 cases. There are random pairs, mutated homolog pairs and flanked pairs, with B/Z/X and N included. For every case:

- the score equals Biopython's exactly, both with traceback and score-only;
- the returned alignment re-scores to that optimum;
- the alignment is a member of Biopython's full set of optimal alignments;
- the alignment is **identical to Biopython's first alignment** (1,017 of 1,017).

Regenerate the fixtures with `npm run oracle`. This creates `tools/oracle/.venv` with Biopython 1.88 on first use.

## Benchmarks

Measured with `npm run bench` on an AMD Ryzen 5 5625U (6 cores / 12 threads), Node 24.18, Windows 11:

| Case | Result | Peak RSS |
|---|---|---|
| Pack Swiss-Prot from .gz (576,073 proteins, 209.2M residues, 288 MB of FASTA text) | **3.8 s** (75 MB/s), 302 MB store | 148 MB |
| Open the store, residues in memory | 175 ms; random sequence access 0.2 µs; scan of all sequences 63 ms | 290 MB |
| Open the store, residues on disk | 63 ms; random sequence access 10 µs | 100 MB |
| Look up a sequence by id | 14 µs | |
| Smith-Waterman score-only, 300×300 BLOSUM62 | 50 MCUPS (million DP cells per second) | 89 MB |
| Smith-Waterman with traceback | 47 MCUPS | |

The plain TypeScript alignment speed (about 50 MCUPS) is the baseline that the Phase 3 WebAssembly SIMD kernel will be measured against.

## Layout

```
src/alphabet     protein and DNA alphabets, encoding tables, auto-detection
src/fasta        streaming scanner, record parser, gzip-aware file source
src/store        .hfdb format, streaming writer, reader, pack pipeline
src/align        substitution matrices, Gotoh/NW/SW, linear-gap traceback, text formatting
src/cli.ts       command-line interface
bench/           benchmark harness (each case runs in its own process for clean peak-RSS numbers)
scripts/         Swiss-Prot download, synthetic DNA generator, sample builder, oracle runner
tools/oracle     Biopython fixture generators (dev only)
test/            node:test suites and fixtures
```

## Data

`test/fixtures/fasta/swissprot-sample.fasta` holds 220 records from UniProtKB/Swiss-Prot (The UniProt Consortium), used under CC BY 4.0. `dna-sample.fasta` is synthetic (`npm run data:synth-dna`). The full Swiss-Prot download lives in `data/`, which is not committed.

## Roadmap

| # | Phase | Status |
|---|---|---|
| 0 | Setup and benchmark harness | done |
| 1 | Sequence storage | done |
| 2 | Alignment algorithms | done |
| 3 | Speed kernel: AssemblyScript to WASM SIMD, striped Smith-Waterman (Farrar), banded alignment | next |
| 4 | K-mer / neighbourhood-word index | |
| 5 | Seed-and-extend search + CLI + BLAST tabular output | |
| 6 | Karlin-Altschul statistics, SEG/DUST masking, accuracy against BLAST+ | |
| 7 | Parallel (worker_threads) and distributed search | |
| 8 | Job engine and API | |
| 9 | Alignment viewer | |
| 10 | Multiple alignment and phylogenetic trees (stretch) | |

BLAST+ and parasail are not installed yet. They become the comparison references in Phases 6 and 3 respectively.

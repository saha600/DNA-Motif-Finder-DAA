# DNA Motif Finder with Pattern Variants (DAA PBL 9)

A Python project that finds exact DNA motifs in synthetic DNA sequences and compares **Knuth–Morris–Pratt (KMP)** with **Rabin–Karp**.

## Features
- Accepts DNA sequences containing only A, T, G and C.
- Generates reproducible synthetic DNA.
- Searches multiple motifs of different lengths.
- Reports all zero-based starting positions, including overlapping matches.
- Compares search time and total matches for both algorithms.
- Includes boundary, repeated-character, absent-pattern and validation tests.

## Requirements
- Python 3.10 or newer
- No third-party packages required

## Run
```bash
python app.py
```
Or run the demonstration:
```bash
python dna_motif_finder.py
```

## Run tests
```bash
python -m unittest discover -s tests -v
```

## Algorithms and complexity
| Algorithm | Preprocessing | Search | Extra space |
|---|---:|---:|---:|
| KMP | O(m) | O(n) | O(m) |
| Rabin–Karp | O(m) | Average O(n+m); worst O(nm) | O(1) rolling-hash space |

Here, `n` is the DNA sequence length and `m` is the motif length. Rabin–Karp verifies characters after a hash match, so hash collisions do not create false matches.

## Example
Sequence: `ATGCGATGCGATGCGTTA`  
Motifs: `ATGC`, `GCG`, `TT`, `AAAA`

The program prints each motif's match positions for KMP and Rabin–Karp, then measures elapsed search time. Timing is machine- and input-dependent; run the benchmark on your own computer and report the observed values.

## GitHub upload
1. Create a new public repository named `DNA-Motif-Finder-DAA`.
2. Upload `dna_motif_finder.py`, `app.py`, `README.md`, and the `tests` folder.
3. Add your report and PPT if your faculty wants them in the same repository.
4. Copy the repository URL and submit it.

## Academic note
This project performs exact string matching on synthetic DNA strings. It is an educational algorithm demonstration, not a clinical or biological diagnosis tool.

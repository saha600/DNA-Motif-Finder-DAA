"""
DNA Motif Finder with Pattern Variants
DAA PBL 9 - Exact DNA motif matching using KMP and Rabin-Karp.
Coordinates returned by both algorithms are zero-based.
"""
from __future__ import annotations
import random
import time
from typing import Iterable

DNA_ALPHABET = "ATGC"

def validate_dna(sequence: str) -> str:
    sequence = sequence.upper().strip()
    invalid = sorted(set(sequence) - set(DNA_ALPHABET))
    if invalid:
        raise ValueError(f"Invalid DNA character(s): {', '.join(invalid)}. Use only A, T, G, C.")
    return sequence

def generate_dna(length: int, seed: int | None = 42) -> str:
    """Generate reproducible synthetic DNA of the requested length."""
    if length < 0:
        raise ValueError("Length cannot be negative.")
    rng = random.Random(seed)
    return "".join(rng.choice(DNA_ALPHABET) for _ in range(length))

def build_lps(pattern: str) -> list[int]:
    """Build the Longest Proper Prefix which is also Suffix array."""
    m = len(pattern)
    lps = [0] * m
    length = 0
    i = 1
    while i < m:
        if pattern[i] == pattern[length]:
            length += 1
            lps[i] = length
            i += 1
        elif length:
            length = lps[length - 1]
        else:
            lps[i] = 0
            i += 1
    return lps

def kmp_search(text: str, pattern: str) -> list[int]:
    """Return every exact occurrence start index using KMP; overlaps are included."""
    text, pattern = validate_dna(text), validate_dna(pattern)
    if not pattern:
        raise ValueError("Pattern cannot be empty.")
    n, m = len(text), len(pattern)
    lps = build_lps(pattern)
    matches, i, j = [], 0, 0
    while i < n:
        if text[i] == pattern[j]:
            i += 1
            j += 1
            if j == m:
                matches.append(i - j)
                j = lps[j - 1]
        elif j:
            j = lps[j - 1]
        else:
            i += 1
    return matches

def rabin_karp_search(text: str, pattern: str, prime: int = 101) -> list[int]:
    """Return every exact occurrence using a rolling hash; verify characters on hash hits."""
    text, pattern = validate_dna(text), validate_dna(pattern)
    if not pattern:
        raise ValueError("Pattern cannot be empty.")
    n, m = len(text), len(pattern)
    if m > n:
        return []
    d, h = 256, 1
    for _ in range(m - 1):
        h = (h * d) % prime
    p_hash = t_hash = 0
    for i in range(m):
        p_hash = (d * p_hash + ord(pattern[i])) % prime
        t_hash = (d * t_hash + ord(text[i])) % prime
    matches = []
    for i in range(n - m + 1):
        if p_hash == t_hash and text[i:i + m] == pattern:
            matches.append(i)
        if i < n - m:
            t_hash = (d * (t_hash - ord(text[i]) * h) + ord(text[i + m])) % prime
            t_hash %= prime
    return matches

def search_variants(sequence: str, patterns: Iterable[str], algorithm: str = "kmp") -> dict[str, list[int]]:
    """Search patterns of different lengths; dictionary maps each pattern to start positions."""
    sequence = validate_dna(sequence)
    searcher = {"kmp": kmp_search, "rabin-karp": rabin_karp_search}.get(algorithm.lower())
    if searcher is None:
        raise ValueError("Algorithm must be 'kmp' or 'rabin-karp'.")
    results = {}
    for pattern in patterns:
        pattern = validate_dna(pattern)
        results[pattern] = searcher(sequence, pattern)
    return results

def benchmark(sequence: str, patterns: Iterable[str]) -> dict[str, dict]:
    """Measure search-only wall-clock time. Results vary by machine and input."""
    output = {}
    for name, func in (("KMP", kmp_search), ("Rabin-Karp", rabin_karp_search)):
        start = time.perf_counter()
        total = 0
        for pattern in patterns:
            total += len(func(sequence, pattern))
        elapsed = time.perf_counter() - start
        output[name] = {"seconds": elapsed, "total_matches": total}
    return output

if __name__ == "__main__":
    dna = "ATGCGATGCGATGCGTTA"
    motifs = ["ATGC", "GCG", "TT", "AAAA"]
    print("DNA sequence:", dna)
    print("KMP results:", search_variants(dna, motifs, "kmp"))
    print("Rabin-Karp results:", search_variants(dna, motifs, "rabin-karp"))
    print("Benchmark:", benchmark(dna, motifs))

"""Simple interactive command-line interface for DNA Motif Finder."""
from dna_motif_finder import generate_dna, search_variants, benchmark, validate_dna

def main():
    print("=" * 58)
    print("          DNA MOTIF FINDER - DAA PBL 9")
    print("=" * 58)
    choice = input("Use (1) sample DNA or (2) generate synthetic DNA? ").strip()
    if choice == "2":
        length = int(input("Enter DNA length: "))
        sequence = generate_dna(length)
    else:
        sequence = "ATGCGATGCGATGCGTTA"
    print("\nDNA sequence:", sequence)
    raw = input("Enter DNA motifs separated by commas (e.g. ATG,GCG,TT): ")
    patterns = [p.strip().upper() for p in raw.split(",") if p.strip()]
    if not patterns:
        print("No motifs entered.")
        return
    algorithm = input("Choose algorithm (kmp/rabin-karp/both) [both]: ").strip().lower() or "both"
    algorithms = ["kmp", "rabin-karp"] if algorithm == "both" else [algorithm]
    for algo in algorithms:
        results = search_variants(sequence, patterns, algo)
        print(f"\n{algo.upper()} RESULTS")
        for motif, positions in results.items():
            print(f"  {motif} (length {len(motif)}): {positions or 'No match'}")
    if len(algorithms) == 2:
        print("\nPERFORMANCE COMPARISON")
        for name, data in benchmark(sequence, patterns).items():
            print(f"  {name}: {data['seconds']:.8f} seconds; {data['total_matches']} total matches")

if __name__ == "__main__":
    main()

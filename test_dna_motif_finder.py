import unittest
from dna_motif_finder import kmp_search, rabin_karp_search, search_variants, generate_dna

class TestDNAMotifFinder(unittest.TestCase):
    def test_basic_match(self):
        self.assertEqual(kmp_search("ATGCATGC", "ATGC"), [0, 4])
        self.assertEqual(rabin_karp_search("ATGCATGC", "ATGC"), [0, 4])

    def test_overlapping_matches(self):
        self.assertEqual(kmp_search("AAAAA", "AAA"), [0, 1, 2])
        self.assertEqual(rabin_karp_search("AAAAA", "AAA"), [0, 1, 2])

    def test_boundary_match(self):
        self.assertEqual(kmp_search("CCCATG", "ATG"), [3])
        self.assertEqual(rabin_karp_search("CCCATG", "ATG"), [3])

    def test_absent_and_long_pattern(self):
        self.assertEqual(kmp_search("ATGC", "GGG"), [])
        self.assertEqual(rabin_karp_search("ATGC", "GGG"), [])
        self.assertEqual(kmp_search("AT", "ATG"), [])

    def test_multiple_lengths(self):
        result = search_variants("ATGCGATG", ["AT", "GCG", "ATGC"], "kmp")
        self.assertEqual(result, {"AT": [0, 5], "GCG": [2], "ATGC": [0]})

    def test_invalid_character(self):
        with self.assertRaises(ValueError):
            kmp_search("ATGX", "AT")

    def test_empty_pattern(self):
        with self.assertRaises(ValueError):
            rabin_karp_search("ATGC", "")

    def test_generator(self):
        self.assertEqual(len(generate_dna(100)), 100)
        self.assertTrue(set(generate_dna(100)).issubset(set("ATGC")))

if __name__ == "__main__":
    unittest.main()

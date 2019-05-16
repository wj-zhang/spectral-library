import unittest
from find_protein_hits import find_protein_hits


class CompareTwoResultsCase(unittest.TestCase):
    """Tests for 'find_protein_hit.py'."""

    def setUp(self):
        self.id_file = '../tmp/library_id.csv'
        self.db_file = '../data/Human_UniProt_Isoforms_database_with_cRAP_and_spiked-in_proteins_08252015.fasta'
        self.hit_file = '../tmp/protein_hit.csv'

    def test_compare_two_results(self):
        find_protein_hits(self.id_file, self.db_file, self.hit_file)


if __name__ == '__main__':
    unittest.main()

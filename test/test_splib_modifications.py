import unittest
import splib_modifications


class SplibModificationsTestCase(unittest.TestCase):
    """Tests for 'splib_modifications.py'."""

    def setUp(self):
        self.peptide = 'n[43]AIDMC[160]PKNM[147]ASYYM[147]'

    def test_convert_modification(self):
        modified_seq = splib_modifications.convert_modification_forward(self.peptide)
        self.assertEqual(modified_seq, 'bAIDMeCPKNdfMASYYdfM')

        peptide = splib_modifications.convert_modification_backward(modified_seq)
        self.assertEqual(peptide, self.peptide)


if __name__ == '__main__':
    unittest.main()

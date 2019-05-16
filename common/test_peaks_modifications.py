import unittest
import peaks_modifications


class SplibModificationsTestCase(unittest.TestCase):
    """Tests for 'peaks_modifications.py'."""

    def setUp(self):
        self.peptide = 'DVAEC(+57.02)GPQQ(+.98)ELDLNSPR'

    def test_convert_modification(self):
        modified_seq = peaks_modifications.convert_modification_forward(self.peptide)
        self.assertEqual(modified_seq, 'DVAEeCGPQhQELDLNSPR')


if __name__ == '__main__':
    unittest.main()

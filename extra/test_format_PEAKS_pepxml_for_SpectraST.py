import unittest
import format_PEAKS_pepxml_for_SpectraST


class FormatPEAKSpepxmlForSpectraSTTestCase(unittest.TestCase):
    """Tests for 'format_PEAKS_pepxml_for_SpectraST.py'."""

    def setUp(self):
        self.in_file = '../data/PEAKS_peptides.pep.xml'
        self.out_file = '../tmp/PEAKS_peptides.ST.pep.xml'

    def test_format_PEAKS_pepxml_for_SpectraST(self):
        format_PEAKS_pepxml_for_SpectraST.run(self.in_file, self.out_file)


if __name__ == '__main__':
    unittest.main()

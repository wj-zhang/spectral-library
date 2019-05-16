import unittest
from pyteomics import mgf
import combine_clone_from_two_mgfs


class CombineCloneFromTwoMgfsTestCase(unittest.TestCase):
    """Tests for 'combine_clone_mgf.py'."""

    def setUp(self):
        self.single_ms_file = '../data/CharmeRT/search/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera.mgf'
        self.clone_ms_file = '../data/CharmeRT/search/HeLa_1ug_isow2_1hgradient_datasetA_rep1_chimera.mgf'
        self.ms_file = '../tmp/CharmeRT/HeLa_1ug_isow2_1hgradient_datasetA_rep1_combined.mgf'

    def test_feature_association(self):
        spectra = combine_clone_from_two_mgfs.combine(self.clone_ms_file, self.single_ms_file)

        with open(self.ms_file, 'w') as destination:
            mgf.write(spectra, destination)


if __name__ == '__main__':
    unittest.main()

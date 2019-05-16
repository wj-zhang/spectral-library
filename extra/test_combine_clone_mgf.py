import unittest
from pyteomics import mgf
import combine_clone_mgf


class CombineCloneMgfTestCase(unittest.TestCase):
    """Tests for 'combine_clone_mgf.py'."""

    def setUp(self):
        self.ms_file = '../tmp/PRG2017/corrected_mass.mgf'

    def test_feature_association(self):
        spectra = combine_clone_mgf.combine(self.ms_file)

        with open(self.ms_file, 'w') as destination:
            mgf.write(spectra, destination)


if __name__ == '__main__':
    unittest.main()

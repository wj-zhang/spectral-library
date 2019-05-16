import unittest
import estimate_mass_shift_from_mgf


class EstimateMassShiftFromMgfTestCase(unittest.TestCase):
    """Tests for 'estimate_mass_shift_from_mgf.py'."""

    def setUp(self):
        self.mgf_file = '../data/workflow/34 200 pg A549 target 20 (15.6ppm)_20190107_1_frac_cali_FDR_0.1.mgf'

    def test_estimate_mass_shift(self):
        estimate_mass_shift_from_mgf.estimate_mass_shift(self.mgf_file)


if __name__ == '__main__':
    unittest.main()

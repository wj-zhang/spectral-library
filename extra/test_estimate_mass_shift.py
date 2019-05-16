import unittest
import estimate_mass_shift


class EstimateMassShiftTestCase(unittest.TestCase):
    """Tests for 'estimate_mass_shift.py'."""

    def setUp(self):
        self.pepxml_file = '../data/A549/peptides.pep.xml'

    def test_estimate_mass_shift(self):
        average_mass_difference = estimate_mass_shift.fragment_ion_mass_shift(self.pepxml_file)
        print('Average ion mass difference: {0:3.1f} ppm'.format(average_mass_difference))


if __name__ == '__main__':
    unittest.main()

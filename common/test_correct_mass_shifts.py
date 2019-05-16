import unittest
import correct_mass_shifts
import estimate_mass_shift_from_mgf
from pyteomics import mgf


class CorrectMassShiftsTestCase(unittest.TestCase):
    """Tests for 'correct_mass_shifts.py'."""

    def setUp(self):
        self.mgf_file = '../data/workflow/NoCali/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15 (no smaple 249_15.6ppm)_NoCali_5_frac_FDR_0.1.mgf'
        self.corrected_file =  '../tmp/corrected.mgf'

    def test_estimate_mass_shift(self):
        precursor_mu, _, ion_mu, _ = estimate_mass_shift_from_mgf.estimate_mass_shift(self.mgf_file, plot=False)

        spectra = correct_mass_shifts.correct(self.mgf_file, precursor_mu, ion_mu)
        with open(self.corrected_file, 'w') as destination:
            mgf.write(spectra, destination)

        estimate_mass_shift_from_mgf.estimate_mass_shift(self.corrected_file, plot=True)


if __name__ == '__main__':
    unittest.main()

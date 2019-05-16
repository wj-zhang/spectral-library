import unittest
from spectrum_alignment import spectrum_alignment


class SpectrumAlignmentTestCase(unittest.TestCase):

    def test_spectrum_alignment(self):
        masses_a = [1, 2, 3]
        intensities_a = [6, 5, 4]

        masses_b = [1, 1.9, 2.1, 3.1]
        intensities_b = [6, 1, 5, 4]

        psm_a = {'m/z array': masses_a, 'intensity array': intensities_a}
        psm_b = {'m/z array': masses_b, 'intensity array': intensities_b}

        alignment = spectrum_alignment(psm_a, psm_b, ion_ppm=0.1)
        product = [intensities_a[i] * intensities_b[j] for i, j in alignment.items()]
        self.assertAlmostEqual(sum(product), 77, places=2)

        masses_b = [1, 1.9, 2.1, 3.3]
        intensities_b = [6, 1, 5, 4]

        psm_b = {'m/z array': masses_b, 'intensity array': intensities_b}

        alignment = spectrum_alignment(psm_a, psm_b, ion_ppm=0.06)
        product = [intensities_a[i] * intensities_b[j] for i, j in alignment.items()]
        self.assertAlmostEqual(sum(product), 61, places=2)


if __name__ == '__main__':
    unittest.main()

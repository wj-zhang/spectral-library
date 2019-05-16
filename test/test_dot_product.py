import unittest
from dot_product import dot_product


class DotProductTestCase(unittest.TestCase):

    def test_dot_product(self):
        masses_a = [1, 2, 3]
        intensities_a = [6, 5, 4]

        masses_b = [1, 1.9, 2.1, 3.1]
        intensities_b = [6, 1, 5, 4]

        psm_a = {'m/z array': masses_a, 'intensity array': intensities_a}
        psm_b = {'m/z array': masses_b, 'intensity array': intensities_b}

        product = dot_product(psm_a, psm_b)
        self.assertAlmostEqual(product, 77, places=2)

        product = dot_product(masses_a, intensities_a, masses_b, intensities_b)
        self.assertAlmostEqual(product, 77, places=2)

        masses_b = [1, 1.9, 2.1, 3.3]
        intensities_b = [6, 1, 5, 4]

        psm_b = {'m/z array': masses_b, 'intensity array': intensities_b}

        product = dot_product(psm_a, psm_b, error=0.2)
        self.assertAlmostEqual(product, 61, places=2)

        product = dot_product(masses_a, intensities_a, masses_b, intensities_b, error=0.2)
        self.assertAlmostEqual(product, 61, places=2)


if __name__ == '__main__':
    unittest.main()

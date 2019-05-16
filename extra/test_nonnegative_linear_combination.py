import unittest
import numpy as np
from nonnegative_linear_combination import nlc


class NonnegativeLinearCombinationTestCase(unittest.TestCase):
    """Tests for 'nonnegative_linear_combination.py'."""

    def setUp(self):
        self.C = np.array([[0.0372, 0.2869], [0.6861, 0.7071], [0.6233, 0.6245], [0.6344, 0.6170]])
        self.d = np.array([0.8587, 0.1781, 0.0747, 0.8405])
        self.options = {'show_progress': False}

    def test_nlc(self):
        w = nlc(self.C, self.d, options=self.options)
        self.assertAlmostEqual(w[1], 0.69, places=2)


if __name__ == '__main__':
    unittest.main()

import unittest
import numpy as np
from preprocessing import log_normalize


class PreprocessingTestCase(unittest.TestCase):
    def setUp(self):
        self.x = np.array([1, 2, 4])

    def test_log_normalize(self):
        x_log_normalize = log_normalize(self.x)
        self.assertAlmostEqual(x_log_normalize[1], 0.4472, places=2)


if __name__ == '__main__':
    unittest.main()

import unittest
import simulate_mixture_spectra
import numpy as np


class SimulateMixtureSpectraTestCase(unittest.TestCase):
    """Tests for 'simulate_mixture_spectra.py'."""

    def setUp(self):
        self.library_file = "../tmp/test_simulation/library.ms2"
        self.sim_our_library_file = "../tmp/test_simulation/simulation_library.ms2"
        self.sim_st_library_file = "../tmp/test_simulation/simulation_library_st.ms2"
        self.sim_query_file = "../tmp/test_simulation/simulation_query.mgf"
        self.id_file = "../tmp/test_simulation/simulation_id.csv"
        self.infor_file = "../tmp/test_simulation/simulation_infor.txt"
        self.random_seed = 0

        self.option = {'ratio': {1: 0.3, 2: 0.3, 3: 0.2, 4: 0.1, 5: 0.1},
                       'weights': {2: [[0.5, 0.5], [0.8, 0.2]],
                                   3: [[1 / 3, 1 / 3, 1 / 3], [0.5, 0.3, 0.2]],
                                   4: [[0.25, 0.25, 0.25, 0.25], [0.5, 0.2, 0.2, 0.1]],
                                   5: [[0.2, 0.2, 0.2, 0.2, 0.2], [0.25, 0.25, 0.2, 0.2, 0.1]]},
                       'mz range': 2
                       }

        self.mz = np.array([0.1, 0.11, 0.11, 0.3, 0.4])
        self.intensity = np.array([0.1, 0.1, 0.1, 0.1, 0.1])

    def test_simulate_mixture_spectra(self):
        infor = simulate_mixture_spectra.generate(self.library_file, self.sim_our_library_file, self.sim_st_library_file,
                                                  self.sim_query_file, self.id_file, self.infor_file, self.option,
                                                  random_state=self.random_seed, is_shuffle=False, write=False)
        print(infor)

    def test_combine_identical_peaks(self):
        mz, intensities = simulate_mixture_spectra.combine_identical_peaks(self.mz, self.intensity, 0.02)
        print(mz)
        print(intensities)


if __name__ == '__main__':
    unittest.main()

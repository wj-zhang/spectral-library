import unittest
import splib
import ms2
from splib_modifications import mod_env, convert_modification_forward
from pyteomics import mass


class SplibTestCase(unittest.TestCase):
    """Tests for 'splib.py'."""

    def setUp(self):
        self.file = '../data/splib/NIST_human_IT_2012-05-30_7AA/NIST_human_IT_2012-05-30_7AA.splib'
        self.file2 = '../tmp/nist.ms2'

    def test_read_splib(self):
        reader = splib.read(self.file)

        with open(self.file2, 'w') as destination:
            for psm in reader:
                ms2.write(destination, psm)

                modified_seq = convert_modification_forward(psm['params']['modified seq'])
                mz = psm['params']['scan'][2]
                z = psm['params']['charge'][0]
                calc_mz = mass.fast_mass2(modified_seq, charge=z, aa_mass=mod_env['aa_mass'])
                self.assertAlmostEqual(mz, calc_mz, places=2)


if __name__ == '__main__':
    unittest.main()

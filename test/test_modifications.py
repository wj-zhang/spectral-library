import unittest
import modifications
from pyteomics import mass


class ModificationsTestCase(unittest.TestCase):
    """Tests for 'modifications.py'."""

    def setUp(self):
        pass

    def test_convert_between_int_and_str(self):
        id = modifications.convert_between_int_and_str(24)
        self.assertEqual(id, 'ce')
        id = modifications.convert_between_int_and_str('ce')
        self.assertEqual(id, 24)

    def test_get_modification_by_record_id(self):
        mod = modifications.get_modification_by_record_id(24)
        self.assertEqual(mod['record_id'], 24)

    def test_unimod(self):
        db = modifications.unimod_db

        aa_comp = dict(mass.std_aa_comp)
        aa_comp['p'] = db.by_title('Phospho')['composition']
        m = mass.calculate_mass('PEpTIDE', aa_comp=aa_comp)
        self.assertAlmostEqual(m, 782.2735, places=2)

        aa_mass = dict(mass.std_aa_mass)
        aa_mass['p'] = mass.calculate_mass(aa_comp['p'])
        m = mass.fast_mass2('PEpTIDE', aa_mass=aa_mass)
        self.assertAlmostEqual(m, 782.2735, places=2)

        aa_comp = dict(mass.std_aa_comp)
        aa_comp['pT'] = db.by_title('Phospho')['composition'] + mass.std_aa_comp['T']
        m = mass.calculate_mass('PEpTIDE', aa_comp=aa_comp)
        self.assertAlmostEqual(m, 782.2735, places=2)

        aa_mass = dict(mass.std_aa_mass)
        aa_comp['pT'] = db.by_title('Phospho')['composition'] + mass.std_aa_comp['T']
        aa_mass['pT'] = mass.calculate_mass(aa_comp['pT'])
        m = mass.fast_mass2('PEpTIDE', aa_mass=aa_mass)
        self.assertAlmostEqual(m, 782.2735, places=2)

    def test_modifications_env(self):
        mod_digits = [4, 35, 7]
        mod_env = modifications.modifications_env(mod_digits)
        self.assertEqual(mod_env['symbols'], ['e', 'df', 'h'])


if __name__ == '__main__':
    unittest.main()

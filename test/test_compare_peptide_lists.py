import unittest
from compare_peptide_lists import compare_peptide_lists


class ComparePeptideListsCase(unittest.TestCase):
    def setUp(self):
        self.file_1 = '../tmp/library_id_1.csv'
        self.file_2 = '../tmp/library_id_2.csv'
        self.out_1 = '../tmp/a_only.csv'
        self.out_2 = '../tmp/b_only.csv'
        self.out_3 = '../tmp/ab_common.csv'

    def test_compare_peptide_lists(self):
        a_only, b_only, ab_common = compare_peptide_lists(self.file_1, self.file_2)
        print('Peptide number: {} (a_only), {} (b_only), {} (ab_common)'.format(len(a_only.index), len(b_only.index),
                                                                                len(ab_common.index)))

        a_only.to_csv(self.out_1, sep=',', index=False)
        b_only.to_csv(self.out_2, sep=',', index=False)
        ab_common.to_csv(self.out_3, sep=',', index=False)


if __name__ == '__main__':
    unittest.main()

import unittest
import plot_venn_from_csv


class PlotVennFromCsvTestCase(unittest.TestCase):
    """Tests for 'plot_venn_from_csv.py'."""

    def setUp(self):
        self.files = {
            'Ours': '../tmp/CharmeRT/real/Ours/HeLa_1ug_isow2_1hgradient_datasetA_rep1_id.csv',
            'MSPLIT': '../tmp/CharmeRT/real/MSPLIT/HeLa_1ug_isow2_1hgradient_datasetA_rep1_id.csv',
            'SpectraST': '../tmp/CharmeRT/real/SpectraST/HeLa_1ug_isow2_1hgradient_datasetA_rep1_id.csv',
        }
        self.column_name = [
            # 'Index',
            # 'Protein',
            'Peptide',
        ]

    def test_plot_venn(self):
        plot_venn_from_csv.plot_venn(self.files, self.column_name)


if __name__ == '__main__':
    unittest.main()

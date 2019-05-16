import unittest
import run_pipeline_MSPLIT


class RunPipelineMSPLITTestCase(unittest.TestCase):
    """Tests for 'run_pipeline_MSPLIT.py'."""

    def setUp(self):
        self.decoy_library_file = '../tmp/CharmeRT/SpectraST/HeLa_1ug_isow2_3hgradient_datasetA_chimera_Cons_Q_Decoy.sptxt'
        self.ms_file = '../data/CharmeRT/search/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera.mgf'
        self.raw_id_file = '../tmp/CharmeRT/MSPLIT/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera.txt'
        self.svm_tmp_dir = '../tmp/CharmeRT/MSPLIT/svm_tmp'
        self.fdr_id_file = '../tmp/CharmeRT/MSPLIT/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera_fdr.txt'
        self.id_file = '../tmp/CharmeRT/MSPLIT/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera.csv'

    def test_run_pipeline(self):
        print('library search ...')
        run_pipeline_MSPLIT.run_library_search(self.decoy_library_file, self.ms_file, self.raw_id_file)

        print('perform fdr ...')
        run_pipeline_MSPLIT.run_fdr(self.raw_id_file, self.fdr_id_file, self.svm_tmp_dir, 0.01)

        print('format id file to csv file ...')
        run_pipeline_MSPLIT.run_format_to_csv_file(self.fdr_id_file, self.id_file)


if __name__ == '__main__':
    unittest.main()

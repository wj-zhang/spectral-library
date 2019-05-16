import unittest
import run_pipeline_MSPLIT


class RunPipelineMSPLITTestCase(unittest.TestCase):
    """Tests for 'run_pipeline_MSPLIT.py'."""

    def setUp(self):
        self.library_file = '../tmp/CharmeRT/simulation/SpectraST/simulation_library.sptxt'
        self.ms_file = '../tmp/CharmeRT/simulation/data/simulation_query.mgf'

        self.svm_tmp_dir = '../tmp/CharmeRT/simulation/MSPLIT/svm_tmp'
        self.raw_id_file = '../tmp/CharmeRT/simulation/MSPLIT/simulation_query.txt'
        self.fdr_id_file = '../tmp/CharmeRT/simulation/MSPLIT/simulation_query_fdr.txt'
        self.id_file = '../tmp/CharmeRT/simulation/MSPLIT/simulation_query.csv'

    def test_run_pipeline(self):
        print('library search ...')
        run_pipeline_MSPLIT.run_library_search(self.library_file, self.ms_file, self.raw_id_file)

        print('perform fdr ...')
        run_pipeline_MSPLIT.run_fdr_by_threshold(self.raw_id_file, self.fdr_id_file, self.svm_tmp_dir, 0, 0)

        print('format id file to csv file ...')
        run_pipeline_MSPLIT.run_format_to_csv_file(self.fdr_id_file, self.id_file)


if __name__ == '__main__':
    unittest.main()

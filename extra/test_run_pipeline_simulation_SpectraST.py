import unittest
import run_pipeline_SpectraST


class RunPipelineSpectraSTTestCase(unittest.TestCase):
    """Tests for 'run_pipeline_SpectraST.py'."""

    def setUp(self):
        self.ms2_library_file = '../tmp/CharmeRT/simulation/data/simulation_library_st.ms2'
        self.ms_file = '../tmp/CharmeRT/simulation/data/simulation_query.mgf'

        self.library_file = '../tmp/CharmeRT/simulation/SpectraST/simulation_library'
        self.raw_id_file = '../tmp/CharmeRT/simulation/SpectraST/simulation_query.xls'

        self.id_file = '../tmp/CharmeRT/simulation/SpectraST/simulation_query.csv'

    def test_run_pipeline(self):
        print('build library from ms2 file ...')
        run_pipeline_SpectraST.run_build_library(self.ms2_library_file, self.library_file)

        print('library search ...')
        run_pipeline_SpectraST.run_library_search(self.library_file, self.ms_file, self.raw_id_file)

        print('fdr ...')
        run_pipeline_SpectraST.run_fdr(self.raw_id_file, self.id_file, fdr=0.01)


if __name__ == '__main__':
    unittest.main()

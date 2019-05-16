import unittest
import run_pipeline


class RunPipelineTestCase(unittest.TestCase):
    """Tests for 'run_pipeline.py'."""

    def setUp(self):
        self.mgf_files = [
            '../data/workflow/34 200 pg A549 target 20 (15.6ppm)_20190107_1_frac_cali_FDR_0.1.mgf',
            '../data/workflow/34 200 pg A549 target 20 (15.6ppm)_20190107_2_frac_cali_FDR_0.1.mgf',
            '../data/workflow/34 200 pg A549 target 20 (15.6ppm)_20190107_3_frac_cali_FDR_0.1.mgf',
            '../data/workflow/200 ng hela_20181229_26_1_frac_cali_FDR_0.1.mgf',
            '../data/workflow/200 ng hela_20181229_27_2_frac_cali_FDR_0.1.mgf',
            '../data/workflow/200 ng hela_20181229_207_5_frac_cali_FDR_0.1.mgf',
            '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15 (no smaple 249_15.6ppm)_index_1_frac_FDR_0.1.mgf',
            '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15 (no smaple 249_15.6ppm)_index_2_frac_FDR_0.1.mgf',
            '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15 (no smaple 249_15.6ppm)_index_3_frac_FDR_0.1.mgf',
            '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15 (no smaple 249_15.6ppm)_index_4_frac_FDR_0.1.mgf',
            '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15 (no smaple 249_15.6ppm)_index_5_frac_FDR_0.1.mgf'
        ]

        self.library_file = '../tmp/workflow/peaks_library.ms2'
        self.decoy_file = '../tmp/workflow/decoy_library.ms2'

        self.ms_file = '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15_Slot1-19_1_249.d (15.6ppm)_1_frac_FDR_5_index.mgf'
        self.id_file = '../tmp/workflow/library_id_decoy.csv'

    def test_run_pipeline(self):
        print('build library from mgf ...')
        run_pipeline.run_build_library_from_mgf(self.mgf_files, self.library_file)
        print('remove conflicting identifications ...')
        run_pipeline.run_remove_conflicting_identifications(self.library_file)
        print('generate decoy library ...')
        run_pipeline.run_generate_decoy_library(self.library_file, self.decoy_file)
        print('search library from PEAKS ...')
        run_pipeline.run_search_library_from_PEAKS(self.ms_file, self.decoy_file, self.id_file)
        print('predict retention time ...')
        run_pipeline.run_predict_retention_time(self.id_file)


if __name__ == '__main__':
    unittest.main()

import unittest
import run_pipeline_two_rounds


class RunPipelineTwoRoundsTestCase(unittest.TestCase):
    """Tests for 'run_pipeline_two_rounds.py'."""

    def setUp(self):
        self.mgf_files = [
            '../data/CharmeRT/library/HeLa_1ug_isow2_3h_all_frac_No_Cali_FDR_0.1%.mgf'
        ]
        self.corrected_file = '../tmp/CharmeRT/mgf/corrected_mass.mgf'

        self.library_file = '../tmp/CharmeRT/mgf/peaks_library.ms2'
        self.decoy_file = '../tmp/CharmeRT/mgf/decoy_library.ms2'

        self.library_target_file = '../tmp/PRG2017/corrected_target.ms2'
        self.library_decoy_file = '../tmp/PRG2017/corrected_decoy.ms2'

        self.ms_file = '../data/PRG2017/10222015_ABRF3_1000ng_68X_all_frac_No_Cali_FDR_1%.mgf'

        self.id_target_file = '../tmp/PRG2017/library_id_target.csv'
        self.id_decoy_file = '../tmp/PRG2017/library_id_decoy.csv'
        self.id_file = '../tmp/PRG2017/library_id.csv'

    def test_run_pipeline(self):
        print('correct mass shifts ...')
        run_pipeline_two_rounds.run_correct_mass_shifts(self.mgf_files, self.corrected_file, mode='library')
        print('build library from mgf ...')
        run_pipeline_two_rounds.run_build_library_from_mgf(self.corrected_file, self.library_file)
        print('remove conflicting identifications ...')
        run_pipeline_two_rounds.run_remove_conflicting_identifications(self.library_file)
        # print('generate decoy library ...')
        # run_pipeline_two_rounds.run_generate_decoy_library(self.library_file, self.decoy_file)
        # print('seperate target and decoy ...')
        # run_pipeline_two_rounds.run_separate_target_decoy(self.decoy_file, self.library_target_file, self.library_decoy_file)
        #

        # print('pre-library search ...')
        # spectra = run_pipeline_two_rounds.run_pre_library_search(self.ms_file, self.decoy_file)
        # print('correct mass shifts ...')
        # run_pipeline_two_rounds.run_correct_mass_shifts((self.ms_file, spectra), self.corrected_file, mode='query')
        # print('combine clone spectra ...')
        # run_pipeline_two_rounds.run_combine_clone_mgf(self.corrected_file)
        # print('main library search (target) ...')
        # run_pipeline_two_rounds.run_main_library_search(self.corrected_file, self.library_target_file, self.id_target_file)
        # print('main library search (decoy) ...')
        # run_pipeline_two_rounds.run_main_library_search(self.corrected_file, self.library_decoy_file, self.id_decoy_file)
        # print('merge target and decoy ids ...')
        # run_pipeline_two_rounds.run_merge_target_decoy_ids(self.id_target_file, self.id_decoy_file, self.id_file)
        # print('predict retention time ...')
        # run_pipeline_two_rounds.run_predict_retention_time(self.id_file)


if __name__ == '__main__':
    unittest.main()

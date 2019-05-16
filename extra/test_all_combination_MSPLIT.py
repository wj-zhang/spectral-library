import unittest
import run_pipeline_MSPLIT
import run_pipeline_two_rounds


class RunPipelineTwoRoundsTestCase(unittest.TestCase):
    """Tests for 'run_pipeline_two_rounds.py'."""

    def setUp(self):
        # self.pepxml_file = '../data/CharmeRT/library/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1.pep.xml'
        # self.library_corrected_dir = '../tmp/CharmeRT/real/ours'
        # self.library_corrected_suffix = '_corrected'

        # self.library_file = '../tmp/CharmeRT/real/ours/peaks_library.ms2'
        # self.decoy_file = '../tmp/CharmeRT/real/ours/decoy_library.ms2'
        self.library_all_file = '../tmp/CharmeRT/real/Ours/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1_all.ms2'
        self.library_target_file = '../tmp/CharmeRT/real/Ours/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1_target.ms2'
        self.library_decoy_file = '../tmp/CharmeRT/real/Ours/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1_decoy.ms2'

        # self.single_ms_file = '../data/CharmeRT/search/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera.mgf'
        # self.clone_ms_file = '../data/CharmeRT/search/HeLa_1ug_isow2_1hgradient_datasetA_rep1_chimera.mgf'

        # self.single_corrected_ms_file = '../tmp/CharmeRT/real/ours/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera_corrected.mgf'
        self.combined_ms_file = '../tmp/CharmeRT/real/Ours_before_bb5a0101/HeLa_1ug_isow8_1hgradient_datasetE_rep1_non_chimera_combined.mgf'
        self.ms_file = '../data/CharmeRT/search/HeLa_1ug_isow8_1hgradient_datasetE_rep1_non_chimera.mgf'
        self.raw_id_target_file = '../tmp/joint_optimization_score/M_library_id_target.txt'
        self.raw_id_decoy_file = '../tmp/joint_optimization_score/M_library_id_decoy.txt'
        self.raw_id_fdr_target_file = '../tmp/joint_optimization_score/M_library_id_target_fdr.txt'
        self.raw_id_fdr_decoy_file = '../tmp/joint_optimization_score/M_library_id_decoy_fdr.txt'
        self.id_target_file = '../tmp/joint_optimization_score/M_library_id_target.csv'
        self.id_decoy_file = '../tmp/joint_optimization_score/M_library_id_decoy.csv'

        # self.id_target_file = '../tmp/joint_optimization_score/library_id_target_1_mean.csv'
        # self.id_decoy_file = '../tmp/joint_optimization_score/library_id_decoy_1_mean.csv'

        self.id_file = '../tmp/joint_optimization_score/null.csv'

        self.fdr = 1/100

        self.id_file_demix_max = '../tmp/joint_optimization_score/library_id_demix_max.csv'
        self.id_file_demix_multi = '../tmp/joint_optimization_score/library_id_demix_multi.csv'
        self.id_file_demix_my = '../tmp/joint_optimization_score/library_id_demix_my.csv'
        self.svm_tmp = '../tmp/joint_optimization_score/svm_tmp'
        # self.corrected_pepxml_file = '../tmp/CharmeRT/real/ours/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1_corrected.pep.xml'

    def test_run_pipeline(self):
        # print('correct mass shifts of library ...')
        # corrected_pepxml_file = run_pipeline_two_rounds.run_correct_mass_shifts_from_pepxml(self.pepxml_file,
        #                                                                                     self.library_corrected_dir,
        #                                                                                     self.library_corrected_suffix)
        # print('build library from pepxml ...')
        # run_pipeline_two_rounds.run_build_library_from_pepxml(corrected_pepxml_file, self.library_file)
        # print('remove conflicting identifications ...')
        # run_pipeline_two_rounds.run_remove_conflicting_identifications(self.library_file)
        # print('generate decoy library ...')
        # run_pipeline_two_rounds.run_generate_decoy_library(self.library_file, self.decoy_file)
        # print('seperate target and decoy ...')
        # run_pipeline_two_rounds.run_separate_target_decoy(self.decoy_file, self.library_target_file, self.library_decoy_file)
        #
        # print('pre-library search ...')
        # spectra = run_pipeline_two_rounds.run_pre_library_search(self.single_ms_file, self.decoy_file)
        # print('correct mass shifts of query ...')
        # run_pipeline_two_rounds.run_correct_mass_shifts((self.single_ms_file, spectra), self.single_corrected_ms_file, mode='query')
        #
        # print('combine clone spectra ...')
        # run_pipeline_two_rounds.run_combine_clone_from_two_mgfs(self.single_corrected_ms_file, self.clone_ms_file, self.combined_ms_file)

        # print('library search demix my ...')
        # run_pipeline_two_rounds.run_library_search_dimix_my(self.combined_ms_file, self.library_all_file, self.id_file_demix_my)

        # print('main library search (target) ...')
        # run_pipeline_MSPLIT.run_library_search(self.library_target_file, self.ms_file, self.raw_id_target_file)

        #
        # print('main library search (decoy) ...')
        # run_pipeline_MSPLIT.run_library_search(self.library_decoy_file, self.ms_file, self.raw_id_decoy_file)

        # run_pipeline_MSPLIT.run_fdr(self.raw_id_target_file, self.raw_id_fdr_target_file, self.svm_tmp, 0)
        # run_pipeline_MSPLIT.run_fdr_by_threshold(self.raw_id_decoy_file, self.raw_id_fdr_decoy_file, self.svm_tmp, -100, -100)

        # print('formate (target) ...')
        run_pipeline_MSPLIT.run_convert_raw_id_to_csv_file(self.raw_id_fdr_target_file, self.id_target_file, False)
        run_pipeline_MSPLIT.run_convert_raw_id_to_csv_file(self.raw_id_fdr_decoy_file, self.id_decoy_file, True)


        # print('main library search (decoy) ...')
        # run_pipeline_MSPLIT.run_library_search(self.library_decoy_file, self.ms_file, self.raw_id_decoy_file)
        #
        # print('main library search (decoy) ...')
        # run_pipeline_two_rounds.run_main_library_search(self.combined_ms_file, self.library_decoy_file, self.id_decoy_file)
        print('merge target and decoy ids ...')
        run_pipeline_two_rounds.run_merge_target_decoy_ids(self.id_target_file, self.id_decoy_file, self.id_file, self.fdr)


if __name__ == '__main__':
    unittest.main()

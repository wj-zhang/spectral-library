import unittest
import run_pipeline_two_rounds


class RunPipelineTwoRoundsTestCase(unittest.TestCase):
    """Tests for 'run_pipeline_two_rounds.py'."""

    def setUp(self):
        # self.pepxml_file = '../data/CharmeRT/library/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1.pep.xml'
        # self.library_corrected_dir = '../tmp/CharmeRT/real/ours'
        # self.library_corrected_suffix = '_corrected'

        self.library_file = '../tmp/CharmeRT/real/ours/peaks_library.ms2'
        self.decoy_file = '../tmp/CharmeRT/real/ours/decoy_library.ms2'
        self.library_all_file = '../tmp/CharmeRT/real/Ours/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1_all.ms2'
        # self.library_target_file = '../tmp/CharmeRT/real/Ours/bk/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1_target.ms2'
        self.library_target_file = '../tmp/CharmeRT/simulation/data/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1_simulation_library.ms2'
        self.library_decoy_file = '../tmp/CharmeRT/real/Ours/bk/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1_decoy_1.ms2'

        # self.single_ms_file = '../data/CharmeRT/search/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera.mgf'
        # self.clone_ms_file = '../data/CharmeRT/search/HeLa_1ug_isow2_1hgradient_datasetA_rep1_chimera.mgf'

        # self.single_corrected_ms_file = '../tmp/CharmeRT/real/ours/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera_corrected.mgf'
        # self.combined_ms_file = '../tmp/CharmeRT/real/Ours_before_bb5a0101/HeLa_1ug_isow8_1hgradient_datasetE_rep1_non_chimera_combined.mgf'
        self.combined_ms_file = '../tmp/CharmeRT/simulation/data/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1_simulation_query.mgf'
        self.id_target_file = '../tmp/joint_optimization_score/library_id_target_new_2_28_max.csv'
        self.id_decoy_file = '../tmp/joint_optimization_score/library_id_decoy_new_2_28_max.csv'

        # self.id_target_file = '../tmp/joint_optimization_score/library_id_target_1_mean.csv'
        # self.id_decoy_file = '../tmp/joint_optimization_score/library_id_decoy_1_mean.csv'

        self.id_file = '../tmp/joint_optimization_score/null.csv'

        self.fdr = 0.1/100
        self.ratio = 1

        self.id_file_demix_max = '../tmp/joint_optimization_score/library_id_demix_max.csv'
        self.id_file_demix_multi = '../tmp/joint_optimization_score/library_id_demix_multi.csv'
        self.id_file_demix_my = '../tmp/joint_optimization_score/library_id_demix_my.csv'
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
        # run_pipeline_two_rounds.run_generate_decoy_library(self.library_target_file, self.library_decoy_file, ratio=self.ratio)

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
        #
        print('main library search (target) ...')
        run_pipeline_two_rounds.run_main_library_search(self.combined_ms_file, self.library_target_file, self.id_target_file)

        # print('main library search (decoy) ...')
        # run_pipeline_two_rounds.run_main_library_search(self.combined_ms_file, self.library_decoy_file, self.id_decoy_file)

        print('merge target and decoy ids ...')
        run_pipeline_two_rounds.run_merge_target_decoy_ids(self.id_target_file, self.id_decoy_file, self.id_file, self.fdr)


if __name__ == '__main__':
    unittest.main()

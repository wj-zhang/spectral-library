import unittest
import run_pipeline_two_rounds


class RunPipelineTwoRoundsTestCase(unittest.TestCase):
    """Tests for 'run_pipeline_two_rounds.py'."""

    def setUp(self):
        self.pepxml_file = '../data/CharmeRT/library/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1.pep.xml'
        self.library_corrected_dir = '../tmp/CharmeRT/'
        self.library_corrected_suffix = '_corrected'

        self.library_file = '../tmp/CharmeRT/peaks_library.ms2'
        self.decoy_file = '../tmp/CharmeRT/decoy_library.ms2'
        self.library_target_file = '../tmp/CharmeRT/corrected_target.ms2'
        self.library_decoy_file = '../tmp/CharmeRT/corrected_decoy.ms2'

        self.single_ms_file = '../data/CharmeRT/search/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera.mgf'
        self.clone_ms_file = '../data/CharmeRT/search/HeLa_1ug_isow2_1hgradient_datasetA_rep1_chimera.mgf'

        self.single_corrected_ms_file = '../tmp/CharmeRT/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera_corrected.mgf'
        self.combined_ms_file = '../tmp/CharmeRT/HeLa_1ug_isow2_1hgradient_datasetA_rep1_combined.mgf'

        self.id_target_file = '../tmp/CharmeRT/library_id_target.csv'
        self.id_decoy_file = '../tmp/CharmeRT/library_id_decoy.csv'
        self.id_file = '../tmp/CharmeRT/library_id.csv'
        self.corrected_pepxml_file = '../tmp/CharmeRT/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1_corrected.pep.xml'

    def test_run_pipeline(self):
        # print('correct mass shifts of library ...')
        # corrected_pepxml_file = run_pipeline_two_rounds.run_correct_mass_shifts_from_pepxml(self.pepxml_file,
        #                                                                                     self.library_corrected_dir,
        #                                                                                     self.library_corrected_suffix)
        print('build library from pepxml ...')
        run_pipeline_two_rounds.run_build_library_from_pepxml(self.corrected_pepxml_file, self.library_file)
        print('remove conflicting identifications ...')
        run_pipeline_two_rounds.run_remove_conflicting_identifications(self.library_file)
        print('generate decoy library ...')
        run_pipeline_two_rounds.run_generate_decoy_library(self.library_file, self.decoy_file)
        print('seperate target and decoy ...')
        run_pipeline_two_rounds.run_separate_target_decoy(self.decoy_file, self.library_target_file, self.library_decoy_file)

        print('pre-library search ...')
        spectra = run_pipeline_two_rounds.run_pre_library_search(self.single_ms_file, self.decoy_file)
        print('correct mass shifts of query ...')
        run_pipeline_two_rounds.run_correct_mass_shifts((self.single_ms_file, spectra), self.single_corrected_ms_file, mode='query')

        print('combine clone spectra ...')
        run_pipeline_two_rounds.run_combine_clone_from_two_mgfs(self.single_corrected_ms_file, self.clone_ms_file, self.combined_ms_file)
        print('main library search (target) ...')
        run_pipeline_two_rounds.run_main_library_search(self.combined_ms_file, self.library_target_file, self.id_target_file)
        print('main library search (decoy) ...')
        run_pipeline_two_rounds.run_main_library_search(self.combined_ms_file, self.library_decoy_file, self.id_decoy_file)
        print('merge target and decoy ids ...')
        run_pipeline_two_rounds.run_merge_target_decoy_ids(self.id_target_file, self.id_decoy_file, self.id_file)


if __name__ == '__main__':
    unittest.main()

import unittest
import run_pipeline_two_rounds


class RunPipelineTwoRoundsTestCase(unittest.TestCase):
    """Tests for 'run_pipeline_two_rounds.py'."""

    def setUp(self):

        self.pepxml_files = '../data/CharmeRT/library/HeLa_1ug_isow2_3hgradient_datasetA_chimera.pep.xml'
        self.library_file = '../tmp/CharmeRT/peaks_library.ms2'
        self.decoy_file = '../tmp/CharmeRT/decoy_library.ms2'
        self.library_target_file = '../tmp/CharmeRT/corrected_target.ms2'
        self.library_decoy_file = '../tmp/CharmeRT/corrected_decoy.ms2'

        self.single_ms_file = '../data/CharmeRT/search/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera.mgf'
        self.clone_ms_file = '../data/CharmeRT/search/HeLa_1ug_isow2_1hgradient_datasetA_rep1_chimera.mgf'
        self.ms_file = '../tmp/CharmeRT/HeLa_1ug_isow2_1hgradient_datasetA_rep1_combined.mgf'

        self.id_target_file = '../tmp/CharmeRT/library_id_target.csv'
        self.id_decoy_file = '../tmp/CharmeRT/library_id_decoy.csv'
        self.id_file = '../tmp/CharmeRT/library_id.csv'

    def test_run_pipeline(self):
        print('build library from pepxml ...')
        run_pipeline_two_rounds.run_build_library_from_pepxml(self.pepxml_files, self.library_file)
        print('remove conflicting identifications ...')
        run_pipeline_two_rounds.run_remove_conflicting_identifications(self.library_file)
        print('generate decoy library ...')
        run_pipeline_two_rounds.run_generate_decoy_library(self.library_file, self.decoy_file)
        print('seperate target and decoy ...')
        run_pipeline_two_rounds.run_separate_target_decoy(self.decoy_file, self.library_target_file, self.library_decoy_file)
        print('combine clone spectra ...')
        run_pipeline_two_rounds.run_combine_clone_from_two_mgfs(self.single_ms_file, self.clone_ms_file, self.ms_file)
        print('main library search (target) ...')
        run_pipeline_two_rounds.run_main_library_search(self.ms_file, self.library_target_file, self.id_target_file)
        print('main library search (decoy) ...')
        run_pipeline_two_rounds.run_main_library_search(self.ms_file, self.library_decoy_file, self.id_decoy_file)
        print('merge target and decoy ids ...')
        run_pipeline_two_rounds.run_merge_target_decoy_ids(self.id_target_file, self.id_decoy_file, self.id_file)


if __name__ == '__main__':
    unittest.main()

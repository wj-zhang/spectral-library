import unittest
import run_pipeline_SpectraST


class RunPipelineSpectraSTTestCase(unittest.TestCase):
    """Tests for 'run_pipeline_SpectraST.py'."""

    def setUp(self):
        self.PEAKS_pepxml_file = '../data/CharmeRT/library/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1.pep.xml'
        self.pepxml_file = '../data/CharmeRT/library/HeLa_1ug_isow2_3hgradient_datasetA_chimera.st.pep.xml'
        self.raw_library_file = '../tmp/CharmeRT/SpectraST/HeLa_1ug_isow2_3hgradient_datasetA_chimera'
        self.consensus_library_file = '../tmp/CharmeRT/SpectraST/HeLa_1ug_isow2_3hgradient_datasetA_chimera_cons'
        self.quality_library_file = '../tmp/CharmeRT/SpectraST/HeLa_1ug_isow2_3hgradient_datasetA_chimera_Cons_Q'
        self.decoy_file = '../tmp/CharmeRT/SpectraST/HeLa_1ug_isow2_3hgradient_datasetA_chimera_Cons_Q_Decoy'

        self.ms_file = '../data/CharmeRT/search/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera.mgf'
        self.raw_id_file = '../tmp/CharmeRT/SpectraST/HeLa_1ug_isow2_1hgradient_datasetA_rep1.xls'
        self.id_file = '../tmp/CharmeRT/SpectraST/HeLa_1ug_isow2_1hgradient_datasetA_rep1_non_chimera.csv'

    def test_run_pipeline(self):
        print('format pepxml file for SpectraST ...')
        run_pipeline_SpectraST.run_format_pepxml(self.PEAKS_pepxml_file, self.pepxml_file)
        print('build raw library from pepxml file ...')
        run_pipeline_SpectraST.run_build_library(self.pepxml_file, self.raw_library_file)

        print('build consensus library ...')
        run_pipeline_SpectraST.run_build_consensus_library(self.raw_library_file, self.consensus_library_file)
        print('perform quality control of library ...')
        run_pipeline_SpectraST.run_quality_control_library(self.consensus_library_file, self.quality_library_file)
        print('generate decoy library ...')
        run_pipeline_SpectraST.run_generate_decoy_library(self.quality_library_file, self.decoy_file)
        print('library search ...')
        run_pipeline_SpectraST.run_library_search(self.decoy_file, self.ms_file, self.raw_id_file)
        print('fdr ...')
        run_pipeline_SpectraST.run_fdr(self.raw_id_file, self.id_file, fdr=0.01)


if __name__ == '__main__':
    unittest.main()

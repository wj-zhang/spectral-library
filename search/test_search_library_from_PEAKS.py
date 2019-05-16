import unittest
import search_library_from_PEAKS
from peaks_modifications import convert_modification_backward


class SearchLibraryFromPeaksTestCase(unittest.TestCase):
    """Tests for 'search_library_from_PEAKS.py'."""

    def setUp(self):
        self.ms_file = '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15_Slot1-19_1_249.d (15.6ppm)_1_frac_FDR_5_index.mgf'
        # self.ms_file = '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15_Slot1-19_1_249.d (15.6ppm).mgf'
        # self.ms_file = '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15_Slot1-19_1_249.d (15.6ppm)_1_frac_FDR_5_index.mgf'
        self.library_file = '../tmp/workflow/decoy_library_test_11_files.ms2'
        self.id_file = '../tmp/workflow/library_id_test_11_files_decoy.csv'

    def test_search_library_from_PEAKS_search(self):
        psms = search_library_from_PEAKS.search(self.ms_file, self.library_file, fdr=0.0)
        psms['Peptide'] = psms['Peptide'].apply(convert_modification_backward)
        psms.to_csv(self.id_file, header=True, index=False)


if __name__ == '__main__':
    unittest.main()

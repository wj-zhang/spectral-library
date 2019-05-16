import unittest
import search_library_from_mgf
from peaks_modifications import convert_modification_backward


class SearchLibraryFromMgfTestCase(unittest.TestCase):
    """Tests for 'search_library_from_mgf.py'."""

    def setUp(self):
        self.ms_file = '../tmp/associated.mgf'
        self.library_file = '../tmp/workflow/decoy_library_5.ms2'
        self.id_file = '../tmp/workflow/library_id_test_11_files_decoy.csv'

    def test_search_library_from_mgf(self):
        psms = search_library_from_mgf.search(self.ms_file, self.library_file, fdr=0.01)
        psms['Peptide'] = psms['Peptide'].apply(convert_modification_backward)
        psms.to_csv(self.id_file, header=True, index=False)


if __name__ == '__main__':
    unittest.main()

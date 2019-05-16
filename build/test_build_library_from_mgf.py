import unittest
import build_library_from_mgf
import ms2


class BuildLibraryFromMgfTestCase(unittest.TestCase):
    """Tests for 'build_library_from_mgf.py'."""

    def setUp(self):
        self.mgf_files = '../data/PEAKS_peptides.mgf'
        self.mgf_files = [
            # '../data/workflow/200 ng hela_20181229_26_1_frac_FDR_0.1.mgf',
            # '../data/workflow/200 ng hela_20181229_27_2_frac_FDR_0.1.mgf',
            '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15 (no smaple 249_15.6ppm)_index_1_frac_FDR_0.1.mgf',
            '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15 (no smaple 249_15.6ppm)_index_2_frac_FDR_0.1.mgf',
            '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15 (no smaple 249_15.6ppm)_index_3_frac_FDR_0.1.mgf',
            '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15 (no smaple 249_15.6ppm)_index_4_frac_FDR_0.1.mgf',
            '../data/workflow/18. cur vs pre2_K0.75-1.2_200ms_1700V_P6_Thre0.5_Tar15 (no smaple 249_15.6ppm)_index_5_frac_FDR_0.1.mgf'
            ]
        self.library_file = '../tmp/workflow/peaks_library.ms2'

    def test_build_library_from_mgf(self):
        header = build_library_from_mgf.create_header(self.mgf_files)
        creator = build_library_from_mgf.create(self.mgf_files)

        with open(self.library_file, 'w') as destination:
            ms2.write_header(destination, header)
            for psm in creator:
                ms2.write(destination, psm)


if __name__ == '__main__':
    unittest.main()

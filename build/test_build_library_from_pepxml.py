import unittest
import build_library_from_pepxml
import ms2


class BuildLibraryFromPepxmlTestCase(unittest.TestCase):
    """Tests for 'build_library_from_pepxml.py'."""

    def setUp(self):
        self.pepxml_files = '../data/workflow/pepxml/peptides.pep.xml'
        # self.pepxml_files = ['../data/PEAKS_peptides.pep.xml', '../data/Cell_peptides.pep.xml']
        self.library_file = '../tmp/workflow/pepxml/peaks_library.ms2'

    def test_build_library_from_PEAKS(self):
        header = build_library_from_pepxml.create_header(self.pepxml_files)
        creator = build_library_from_pepxml.create(self.pepxml_files)

        with open(self.library_file, 'w') as destination:
            ms2.write_header(destination, header)
            for psm in creator:
                ms2.write(destination, psm)


if __name__ == '__main__':
    unittest.main()

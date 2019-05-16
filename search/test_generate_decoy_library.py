import unittest
import generate_decoy_library


class GenerateDecoyLibraryTestCase(unittest.TestCase):
    """Tests for 'generate_decoy_library.py'."""

    def setUp(self):
        self.target_file = '../tmp/workflow/peaks_library_test_11_files.ms2'
        self.decoy_file = '../tmp/workflow/decoy_library_test_11_files.ms2'

    def test_generate_decoy_library_generate(self):
        generate_decoy_library.generate(self.target_file, self.decoy_file, combined=True)


if __name__ == '__main__':
    unittest.main()

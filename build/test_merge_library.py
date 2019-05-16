import unittest
import merge_library
import ms2


class MergeLibraryTestCase(unittest.TestCase):
    """Tests for 'merge_library.py'."""

    def setUp(self):
        self.library_files = [
            '../tmp/peaks_library.ms2',
            '../tmp/workflow/peaks_library.ms2'
            ]
        self.library_file = '../tmp/merged_peaks_library.ms2'

    def test_merge_library(self):
        merger = merge_library.merge(self.library_files)

        with open(self.library_file, 'w') as destination:
            for psm in merger:
                ms2.write(destination, psm)


if __name__ == '__main__':
    unittest.main()

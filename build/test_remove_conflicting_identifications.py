import unittest
import remove_conflicting_identifications
import ms2


class RemoveConflictingIdentificationsTestCase(unittest.TestCase):
    """Tests for 'remove_conflicting_identifications.py'."""

    def setUp(self):
        self.library_file = '../tmp/workflow/peaks_library.ms2'

    def test_remove_conflicting_identifications(self):
        header = ms2.read_header(self.library_file)
        spectrum_list = remove_conflicting_identifications.remove(self.library_file)

        with open(self.library_file, 'w') as destination:
            ms2.write_header(destination, header)
            for psm in spectrum_list:
                ms2.write(destination, psm)


if __name__ == '__main__':
    unittest.main()

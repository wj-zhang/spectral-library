import unittest
import build_library_from_splib
import ms2


class BuildLibraryFromSplibTestCase(unittest.TestCase):
    """Tests for 'library_from_splib.py'."""

    def setUp(self):
        self.splib_file = '../data/splib/A01.splib'
        self.library_file = '../tmp/splib_library.ms2'

    def test_build_library_from_PEAKS(self):
        header = build_library_from_splib.build_header(self.splib_file)
        builder = build_library_from_splib.build(self.splib_file)

        with open(self.library_file, 'w') as destination:
            ms2.write_header(destination, header)
            for psm in builder:
                ms2.write(destination, psm)


if __name__ == '__main__':
    unittest.main()

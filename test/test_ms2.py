import unittest
import ms2


class Ms2TestCase(unittest.TestCase):
    """Tests for 'ms2.py'."""

    def setUp(self):
        self.file = '../data/example.ms2'
        self.file2 = '../tmp/example2.ms2'

    def test_read_and_write(self):
        header = ms2.read_header(self.file)
        reader = ms2.read(self.file)

        with open(self.file2, 'w') as destination:
            ms2.write_header(destination, header)
            for psm in reader:
                ms2.write(destination, psm)


if __name__ == '__main__':
    unittest.main()

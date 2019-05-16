import unittest
import write_id_pepxml


class WriteIdPepxmlCase(unittest.TestCase):
    """Tests for 'write_id_pepxml.py'."""

    def setUp(self):
        self.id_file = '../data/Peaks_peptides.csv'
        self.pepxml_file = '../tmp/peaks_id.pepxml'

    def test_write_id_pepxml(self):
        write_id_pepxml.write(self.id_file, self.pepxml_file)


if __name__ == '__main__':
    unittest.main()

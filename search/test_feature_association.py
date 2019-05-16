import unittest
from pyteomics import mgf
import feature_association


class FeatureAssociationTestCase(unittest.TestCase):
    """Tests for 'feature_association.py'."""

    def setUp(self):
        self.ms_file = '../data/10222015_ABRF3_1000ng_68X_20181227.mgf'
        self.feature_file = '../data/10222015_ABRF3_1000ng_68X_features.csv'
        self.associated_ms_file = '../tmp/associated.mgf'

    def test_feature_association(self):
        spectra = feature_association.associate(self.ms_file, self.feature_file)

        with open(self.associated_ms_file, 'w') as destination:
            mgf.write(spectra, destination)


if __name__ == '__main__':
    unittest.main()

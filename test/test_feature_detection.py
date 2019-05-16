import unittest
import feature_detection
from pyopenms import FeatureXMLFile


class FeatureDetectionCase(unittest.TestCase):
    """Tests for 'feature_detection.py'."""

    def setUp(self):
        self.mzml_file = '../data/10222015_ABRF3_1000ng_68X.mzML'
        self.feature_file = '../tmp/10222015_ABRF3_1000ng_68X.featureXML'

    def test_feature_detection(self):
        features = feature_detection.run(self.mzml_file)

        features.setUniqueIds()
        fh = FeatureXMLFile()
        fh.store(self.feature_file, features)
        print("Found", features.size(), "features")


if __name__ == '__main__':
    unittest.main()

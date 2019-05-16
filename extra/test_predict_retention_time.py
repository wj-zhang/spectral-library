import unittest
import pandas as pd
import predict_retention_time


class PredictRetentionTimeTestCase(unittest.TestCase):
    """Tests for 'predict_retention_time.py'."""

    def setUp(self):
        self.in_file = '../tmp/workflow/library_id_test_11_files_decoy.csv'
        self.out_file = '../tmp/workflow/library_rt_id_test_11_files_decoy.csv'

    def test_predict_retention_time(self):
        ids = pd.read_csv(self.in_file)
        ids = predict_retention_time.predict(ids, fdr=0.0)
        ids.to_csv(self.out_file, header=True, index=False)


if __name__ == '__main__':
    unittest.main()

import unittest
import run_pipeline_two_rounds


class RunPipelineTwoRoundsTestCase(unittest.TestCase):
    """Tests for 'run_pipeline_two_rounds.py'."""

    def setUp(self):

        self.library_file = '../tmp/CharmeRT/simulation/data/simulation_library.ms2'
        self.ms_file = '../tmp/CharmeRT/simulation/data/simulation_query.mgf'
        self.id_file = '../tmp/CharmeRT/simulation/Ours/library_id_my_my.csv'

    def test_run_pipeline(self):
        print('library search ...')
        run_pipeline_two_rounds.run_main_library_search(self.ms_file, self.library_file, self.id_file, mode='SIMULATION_SEARCH')


if __name__ == '__main__':
    unittest.main()

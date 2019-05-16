import unittest
import evaluate_simulation_results


class EvaluateSimulationResultsTestCase(unittest.TestCase):

    def setUp(self):
        self.true_id_file = '../tmp/CharmeRT/simulation/data/HeLa_1ug_isow2_3hgradient_datasetA_chimera_0.1_simulation_id.csv'
        self.id_files = {
            # 'Ours': '../tmp/CharmeRT/simulation/Ours/library_id.csv',
            'Ours': '../tmp/joint_optimization_score/library_id_target_new_2_28_max.csv',
                        #  'SpectraST': '../tmp/CharmeRT/simulation/SpectraST/simulation_query.csv',
                        # 'MSPLIT': '../tmp/CharmeRT/simulation/MSPLIT/simulation_query.csv',
                         }

    def test_evaluate_simulation_results(self):
        for method, id_file in self.id_files.items():
            result = evaluate_simulation_results.evaluate(self.true_id_file, id_file)
            print(method)
            print(result)


if __name__ == '__main__':
    unittest.main()

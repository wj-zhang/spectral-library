import unittest
import test_build_library_from_splib
import test_generate_decoy_library
import test_search_library_from_PEAKS


def pipeline_unittest():
    suite = unittest.TestLoader().loadTestsFromModule(test_build_library_from_splib)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromModule(test_generate_decoy_library)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromModule(test_search_library_from_PEAKS)
    unittest.TextTestRunner(verbosity=2).run(suite)


if __name__ == '__main__':
    pipeline_unittest()

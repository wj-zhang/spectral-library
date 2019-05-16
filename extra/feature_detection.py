from pyopenms import *


def run(mzml_file):
    options = PeakFileOptions()
    options.setMSLevels([1])
    fh = MzMLFile()
    fh.setOptions(options)

    input_map = MSExperiment()
    fh.load(mzml_file, input_map)
    input_map.updateRanges()

    ff = FeatureFinder()
    ff.setLogType(LogType.CMD)

    name = "centroided"
    features = FeatureMap()
    seeds = FeatureMap()
    params = FeatureFinder().getParameters(name)
    ff.run(name, input_map, features, params, seeds)

    return features
import multiprocessing
import peaks_features
from itertools import repeat
from pyteomics import mgf, auxiliary


def associate(ms_file, feature_file, delta=1):
    features = peaks_features.read(feature_file)
    reader = mgf.read(ms_file)

    pool = multiprocessing.Pool(processes=8)
    spectra = pool.map(worker, zip(reader, repeat(features), repeat(delta)))

    return spectra


def worker(args):
    spectrum, features, delta = args
    mz, z = spectrum['params']['pepmass'][0], spectrum['params']['charge'][0].real
    rt = float(spectrum['params']['rtinseconds']) / 60

    mx, my = mz - delta, mz + delta
    features_associated = features[
        (features['RT Begin'] < rt) & (features['RT End'] > rt) & (features['m/z'] > mx) & (features['m/z'] < my)]

    mx, my = mz - 0.001, mz + 0.001
    features_associated.drop(features_associated[
                                 (features_associated['m/z'] > mx) & (features_associated['m/z'] < my) & (
                                     features_associated['z'] == z)].index, inplace=True)

    mzs, zs = list(features_associated['m/z']), list(features_associated['z'])
    mzs.append(mz)
    zs.append(z)

    spectrum['params']['pepmass'] = mzs
    spectrum['params']['charge'] = auxiliary.ChargeList([])
    for z in zs:
        spectrum['params']['charge'].append(auxiliary.Charge(z))

    return spectrum

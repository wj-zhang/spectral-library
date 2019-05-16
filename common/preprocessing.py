import numpy as np
from sklearn.preprocessing import normalize
from scipy.stats import rankdata
from fragments import ion_annotation
from pyteomics import mass, parser
from modifications import modifications_env


def preprocess(ms, noise_level=0.01, peptide=None, mod_env=None, ion_ppm=0.00002):
    if mod_env is None:
        mod_env = modifications_env([])

    masses = ms['m/z array']
    intensities = ms['intensity array']

    masses, intensities = mass_range(masses, intensities, ms)
    if not valid_spectrum(masses, peptide):
        return None, None

    masses, intensities = noise_remove(masses, intensities, peptide, noise_level)
    if not valid_spectrum(masses, peptide):
        return None, None

    if peptide:
        intensity_scale(masses, intensities, peptide, aa_mass=mod_env['aa_mass'], ion_ppm=ion_ppm)

    intensities = rank_normalize(intensities)

    return masses, intensities


def valid_spectrum(masses, peptide=None):
    peptide_length = 10 if peptide is None else parser.length(peptide)

    if len(masses) < peptide_length:
        return False
    if masses[-1] - masses[0] < 250:
        return False
    return True


def intensity_scale(masses, intensities, peptide, aa_mass=mass.std_aa_mass, ion_ppm=0.00002):
    ions = ion_annotation(masses, intensities, peptide, aa_mass=aa_mass, ion_ppm=ion_ppm)
    for i in ions.values():
        intensities[i] *= 5


def mass_range(masses, intensities, ms):
    mass_lower = 150
    try:
        mzs, zs = ms['params']['pepmass'], ms['params']['charge']
        mzs, zs = [mz for mz in mzs if mz is not None], [z.real for z in zs]
        mzz = [mz*z for mz, z in zip(mzs, zs)]
        mass_upper = max(mzz)
    except KeyError:
        mz = float(ms['params']['scan'][2])
        z = ms['params']['charge'][0]
        mass_upper = mz * z
    indices = (mass_lower < masses) & (masses < mass_upper)
    return masses[indices], intensities[indices]


def noise_remove(masses, intensities, peptide=None, noise_level=0.01):
    # max_num_peaks = 100 if peptide is None else 5*parser.length(peptide)

    base_peak_intensity = intensities.max(axis=0)
    relative_intensities = intensities / base_peak_intensity
    indices = relative_intensities > noise_level
    masses, intensities = masses[indices], intensities[indices]

    # if len(masses) > max_num_peaks:
    #     indices = np.sort(np.argpartition(intensities, -max_num_peaks)[-max_num_peaks:])
    #     masses, intensities = masses[indices], intensities[indices]

    return masses, intensities


def log_normalize(x):
    x = np.log(x)
    x = x.reshape(-1, 1)
    x = normalize(x, axis=0)
    return x.reshape(-1, )


def rank_normalize(x):
    x = rankdata(x)
    x = x.reshape(-1, 1)
    x = normalize(x, axis=0)
    return x.reshape(-1, )


def is_impure_spectrum(masses, intensities, peptide, aa_mass=mass.std_aa_mass, ion_ppm=0.00002, top_num=20, threshold=0.4):
    top_num = min(top_num, len(intensities))
    top_idx = np.argpartition(intensities, -top_num)[-top_num:]

    ions = ion_annotation(masses, intensities, peptide, aa_mass=aa_mass, ion_ppm=ion_ppm)
    annotated_idx = ions.values()

    annotated_top_idx = list(set(annotated_idx) & set(top_idx))
    return sum(intensities[annotated_top_idx]) / sum(intensities[top_idx]) < threshold

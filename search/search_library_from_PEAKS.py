import ms2
import pandas as pd
import configparser
from pyteomics import mgf, auxiliary
from dot_product import dot_product
from bisect import bisect_left, bisect_right
from preprocessing import preprocess
from modifications import modifications_env


def search(ms_file, library_file, fdr=0.01, mode='MAIN_SEARCH'):
    is_pre_search = False
    if 'PRE' in mode:
        is_pre_search = True

    configs = configparser.ConfigParser()
    configs.read('../common/settings.ini')
    precursor_ppm = configs.getfloat(mode, 'precursor_mass_error_tolerance_in_ppm') / 1000000.0
    ion_ppm = configs.getfloat(mode, 'fragment_ion_mass_error_tolerance_in_ppm') / 1000000.0

    spectrum_list, precursor_mz_list, mod_env = get_library_spectra(library_file)

    psms = []
    if is_pre_search:
        ms_dict = {}

    with mgf.read(ms_file) as ms_reader:
        for ms in ms_reader:
            ms, hit, score = spectrum_matching(ms, spectrum_list, precursor_mz_list, precursor_ppm, ion_ppm)

            if hit:
                decoy = hit['params']['seq'].startswith('#')
                peptide = hit['params']['modified seq']
                scan = int(ms['params']['scans'])
                index = int(ms['params']['title'][6:]) + 1
                rt = float(ms['params']['rtinseconds']) / 60
                library = hit['params']['scan'][0]
                protein = hit['params']['Protein']
                psms.append((index, scan, peptide, rt, decoy, score, library, protein))

                if is_pre_search:
                    ms['params']['sequence'] = peptide
                    ms_dict[scan] = ms

    psms = pd.DataFrame(psms, columns=['Index', 'Scan', 'Peptide', 'RT', 'Decoy', 'Score', 'Library', 'Protein'])
    if fdr > 0:
        psms = auxiliary.filter(psms, fdr=fdr, key='Score', reverse=True, is_decoy='Decoy')

    if is_pre_search:
        spectra = [ms_dict[scan] for scan in psms['Scan']]
        return spectra
    else:
        return psms


def get_library_spectra(library_file):
    header = ms2.read_header(library_file)
    mods = eval(header['Modifications'])
    mod_env = modifications_env(mods)

    library_reader = ms2.read(library_file)
    spectrum_list = list(library_reader)

    spectrum_list.sort(key=lambda x: float(x['params']['scan'][2]), reverse=False)
    precursor_mz_list = [float(x['params']['scan'][2]) for x in spectrum_list]

    return spectrum_list, precursor_mz_list, mod_env


def spectrum_matching(ms, spectrum_list, precursor_mz_list, precursor_ppm, ion_ppm):
    masses_ms, intensities_ms = preprocess(ms, ion_ppm=ion_ppm)
    if masses_ms is None:
        return None, None, None

    spectrum_candidates = get_spectrum_candidates(ms, spectrum_list, precursor_mz_list, precursor_ppm)
    if not spectrum_candidates:
        return None, None, None

    hit, score = score_matching_by_dot_product(masses_ms, intensities_ms, spectrum_candidates, ion_ppm)

    return ms, hit, score


def get_spectrum_candidates(ms, spectrum_list, precursor_mz_list, precursor_ppm):
    mz = ms['params']['pepmass'][0]
    z = ms['params']['charge'][0].real

    error = mz * precursor_ppm
    begin = bisect_left(precursor_mz_list, mz - error)
    end = bisect_right(precursor_mz_list, mz + error)

    spectrum_candidates = []
    for library in spectrum_list[begin: end]:
        if z == int(library['params']['charge'][0]):
            spectrum_candidates.append(library)

    return spectrum_candidates


def score_matching_by_dot_product(masses_ms, intensities_ms, spectrum_candidates, ion_ppm):
    score, hit, second_score = 0, None, 0
    for library in spectrum_candidates:
        masses_library = library['m/z array']
        intensities_library = library['intensity array']

        score_ = dot_product(masses_ms, intensities_ms, masses_library, intensities_library, ion_ppm=ion_ppm)
        if score_ > score:
            score, hit, second_score = score_, library, score
        elif score_ > second_score:
            second_score = score_

    score = adjust_score_by_second_hit(score, second_score)

    return hit, score


def adjust_score_by_second_hit(score, second_score):
    score = score - 0.4 * second_score
    return score

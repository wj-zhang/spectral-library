import ms2
import pandas as pd
import configparser
from pyteomics import mgf, auxiliary
import score_matching_by_nlc_branch_and_bound
from bisect import bisect_left, bisect_right
from preprocessing import preprocess
from modifications import modifications_env
import numpy as np


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
            ms, hits, score, weights, norms = spectrum_matching(ms, spectrum_list, precursor_mz_list, precursor_ppm, ion_ppm, mode)

            if hits:
                decoy = hits[0]['params']['seq'].startswith('#')
                peptide = ';'.join(hit['params']['modified seq'] for hit in hits)
                charge = ';'.join(str(int(hit['params']['charge'][0])) for hit in hits)
                scan = int(ms['params']['scans'])
                index = int(ms['params']['title'][6:]) + 1
                rt = float(ms['params']['rtinseconds']) / 60
                library = ';'.join(hit['params']['scan'][0] for hit in hits)
                protein = ';'.join(hit['params']['Protein'] for hit in hits)

                weights_str = ';'.join(str(weight) for weight in weights)
                norms_str = ';'.join(str(norm) for norm in norms)

                psms.append((index, scan, peptide, charge, rt, decoy, score, library, protein, weights_str, norms_str))

                if is_pre_search:
                    ms['params']['sequence'] = peptide
                    ms_dict[scan] = ms

    psms = pd.DataFrame(psms, columns=['Index', 'Scan', 'Peptide', 'Charge', 'RT', 'Decoy', 'Score', 'Library', 'Protein', 'Weights', 'Norms'])
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


def spectrum_matching(ms, spectrum_list, precursor_mz_list, precursor_ppm, ion_ppm, mode):

    if 'SIM' in mode:
        masses_ms = ms['m/z array']
        intensities_ms = ms['intensity array']
        intensities_ms /= np.linalg.norm(intensities_ms)
    else:
        masses_ms, intensities_ms = preprocess(ms, ion_ppm=ion_ppm)

    if masses_ms is None:
        return None, None, None, None, None
    ms['m/z array'], ms['intensity array'] = masses_ms, intensities_ms

    spectrum_candidates = get_spectrum_candidates(ms, spectrum_list, precursor_mz_list, precursor_ppm)

    if not spectrum_candidates:
        return None, None, None, None, None

    hits, score, weights, norms = score_matching_by_nlc_branch_and_bound.scoring(ms, spectrum_candidates, ion_ppm)

    if not hits:
        return None, None, None, None, None

    return ms, hits, score, weights, norms


def get_spectrum_candidates(ms, spectrum_list, precursor_mz_list, precursor_ppm):
    mzs, zs = ms['params']['pepmass'], ms['params']['charge']
    mzs, zs = [mz for mz in mzs if mz is not None], [z.real for z in zs]

    spectrum_candidates_list = []
    for mz, z in zip(mzs, zs):
        error = mz * precursor_ppm
        begin = bisect_left(precursor_mz_list, mz - error)
        end = bisect_right(precursor_mz_list, mz + error)
        spectrum_candidates = []
        for library in spectrum_list[begin: end]:
            if z == int(library['params']['charge'][0]):
                spectrum_candidates.append(library)
        if spectrum_candidates:
            spectrum_candidates_list.append(spectrum_candidates)

    return spectrum_candidates_list

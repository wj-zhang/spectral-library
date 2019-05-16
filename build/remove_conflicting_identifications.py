import ms2
import math
import configparser
from dot_product import dot_product
from bisect import bisect_left, bisect_right
from modifications import modifications_env
from pyteomics import mass
from fragments import ion_masses


def remove(library_file, similarity_cutoff=0.7, homology_cutoff=0.3, precursor_ppm=None):
    configs = configparser.ConfigParser()
    configs.read('../common/settings.ini')
    if not precursor_ppm:
        precursor_ppm = configs.getfloat('MAIN_SEARCH', 'precursor_mass_error_tolerance_in_ppm') / 1000000.0
    else:
        precursor_ppm /= 1000000.0

    ion_ppm = configs.getfloat('MAIN_SEARCH', 'fragment_ion_mass_error_tolerance_in_ppm') / 1000000.0

    spectrum_list, precursor_mz_list, mod_env = get_library_spectra(library_file)

    del_flag = [False] * len(spectrum_list)
    progress_last = 0
    for idx, ms in enumerate(spectrum_list):
        progress = int(idx / (len(spectrum_list)-1) * 100)
        if not progress == progress_last:
            print(str(progress), end='% ')
            progress_last = progress

        if del_flag[idx]:
            continue

        mz = float(ms['params']['scan'][2])
        z = int(ms['params']['charge'][0])

        error = mz * precursor_ppm
        begin = bisect_left(precursor_mz_list, mz - error)
        end = bisect_right(precursor_mz_list, mz + error)

        masses, intensities = ms['m/z array'], ms['intensity array']
        for idx_lib in range(begin, end):
            if idx == idx_lib:
                continue

            library = spectrum_list[idx_lib]
            if z == int(library['params']['charge'][0]):
                masses_lib, intensities_lib = library['m/z array'], library['intensity array']
                score = dot_product(masses, intensities, masses_lib, intensities_lib, ion_ppm=ion_ppm)
                if score > similarity_cutoff:
                    masses_theor, intensities_theor = \
                        theoretical_spectrum(ms['params']['modified seq'], mod_env['aa_mass'])
                    masses_theor_lib, intensities_theor_lib = \
                        theoretical_spectrum(library['params']['modified seq'], mod_env['aa_mass'])
                    score = dot_product(masses_theor, intensities_theor, masses_theor_lib, intensities_theor_lib,
                                        ion_ppm=ion_ppm)
                    if score < homology_cutoff:
                        if ms['params']['NReplicates'] < library['params']['NReplicates']:
                            del_flag[idx] = True
                        elif ms['params']['NReplicates'] > library['params']['NReplicates']:
                            del_flag[idx_lib] = True
                        else:
                            dot = dot_product(masses, intensities, masses_theor, intensities_theor, ion_ppm = ion_ppm)
                            dot_lib = dot_product(masses_lib, intensities_lib, masses_theor_lib, intensities_theor_lib, ion_ppm = ion_ppm)
                            if dot < dot_lib:
                                del_flag[idx] = True
                            else:
                                del_flag[idx_lib] = True

            if del_flag[idx]:
                break

    print('')
    return [ms for idx, ms in enumerate(spectrum_list) if del_flag[idx] is False]


def theoretical_spectrum(peptide, aa_mass=mass.std_aa_mass):
    b_ion_masses, y_ion_masses = ion_masses(peptide, aa_mass)
    masses = sorted(b_ion_masses + y_ion_masses)
    intensities = [1/math.sqrt(len(masses))] * len(masses)
    return masses, intensities


def get_library_spectra(library_file):
    header = ms2.read_header(library_file)
    mods = eval(header['Modifications'])
    mod_env = modifications_env(mods)

    library_reader = ms2.read(library_file)
    spectrum_list = list(library_reader)

    spectrum_list.sort(key=lambda x: float(x['params']['scan'][2]), reverse=False)
    precursor_mz_list = [float(x['params']['scan'][2]) for x in spectrum_list]

    return spectrum_list, precursor_mz_list, mod_env

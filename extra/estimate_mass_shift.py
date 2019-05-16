import os
import re
import peaks_id
from pyteomics import mgf
from preprocessing import preprocess
from peaks_modifications import mod_env, convert_modification_forward
from fragments import ion_annotation, ion_masses
import matplotlib.pyplot as plt


def fragment_ion_mass_shift(pepxml_file):
    mass_differences = []

    reader = read_psms(pepxml_file)
    for psm in reader:
        mass_differences += ion_mass_difference(psm, mod_env=mod_env)

    plt.hist(mass_differences, bins='auto', range=(-100, 100))
    plt.show()

    return sum(mass_differences)/len(mass_differences)


def ion_mass_difference(psm, mod_env):
    masses, intensities, peptide = psm['m/z array'], psm['intensity array'], psm['params']['modified seq']
    aa_mass = mod_env['aa_mass']
    ions = ion_annotation(masses, intensities, peptide, ion_ppm=0.00005, aa_mass=aa_mass)
    b_ion_masses, y_ion_masses = ion_masses(peptide, aa_mass=aa_mass)

    mass_differences = []
    for ion, idx in ions.items():
        ion_idx = int(ion[1:]) - 1
        if 'b' in ion:
            mass_differences.append((1 - b_ion_masses[ion_idx]/masses[idx])*1000000)
        elif 'y' in ion:
            mass_differences.append((1 - y_ion_masses[ion_idx]/masses[idx])*1000000)

    return mass_differences


def read_psms(pepxml_file):
    psms = peaks_id.read(pepxml_file, keep='first')
    psms = psms.groupby(['spectrum'])

    for ms_file, psm_df in psms:
        ms_file = re.sub(r'\.raw$', '', ms_file) + '.mgf'
        ms_file = os.path.join(os.path.dirname(pepxml_file), ms_file)
        ms2_spectra = get_ms2_spectra(ms_file)

        for idx, row in psm_df.iterrows():
            start_scan = row['start_scan']
            index = row['index']
            ms = ms2_spectra[index]

            mz = ms['params']['pepmass'][0]
            z = ms['params']['charge'][0].real
            rt = float(ms['params']['rtinseconds']) / 60
            modified_seq = convert_modification_forward(row['modified_peptide'])

            masses, intensities = preprocess(ms, peptide=modified_seq, mod_env=mod_env)

            if masses is None:
                continue

            params = {'scan': tuple([start_scan, start_scan, mz]),
                      'charge': tuple([z, mz]),
                      'RTime': rt,
                      'seq': row['peptide'],
                      'modified seq': modified_seq
                      }

            yield {'params': params, 'm/z array': masses, 'intensity array': intensities}


def get_ms2_spectra(mgf_file):
    ms2_spectra = {}
    with mgf.read(mgf_file) as reader:
        for spectrum in reader:
            index = int(spectrum['params']['title'][6:]) + 1
            ms2_spectra[index] = spectrum
    return ms2_spectra

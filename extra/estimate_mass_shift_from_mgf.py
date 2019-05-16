import numpy as np
from pyteomics import mgf, mass
from peaks_modifications import mod_env, convert_modification_forward
from fragments import ion_annotation, ion_masses
import matplotlib.pyplot as plt


def estimate_mass_shift(mgf_file, plot=True, psm_df=None):
    precursor_mass_differences = []
    ion_mass_differences = []

    if isinstance(mgf_file, str):
        reader = read_psms(mgf_file)
    else:
        if psm_df is None:
            reader = mgf_file
        else:
            reader = read_psms_by_df(mgf_file, psm_df)

    for psm in reader:
        precursor_mass_differences.append(precursor_mass_difference(psm, mod_env=mod_env))
        ion_mass_differences += ion_mass_difference(psm, mod_env=mod_env)

    if plot:
        plot_mass_difference(precursor_mass_differences, ion_mass_differences)

    mass_differences = np.asanyarray(precursor_mass_differences)
    precursor_mu, precursor_sigma = mass_differences.mean(), mass_differences.std()
    mass_differences = np.asanyarray(ion_mass_differences)
    ion_mu, ion_sigma = mass_differences.mean(), mass_differences.std()

    return precursor_mu, precursor_sigma, ion_mu, ion_sigma


def precursor_mass_difference(psm, mod_env):
    try:
        peptide, mz, z = psm['params']['modified seq'], psm['params']['charge'][1], psm['params']['charge'][0]
    except KeyError:
        peptide, mz, z = psm['params']['sequence'], psm['params']['pepmass'][0], psm['params']['charge'][0].real
    aa_mass = mod_env['aa_mass']
    theoretical_mz = mass.fast_mass2(peptide, aa_mass=aa_mass, charge=z)
    mass_difference = (1 - theoretical_mz / mz) * 1000000
    return mass_difference


def ion_mass_difference(psm, mod_env):
    try:
        masses, intensities, peptide = psm['m/z array'], psm['intensity array'], psm['params']['modified seq']
    except KeyError:
        masses, intensities, peptide = psm['m/z array'], psm['intensity array'], psm['params']['sequence']
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


def read_psms(mgf_file):
    with mgf.read(mgf_file) as reader:
        for spectrum in reader:
            if not spectrum['params']['sequence']:
                continue

            mz = spectrum['params']['pepmass'][0]
            z = spectrum['params']['charge'][0].real
            sequence = spectrum['params']['sequence']
            modified_seq = convert_modification_forward(sequence)
            masses = spectrum['m/z array']
            intensities = spectrum['intensity array']

            params = {'scan': tuple([None, None, mz]),
                      'charge': tuple([z, mz]),
                      'modified seq': modified_seq
                      }

            yield {'params': params, 'm/z array': masses, 'intensity array': intensities}


def plot_mass_difference(precursor_mass_differences, ion_mass_differences):
    mass_differences = np.asanyarray(precursor_mass_differences)
    mu, sigma = mass_differences.mean(), mass_differences.std()
    mu_sigma_str = '\n'.join([r'Precursor:', r'$\mu=%.2f$' % (mu, ), r'$\sigma=%.2f$' % (sigma, )])

    ax = plt.subplot(2, 1, 1)
    ax.hist(mass_differences, bins='auto', range=(-100, 100))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, mu_sigma_str, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    mass_differences = np.asanyarray(ion_mass_differences)
    mu, sigma = mass_differences.mean(), mass_differences.std()
    mu_sigma_str = '\n'.join([r'Fragment:', r'$\mu=%.2f$' % (mu, ), r'$\sigma=%.2f$' % (sigma, )])

    ax = plt.subplot(2, 1, 2)
    ax.hist(mass_differences, bins='auto', range=(-100, 100))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, mu_sigma_str, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    plt.show()


def read_psms_by_df(ms2_spectra, psm_df):
    for idx, row in psm_df.iterrows():
        index = row['index']
        spectrum = ms2_spectra[index]

        mz = spectrum['params']['pepmass'][0]
        z = spectrum['params']['charge'][0].real
        sequence = row['modified_peptide']
        modified_seq = convert_modification_forward(sequence)
        masses = spectrum['m/z array']
        intensities = spectrum['intensity array']

        params = {'scan': tuple([None, None, mz]),
                  'charge': tuple([z, mz]),
                  'modified seq': modified_seq
                  }

        yield {'params': params, 'm/z array': masses, 'intensity array': intensities}

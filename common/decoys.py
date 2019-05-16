import random
import numpy as np
from fragments import ion_masses, ion_annotation
from pyteomics.parser import parse, std_labels
from modifications import modifications_env


def shuffle(peptide, aa_labels=std_labels, keep_nterm=False, keep_cterm=False):
    peptide = parse(peptide, labels=aa_labels)

    start = 1 if keep_nterm else 0
    end = len(peptide)-1 if keep_cterm else len(peptide)

    _peptide_ = peptide[start:end]
    random.shuffle(_peptide_)
    return ''.join(peptide[:start]) + ''.join(_peptide_) + ''.join(peptide[end:])


def generate_decoy_spectrum(masses, intensities, peptide, mod_env=None, ion_ppm=0.00002):
    if mod_env is None:
        mod_env = modifications_env([])

    ion_flags = ion_annotation(masses, intensities, peptide, aa_mass=mod_env['aa_mass'], ion_ppm=ion_ppm)
    decoy_peptide = shuffle(peptide, aa_labels=mod_env['aa_labels'], keep_cterm=True, keep_nterm=False)
    b_ion_masses, y_ion_masses = ion_masses(decoy_peptide, aa_mass=mod_env['aa_mass'])

    idx_list = list(ion_flags.values())
    decoy_masses, decoy_intensities = np.delete(masses, idx_list).tolist(),  np.delete(intensities, idx_list).tolist()
    random.shuffle(decoy_intensities)

    decoy_ion_masses, decoy_ion_intensities = [], []
    for ion, idx in ion_flags.items():
        ion_type, ion_num = ion[0], int(ion[1:])-1
        decoy_ion_intensities.append(intensities[idx])
        if ion_type == 'b':
            decoy_ion_masses.append(b_ion_masses[ion_num])
        elif ion_type == 'y':
            decoy_ion_masses.append(y_ion_masses[ion_num])
        else:
            raise KeyError

    decoy_masses.extend(decoy_ion_masses), decoy_intensities.extend(decoy_ion_intensities)
    decoy_masses_intensities = sorted(zip(decoy_masses, decoy_intensities))
    decoy_masses_intensities = [list(t) for t in zip(*decoy_masses_intensities)]
    decoy_masses, decoy_intensities = decoy_masses_intensities[0], decoy_masses_intensities[1]

    return decoy_masses, decoy_intensities, decoy_peptide

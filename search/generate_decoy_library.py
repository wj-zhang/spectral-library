import re
import ms2
import configparser
from decoys import generate_decoy_spectrum
from modifications import modifications_env
import random


def generate(target_file, decoy_file, combined=False, random_state=0, ratio=1):
    configs = configparser.ConfigParser()
    configs.read('../common/settings.ini')
    ion_ppm = configs.getfloat('MAIN_SEARCH', 'fragment_ion_mass_error_tolerance_in_ppm') / 1000000.0

    header = ms2.read_header(target_file)
    mods = eval(header['Modifications'])
    mod_env = modifications_env(mods)

    reader = ms2.read(target_file)

    random.seed(random_state)
    with open(decoy_file, 'w') as destination:
        # if combined:
        ms2.write_header(destination, header)

        for psm in reader:
            if combined:
                ms2.write(destination, psm)

            decoy_peptides = []
            for i in range(ratio):
                while True:
                    masses, intensities, peptide = psm['m/z array'], psm['intensity array'], psm['params']['modified seq']
                    decoy_masses, decoy_intensities, decoy_peptide = generate_decoy_spectrum(masses, intensities, peptide,
                                                                                             mod_env=mod_env, ion_ppm=ion_ppm)
                    psm['params']['seq'] = '#' + ''.join(re.findall('[A-Z]+', decoy_peptide))
                    psm['params']['modified seq'] = decoy_peptide
                    decoy_psm = {'params': psm['params'], 'm/z array': decoy_masses, 'intensity array': decoy_intensities}

                    if decoy_peptide not in decoy_peptides:
                        break
                ms2.write(destination, decoy_psm)
                decoy_peptides.append(decoy_peptide)

import os
import re
import datetime
import string
import pandas as pd
import peaks_id
import configparser
import consensus_spectrum
from pyteomics import mgf
from statistics import median
from peaks_modifications import mod_env, convert_modification_forward
from preprocessing import preprocess


def create(pepxml_files, remove_duplicates=False):
    configs = configparser.ConfigParser()
    configs.read('../common/settings.ini')
    ion_ppm = configs.getfloat('MAIN_SEARCH', 'fragment_ion_mass_error_tolerance_in_ppm') / 1000000.0

    if isinstance(pepxml_files, str):
        pepxml_files = [pepxml_files]

    psms = None
    for pepxml_file in pepxml_files:
        psms_new = peaks_id.read(pepxml_file, keep=None, remove_duplicates=remove_duplicates)
        psms_new['dirname'] = os.path.dirname(pepxml_file)
        psms = pd.concat([psms, psms_new])

    psms = psms.groupby(['modified_peptide', 'assumed_charge'])

    translate_table = str.maketrans(dict.fromkeys(string.ascii_lowercase))

    ms_file_dict = {}
    for peptide, psm_df in psms:
        modified_seq = convert_modification_forward(peptide[0])
        seq = modified_seq.translate(translate_table)
        assumed_charge = peptide[1]
        proteins = psm_df['protein'].values[0]

        ms2_list = []
        for idx, row in psm_df.iterrows():
            ms_file = row['spectrum']
            ms_file = re.sub(r'\.raw$', '', ms_file) + '.mgf'
            ms_file = os.path.join(row['dirname'], ms_file)

            try:
                ms2_spectra = ms_file_dict[ms_file]
            except KeyError:
                ms2_spectra = get_ms2_spectra(ms_file)
                ms_file_dict[ms_file] = ms2_spectra

            index = row['index']
            ms = ms2_spectra[index]
            ms2_list.append(ms)

        ms_consensus = consensus_spectrum.create(ms2_list, ion_ppm)
        mz = median([ms['params']['pepmass'][0] for ms in ms2_list])
        ms_consensus['params']['pepmass'] = (mz, None)

        masses, intensities = preprocess(ms_consensus, peptide=modified_seq, mod_env=mod_env, ion_ppm=ion_ppm)
        if masses is None:
            continue

        start_scan = ';'.join([ms['params']['scans'] for ms in ms2_list])
        rt = median([float(ms['params']['rtinseconds']) / 60 for ms in ms2_list])
        params = {'scan': tuple([start_scan, start_scan, mz]),
                  'charge': tuple([assumed_charge, mz]),
                  'RTime': rt,
                  'seq': seq,
                  'modified seq': modified_seq,
                  'NReplicates': ms_consensus['params']['NReplicates'],
                  'Protein': proteins
                  }

        yield {'params': params, 'm/z array': masses, 'intensity array': intensities}


def create_header(pepxml_files):
    if isinstance(pepxml_files, str):
        pepxml_files = [pepxml_files]

    header = {
        'CreationDate': datetime.datetime.now().strftime("%c"),
        'CreationFile': ';'.join(pepxml_files),
        'Modifications': str(mod_env['digits'])
    }

    return header


def get_ms2_spectra(mgf_file):
    ms2_spectra = {}
    with mgf.read(mgf_file) as reader:
        for spectrum in reader:
            index = int(spectrum['params']['title'][6:]) + 1
            ms2_spectra[index] = spectrum
    return ms2_spectra

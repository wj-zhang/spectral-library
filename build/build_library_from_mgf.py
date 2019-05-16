import datetime
import string
import pandas as pd
import configparser
import consensus_spectrum
from pyteomics import mgf
from statistics import median
from peaks_modifications import mod_env, convert_modification_forward
from preprocessing import preprocess


def create(mgf_files):
    configs = configparser.ConfigParser()
    configs.read('../common/settings.ini')
    ion_ppm = configs.getfloat('MAIN_SEARCH', 'fragment_ion_mass_error_tolerance_in_ppm') / 1000000.0

    if isinstance(mgf_files, str):
        mgf_files = [mgf_files]

    index, ms2_spectra, psms = 0, {}, []
    for mgf_file in mgf_files:
        index, psms = get_ms2_spectra(mgf_file, index, ms2_spectra, psms)

    psms = pd.DataFrame(psms, columns=['index', 'modified_peptide', 'assumed_charge'])
    psms = psms.groupby(['modified_peptide', 'assumed_charge'])

    translate_table = str.maketrans(dict.fromkeys(string.ascii_lowercase))

    for peptide, psm_df in psms:
        modified_seq = convert_modification_forward(peptide[0])
        seq = modified_seq.translate(translate_table)
        assumed_charge = peptide[1]
        proteins = 'XXX'

        ms2_list = []
        for idx, row in psm_df.iterrows():
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
        rt = median([float(ms['params']['rtinseconds'])/60 for ms in ms2_list])
        try:
            irt = median([float(ms['params']['irt']) for ms in ms2_list])
            ccsvalue = median([float(ms['params']['ccsvalue']) for ms in ms2_list])
        except KeyError:
            irt, ccsvalue = 0, 0

        params = {'scan': tuple([start_scan, start_scan, mz]),
                  'charge': tuple([assumed_charge, mz]),
                  'RTime': rt,
                  'seq': seq,
                  'modified seq': modified_seq,
                  'NReplicates': ms_consensus['params']['NReplicates'],
                  'Protein': proteins,
                  'iRT': irt,
                  'CCS': ccsvalue
                  }

        yield {'params': params, 'm/z array': masses, 'intensity array': intensities}


def create_header(mgf_files):
    if isinstance(mgf_files, str):
        mgf_files = [mgf_files]

    header = {
        'CreationDate': datetime.datetime.now().strftime("%c"),
        'CreationFile': ';'.join(mgf_files),
        'Modifications': str(mod_env['digits'])
    }

    return header


def get_ms2_spectra(mgf_file, index, ms2_spectra, psms):
    with mgf.read(mgf_file) as reader:
        for spectrum in reader:
            if not spectrum['params']['sequence']:
                continue

            ms2_spectra[index] = spectrum
            psms.append((index, spectrum['params']['sequence'], spectrum['params']['charge'][0].real))
            index += 1

    print('Totally {} psms loaded after {}'.format(index, mgf_file))
    return index, psms

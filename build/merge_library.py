import ms2
import string
import pandas as pd
import consensus_spectrum
from statistics import median
from preprocessing import preprocess
from modifications import modifications_env
from peaks_modifications import mod_env, convert_modification_forward


def merge(library_files):
    translate_table = str.maketrans(dict.fromkeys(string.ascii_lowercase))

    library_header, library_spectra, library_index = None, None, 0
    for file in library_files:
        header, spectra, psms, index = get_library_spectra(file, library_index)

        if library_spectra is None:
            library_header, library_spectra, library_psms, library_index = header, spectra, psms, index
            continue

        merged_psms = library_psms + psms

        merged_psms = pd.DataFrame(merged_psms, columns=['index', 'modified_peptide', 'assumed_charge'])
        merged_psms = merged_psms.groupby(['modified_peptide', 'assumed_charge'])

        for peptide, psm_df in merged_psms:
            if len(psm_df.index) == 1:
                idx = psm_df['index'].iloc[0]
                spectrum = library_spectra[idx] if idx < library_index else spectra[idx-library_index]
                yield spectrum

            elif len(psm_df.index) == 2:
                modified_seq = peptide[0]
                seq = modified_seq.translate(translate_table)
                assumed_charge = peptide[1]
                proteins = 'XXX'

                ms2_list = []
                for idx, row in psm_df.iterrows():
                    idx = row['index']
                    spectrum = library_spectra[idx] if idx < library_index else spectra[idx - library_index]
                    ms2_list.append(spectrum)

                start_scan = ';'.join([ms['params']['scan'][0] for ms in ms2_list])
                mz = median([float(ms['params']['scan'][2]) for ms in ms2_list])
                rt = median([float(ms['params']['RTime']) for ms in ms2_list])
                irt = median([float(ms['params']['iRT']) for ms in ms2_list])
                ccsvalue = median([float(ms['params']['CCS']) for ms in ms2_list])

                ion_ppm = 25 /1000000.0 ###
                ms_consensus = consensus_spectrum.create(ms2_list, ion_ppm)
                ms_consensus['params']['scan'] = tuple([start_scan, start_scan, mz])
                ms_consensus['params']['charge'] = tuple([float(assumed_charge), mz])

                masses, intensities = preprocess(ms_consensus, peptide=modified_seq, mod_env=mod_env, ion_ppm=ion_ppm)
                if masses is None:
                    continue

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
            else:
                raise ValueError


def get_library_spectra(library_file, library_index):
    header = ms2.read_header(library_file)
    spectra = list(ms2.read(library_file))

    psms, index = [], library_index
    for spectrum in spectra:
        psms.append((index, spectrum['params']['modified seq'], spectrum['params']['charge'][0]))
        index += 1

    return header, spectra, psms, index - library_index

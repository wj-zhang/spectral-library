from pyteomics import mgf
import pandas as pd


def combine(in_file):
    spectra = []
    with mgf.read(in_file) as ms_reader:
        for ms in ms_reader:
            spectra.append([ms['params']['scans'], ms])
    spectra_df = pd.DataFrame(spectra, columns=['scan', 'ms'])

    spectra_df = spectra_df.groupby(['scan'])
    combined_spectra = []
    for _, group in spectra_df:
        combined_ms = None
        for _, row in group.iterrows():
            ms = row['ms']
            if combined_ms is None:
                combined_ms = ms
                pepmass = list(combined_ms['params']['pepmass'])
                pepmass.remove(None)
                combined_ms['params']['pepmass'] = tuple(pepmass)
            else:
                combined_ms['params']['charge'].append(ms['params']['charge'][0])
                combined_ms['params']['pepmass'] += (ms['params']['pepmass'][0],)
        combined_spectra.append(combined_ms)

    return combined_spectra

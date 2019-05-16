from pyteomics import mgf, auxiliary
import pandas as pd


def combine(clone_ms_file, single_ms_file):
    spectra = []
    with mgf.read(clone_ms_file) as ms_reader:
        for ms in ms_reader:
            spectra.append([ms['params']['scans'], ms['params']['pepmass'][0], ms['params']['charge'][0]])
    spectra_df = pd.DataFrame(spectra, columns=['scan', 'pepmass', 'charge'])
    spectra_df = spectra_df.groupby(['scan'])

    combined_spectra = []
    for _, group in spectra_df:
        combined_ms = []
        for _, row in group.iterrows():
            if not combined_ms:
                combined_ms.append(row['scan'])
                combined_ms.append(tuple([row['pepmass']]))
                combined_ms.append(auxiliary.ChargeList([]))
                combined_ms[2].append(row['charge'])
            else:
                combined_ms[2].append(row['charge'])
                combined_ms[1] += (row['pepmass'],)
        combined_spectra.append(combined_ms)

    combined_spectra_df = pd.DataFrame(combined_spectra, columns=['scan', 'pepmass', 'charge'])

    spectra = []
    with mgf.read(single_ms_file) as ms_reader:
        for ms in ms_reader:
            spectra.append([ms['params']['scans'], ms])
    spectra_df = pd.DataFrame(spectra, columns=['scan', 'ms'])

    spectra_df = spectra_df.join(combined_spectra_df.set_index('scan'), on='scan')
    combined_spectra = []
    for _, row in spectra_df.iterrows():
        if len(row['charge']) > 1:
            row['ms']['params']['pepmass'] = row['pepmass']
            row['ms']['params']['charge'] = row['charge']
        elif len(row['charge']) == 1:
            assert abs(row['ms']['params']['charge'][0] - row['charge'][0]) < 1e-4
        combined_spectra.append(row['ms'])

    return combined_spectra

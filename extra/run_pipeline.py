import pandas as pd
import ms2
import build_library_from_mgf
import remove_conflicting_identifications
import generate_decoy_library
import search_library_from_PEAKS
import predict_retention_time
from peaks_modifications import convert_modification_backward


def run_build_library_from_mgf(mgf_files, library_file):
    header = build_library_from_mgf.create_header(mgf_files)
    creator = build_library_from_mgf.create(mgf_files)

    with open(library_file, 'w') as destination:
        ms2.write_header(destination, header)
        for psm in creator:
            ms2.write(destination, psm)


def run_remove_conflicting_identifications(library_file):
    header = ms2.read_header(library_file)
    spectrum_list = remove_conflicting_identifications.remove(library_file)

    with open(library_file, 'w') as destination:
        ms2.write_header(destination, header)
        for psm in spectrum_list:
            ms2.write(destination, psm)


def run_generate_decoy_library(target_file, decoy_file):
    generate_decoy_library.generate(target_file, decoy_file, combined=True)


def run_search_library_from_PEAKS(ms_file, library_file, id_file):
    psms = search_library_from_PEAKS.search(ms_file, library_file, fdr=0.0)
    psms['Peptide'] = psms['Peptide'].apply(convert_modification_backward)
    psms.to_csv(id_file, header=True, index=False)


def run_predict_retention_time(id_file):
    ids = pd.read_csv(id_file)
    ids = predict_retention_time.predict(ids, fdr=0.0)
    ids.to_csv(id_file, header=True, index=False)

import pandas as pd
import ms2
import build_library_from_mgf
import remove_conflicting_identifications
import generate_decoy_library
import search_library_from_PEAKS
import search_library_from_mgf
import predict_retention_time
from peaks_modifications import convert_modification_backward
import correct_mass_shifts
import estimate_mass_shift_from_mgf
from pyteomics import mgf, auxiliary, pepxml
import combine_clone_mgf
import build_library_from_pepxml
import combine_clone_from_two_mgfs
import re
import os
import numpy as np


def run_correct_mass_shifts(mgf_files, corrected_file, mode='library'):
    if 'library' in mode:
        with open(corrected_file, 'w') as destination:
            for mgf_file in mgf_files:
                precursor_mu, _, ion_mu, _ = estimate_mass_shift_from_mgf.estimate_mass_shift(mgf_file, plot=False)
                spectra = correct_mass_shifts.correct(mgf_file, precursor_mu, ion_mu, mode=mode)
                mgf.write(spectra, destination)
    elif 'query' in mode:
        mgf_file, spectra = mgf_files
        with open(corrected_file, 'w') as destination:
            precursor_mu, _, ion_mu, _ = estimate_mass_shift_from_mgf.estimate_mass_shift(spectra, plot=False)
            spectra = correct_mass_shifts.correct(mgf_file, precursor_mu, ion_mu, mode=mode)
            mgf.write(spectra, destination)


def run_build_library_from_mgf(mgf_files, library_file):
    header = build_library_from_mgf.create_header(mgf_files)
    creator = build_library_from_mgf.create(mgf_files)

    with open(library_file, 'w') as destination:
        ms2.write_header(destination, header)
        for psm in creator:
            ms2.write(destination, psm)


def run_remove_conflicting_identifications(library_file, library_file_new=None, precursor_ppm=None):
    if not library_file_new:
        library_file_new = library_file

    header = ms2.read_header(library_file)
    spectrum_list = remove_conflicting_identifications.remove(library_file, precursor_ppm=precursor_ppm)

    with open(library_file_new, 'w') as destination:
        ms2.write_header(destination, header)
        for psm in spectrum_list:
            ms2.write(destination, psm)


def run_generate_decoy_library(target_file, decoy_file, random_state=0, ratio=1, combined=False):
    generate_decoy_library.generate(target_file, decoy_file, combined=combined, random_state=random_state, ratio=ratio)


def run_pre_library_search(ms_file, library_file):
    spectra = search_library_from_PEAKS.search(ms_file, library_file, fdr=0.001, mode='PRE_SEARCH')
    return spectra


def run_separate_target_decoy(library_file, library_target_file, library_decoy_file):
    header = ms2.read_header(library_file)
    library_reader = ms2.read(library_file)
    spectrum_list = list(library_reader)

    spectra_target, spectra_decoy = [], []
    for spectrum in spectrum_list:
        if spectrum['params']['seq'].startswith('#'):
            spectra_decoy.append(spectrum)
        else:
            spectra_target.append(spectrum)

    with open(library_target_file, 'w') as destination:
        ms2.write_header(destination, header)
        for spectrum in spectra_target:
            ms2.write(destination, spectrum)

    with open(library_decoy_file, 'w') as destination:
        ms2.write_header(destination, header)
        for spectrum in spectra_decoy:
            ms2.write(destination, spectrum)


def run_combine_clone_mgf(ms_file):
    spectra = combine_clone_mgf.combine(ms_file)

    with open(ms_file, 'w') as destination:
        mgf.write(spectra, destination)


def run_main_library_search(ms_file, library_file, id_file, mode='MAIN_SEARCH'):
    psms = search_library_from_mgf.search(ms_file, library_file, fdr=0.0, mode=mode)
    psms['Peptide'] = psms['Peptide'].apply(convert_modification_backward)
    psms.to_csv(id_file, header=True, index=False)


def run_merge_target_decoy_ids(id_target_file, id_decoy_file, id_file, fdr=0.01):
    target_ids, decoy_ids = pd.read_csv(id_target_file), pd.read_csv(id_decoy_file)
    ids = target_ids.append(decoy_ids, ignore_index=True)
    ids.sort_values(by=['Score'], ascending=False, inplace=True)
    ids.drop_duplicates(subset=['Scan'], keep='first', inplace=True)
    ids = auxiliary.filter(ids, fdr=fdr, key='Score', reverse=True, is_decoy='Decoy')
    psm_num = 0
    peptides = []
    for idx, id in ids.iterrows():
        peptides_id = list(zip(id['Peptide'].split(';'), id['Charge'].split(';')))
        psm_num += len(peptides_id)
        peptides += peptides_id

    print(len(ids))
    print(str(psm_num))
    print(len(set(peptides)))

    ids.to_csv(id_file, header=True, index=False)


def run_cut_off_ids(raw_id_file, id_file, threshold=0.15):
    raw_ids = pd.read_csv(raw_id_file)
    for idx, id in raw_ids.iterrows():
        weights = [float(weight) for weight in id['Weights'].split(';')]
        w_candidates = list(zip(weights, list(range(len(weights)))))
        w_candidates.sort(key=lambda x: x[0], reverse=True)
        w_candidates = list(zip(*w_candidates))
        w, candidates = w_candidates[0], list(w_candidates[1])
        w_percent = np.array(w) / np.max(w)
        idx_cut = np.argmax(w_percent < threshold)
        idx_cut = len(w_percent) if idx_cut == 0 else idx_cut
        hits = candidates[0:idx_cut]

        peptides = id['Peptide'].split(';')
        peptides = [peptides[hit] for hit in hits]

        charges = id['Charge'].split(';')
        charges = [charges[hit] for hit in hits]

        proteins = id['Protein'].split(';')
        proteins = [proteins[hit] for hit in hits]

        raw_ids.at[idx, 'Peptide'] = ';'.join(peptides)
        raw_ids.at[idx, 'Charge'] = ';'.join(charges)
        raw_ids.at[idx, 'Protein'] = ';'.join(proteins)

    raw_ids = raw_ids[['Index', 'Scan', 'Peptide', 'RT', 'Charge', 'Decoy', 'Score', 'Library', 'Protein']]
    raw_ids.to_csv(id_file, header=True, index=False)


def run_format_to_PEAKS_id(id_file, peaks_id_file, clone_ms_file):
    scan_index_feature = []
    with mgf.read(clone_ms_file) as ms_reader:
        for ms in ms_reader:
            scan_index_feature.append([ms['params']['scans'], int(ms['params']['title'][6:]) + 1])

    scan_index_feature_df = pd.DataFrame(scan_index_feature, columns=['Scan', 'Index'])
    scan_index_feature_df = scan_index_feature_df.groupby(['Scan'])

    scan_index_feature = []
    for scan, group in scan_index_feature_df:
        indexs = []
        for _, row in group.iterrows():
            indexs.append(str(row['Index']))
        scan_index_feature.append([int(scan), ';'.join(indexs)])

    scan_index_feature_df = pd.DataFrame(scan_index_feature, columns=['Scan', 'Index_feature'])
    raw_ids = pd.read_csv(id_file)
    raw_ids = raw_ids.join(scan_index_feature_df.set_index('Scan'), on='Scan')

    peaks_ids = []
    for idx, id in raw_ids.iterrows():
        peptides = id['Peptide'].split(';')
        indexs = id['Index_feature'].split(';')
        for i in range(len(peptides)):
            peaks_ids.append([indexs[i], id['Scan'], peptides[i], id['RT'], id['Decoy'], id['Score'], None, None])

    peaks_ids = pd.DataFrame(peaks_ids, index=None, columns=['Index', 'Scan', 'Peptide', 'RT', 'Decoy', 'Score', 'Library', 'Protein'])
    peaks_ids.to_csv(peaks_id_file, header=True, index=False)


def run_predict_retention_time(id_file):
    ids = pd.read_csv(id_file)
    ids = predict_retention_time.predict(ids, fdr=0.0)
    ids.to_csv(id_file, header=True, index=False)


def run_build_library_from_pepxml(pepxml_files, library_file, remove_duplicates=False):
    header = build_library_from_pepxml.create_header(pepxml_files)
    creator = build_library_from_pepxml.create(pepxml_files, remove_duplicates=remove_duplicates)

    with open(library_file, 'w') as destination:
        ms2.write_header(destination, header)
        for psm in creator:
            ms2.write(destination, psm)


def run_combine_clone_from_two_mgfs(single_ms_file, clone_ms_file, ms_file):
    spectra = combine_clone_from_two_mgfs.combine(clone_ms_file, single_ms_file)

    with open(ms_file, 'w') as destination:
        mgf.write(spectra, destination)


def run_correct_mass_shifts_from_pepxml(pepxml_file, corrected_dir, corrected_suffix, remove_duplicates=False):

    with open(pepxml_file, 'r') as reader:
        content = reader.read()

    psms = pepxml.DataFrame(pepxml_file)
    if remove_duplicates:
        psms.drop_duplicates(subset=['spectrum', 'start_scan'], keep=False, inplace=True)

    psms = psms.groupby(['spectrum'])

    for ms_file, psm_df in psms:
        filename_without_ext = re.sub(r'\.raw$', '', ms_file)
        mgf_file = os.path.join(os.path.dirname(pepxml_file), filename_without_ext + '.mgf')
        spectra = _get_ms2_spectra(mgf_file)
        # psm_df_new = psm_df.sort_values(['-10lgP'], ascending=False, inplace=False)
        # psm_df_new.drop_duplicates(subset=['modified_peptide', 'assumed_charge'], keep='first', inplace=True)
        precursor_mu, _, ion_mu, _ = estimate_mass_shift_from_mgf.estimate_mass_shift(spectra, plot=False, psm_df=psm_df)
        spectra = correct_mass_shifts.correct_from_pepxml(spectra, psm_df['index'].values, precursor_mu, ion_mu)
        corrected_file = os.path.join(corrected_dir, filename_without_ext + corrected_suffix + '.mgf')
        with open(corrected_file, 'w') as destination:
            mgf.write(spectra, destination)

        content = re.sub(ms_file, filename_without_ext + corrected_suffix + '.raw', content)

    correct_pepxml_file = os.path.join(corrected_dir, re.sub(r'\.pep\.xml', '', os.path.basename(pepxml_file))
                                       + corrected_suffix + '.pep.xml')
    with open(correct_pepxml_file, 'w') as writer:
        writer.write(content)

    return correct_pepxml_file


def _get_ms2_spectra(mgf_file):
    ms2_spectra = {}
    with mgf.read(mgf_file) as reader:
        for spectrum in reader:
            index = int(spectrum['params']['title'][6:]) + 1
            ms2_spectra[index] = spectrum
    return ms2_spectra

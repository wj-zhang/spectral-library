import subprocess
import csv
import pandas
import splib_modifications
import peaks_modifications
import os


def run_build_library(st_library_file):
    library_file = st_library_file.replace('.sptxt', '.map')
    if os.path.isfile(library_file):
        os.remove(library_file)
    # subprocess.check_call(["java", '-jar', 'SpectrumLib.jar', library_file, "./null", "2", './null'])


def run_library_search(library_file, ms_file, raw_id_file, mass_tol=2):
    subprocess.check_call(["java", '-jar', 'SpectrumLib.jar', library_file, ms_file, str(mass_tol), raw_id_file])


def run_fdr(raw_id_file, fdr_id_file, svm_tmp_dir, fdr=0.01):
    if not os.path.isdir(svm_tmp_dir):
        os.mkdir(svm_tmp_dir)
    subprocess.check_call(["java", '-jar', 'MixtureTDA.jar', raw_id_file, fdr_id_file, str(fdr), svm_tmp_dir])


def run_fdr_by_threshold(raw_id_file, fdr_id_file, svm_tmp_dir, threshold1=0, threshold2=0):
    if not os.path.isdir(svm_tmp_dir):
        os.mkdir(svm_tmp_dir)
    subprocess.check_call(["java", '-jar', 'MixtureTDA.jar', raw_id_file, fdr_id_file, str(threshold1), str(threshold2), svm_tmp_dir])


def run_convert_raw_id_to_csv_file(raw_id_file, csv_file, decoy, column_name='svm1-score'):
    with open(raw_id_file, 'r') as file:
        id_reader = csv.reader(file, delimiter='\t')
        id_list = []
        columns = dict()
        for id in id_reader:
            if id[0].startswith('#'):
                for idx, key in enumerate(id):
                    columns[key] = idx
                continue

            index = int(id[columns['Scan#']].split(';')[0][1:]) + 1
            score = float(id[columns[column_name]])
            id_list.append([index, score, decoy, float(id[columns['alpha']])])

    id_pd = pandas.DataFrame(id_list, columns=['Scan', 'Score', 'Decoy', 'Alpha'])
    id_pd.to_csv(csv_file, header=True, index=False)


def run_format_to_csv_file(fdr_id_file, id_file):
    with open(fdr_id_file, 'r') as csv_file:
        id_reader = csv.reader(csv_file, delimiter='\t')
        id_list = []
        columns = dict()
        for id in id_reader:
            if id[0].startswith('#'):
                for idx, key in enumerate(id):
                    columns[key] = idx
                continue

            index = int(id[columns['Scan#']].split(';')[0][1:]) + 1
            peptides = id[columns['Annotation']].split('!')
            charges = id[columns['Charge']].split('!')
            assert len(peptides) == len(charges)

            proteins = id[columns['Protein']].split('!')
            proteins = ['[' + protein[2:].replace('/', ',') + ']' for protein in proteins]

            for idx in range(len(peptides)):
                peptides[idx] = peaks_modifications.convert_modification_backward(
                    splib_modifications.convert_modification_forward(peptides[idx]))

            id_list.append([index, ';'.join(peptides), ';'.join(charges), ';'.join(proteins)])

    id_pd = pandas.DataFrame(id_list, columns=['Index', 'Peptide', 'Charge', 'Protein'])
    id_pd.to_csv(id_file, header=True, index=False)


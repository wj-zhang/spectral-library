import format_PEAKS_pepxml_for_SpectraST
import subprocess
import csv
import pandas
from pyteomics import auxiliary
import os
import splib_modifications
import peaks_modifications


def run_format_pepxml(in_file, out_file):
    format_PEAKS_pepxml_for_SpectraST.run(in_file, out_file)


def run_build_library(pepxml_file, library_file):
    subprocess.check_call(["spectrast", "-cN" + library_file, pepxml_file])


def run_build_consensus_library(library_file, out_library_file):
    subprocess.check_call(["spectrast", "-cJU", "-cAC", "-cN" + out_library_file, library_file + ".splib"])


def run_quality_control_library(library_file, out_library_file):
    subprocess.check_call(["spectrast", "-cAQ", "-cN" + out_library_file, library_file + ".splib"])


def run_generate_decoy_library(library_file, out_library_file, ratio=1):
    subprocess.check_call(["spectrast", "-cAD", "-cc", "-cy" + str(ratio), "-cN" + out_library_file, library_file + ".splib"], )


def run_library_search(library_file, ms_file, id_file, mass_tol=3):
    id_dirname = os.path.dirname(id_file)
    id_filename = os.path.basename(id_file)
    id_extension = os.path.splitext(id_filename)[1][1:]
    subprocess.check_call(
        ["spectrast", '-sM' + str(mass_tol), '-sE' + id_extension, '-sO' + id_dirname, "-sL" + library_file + ".splib", ms_file])
    st_id_file = id_dirname + '/' + os.path.splitext(os.path.basename(ms_file))[0] + '.' + id_extension

    if os.path.isfile(id_file):
        os.remove(id_file)
    os.rename(st_id_file, id_file)


def run_fdr(id_file, fdr_id_file, fdr=0.01):
    with open(id_file, 'r') as csv_file:
        id_reader = csv.reader(csv_file, delimiter='\t')
        id_list = []
        columns = dict()
        for id in id_reader:
            if id[0].startswith('#'):
                for idx, key in enumerate(id):
                    columns[key] = idx
                continue

            index = int(id[columns['### Query']].split('=')[1]) + 1
            peptide, charge = id[columns['ID']].split('/')
            peptide = peaks_modifications.convert_modification_backward(
                splib_modifications.convert_modification_forward(peptide))
            charge = str(int(charge))
            score = float(id[columns['Fval']])
            decoy = id[columns['Proteins']].startswith('DECOY')
            protein = '[' + id[columns['Proteins']].replace(';', ',') + ']'
            id_list.append([index, peptide, charge, decoy, score, protein])

    id_pd = pandas.DataFrame(id_list, columns=['Index', 'Peptide', 'Charge', 'Decoy', 'Score', 'Protein'])
    ids = auxiliary.filter(id_pd, fdr=fdr, key='Score', reverse=True, is_decoy='Decoy')
    ids.to_csv(fdr_id_file, header=True, index=False)

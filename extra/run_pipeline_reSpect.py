import subprocess
from pyteomics import pepxml, auxiliary
import peaks_modifications
import splib_modifications
import pandas


def run_respect(pepxml_file):
    subprocess.check_call(['RespectParser', pepxml_file])


def run_peptide_i_prophet(in_pepxml_file, out_pepxml_file):
    subprocess.check_call(['xinteract', '-N'+out_pepxml_file, '-p0.05', '-l6', '-PPM', '-O', '-dDECOY', '-i', in_pepxml_file])


def run_pepxml_fdr(pepxml_file, id_file, fdr=0.01, simulation=False):
    id_list = []
    psmReader = pepxml.DataFrame(pepxml_file)
    for _, row in psmReader.iterrows():
        peptide = peaks_modifications.convert_modification_backward(splib_modifications.convert_modification_forward(row['modified_peptide']))
        charge = str(int(row['charge']))
        score = row['fval']
        decoy = row['lib_remark'].startswith('DECOY')
        protein = '[' + ','.join(row['protein']) + ']'
        if simulation:
            index = int(row['retention_time_sec'])
        else:
            index = row['start_scan']
        id_list.append([index, peptide, charge, decoy, score, protein])

    id_pd = pandas.DataFrame(id_list, columns=['Index', 'Peptide', 'Charge', 'Decoy', 'Score', 'Protein'])
    ids = auxiliary.filter(id_pd, fdr=fdr, key='Score', reverse=True, is_decoy='Decoy')
    ids.to_csv(id_file, header=True, index=False)


def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start

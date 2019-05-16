import os
import pandas as pd
from subprocess import Popen, PIPE
from pyteomics import auxiliary


def predict(ids, fdr):
    rt_file = '../tmp/rt_tmp.csv'
    export_for_rt_prediction(ids, rt_file)
    run_rt_prediction(rt_file)
    ids = import_ids_from_rt_prediction(rt_file)

    if fdr > 0:
        ids = auxiliary.filter(ids, fdr=fdr, key='score', reverse=True, is_decoy='IsDecoy')
    else:
        ids = import_ids_to_Peaks(ids)

    return ids


def export_for_rt_prediction(ids, rt_file):

    def format_modification(row):
        seq = row['Peptide']
        seq = seq.replace('C(+57.02)', 'C(Carbamidomethylation)')
        seq = seq.replace('M(+15.99)', 'M(Oxidation (M))')
        seq = seq.replace('N(+.98)', 'N(Deamidation (NQ))')
        seq = seq.replace('Q(+.98)', 'Q(Deamidation (NQ))')
        return seq

    ids['Peptide'] = ids.apply(format_modification, axis=1)

    ids.rename(columns={'Score': '-10lgP', 'Decoy': 'IsDecoy'}, inplace=True)
    ids.to_csv(rt_file, header=True, index=False)


def run_rt_prediction(rt_file):
    rt_prediction_exe = 'C:/Users/xchen/PycharmProjects/RetentionTimePred/dist/main/main.exe'
    rt_file = os.path.abspath(rt_file).replace('\\', '/')

    command_line = [rt_prediction_exe, '--train_regression_file', rt_file, '--pred_rt_file', rt_file]
    print(command_line)
    process = Popen(command_line, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    print(' STDOUT: ' + stdout.decode('UTF-8'))
    print(' STDERR: ' + stderr.decode('UTF-8'))


def import_ids_from_rt_prediction(rt_file):
    rt_file = rt_file.replace('.csv', '-RTpred-RTldf.csv')
    ids = pd.read_csv(rt_file)
    return ids


def import_ids_to_Peaks(ids):
    ids = ids[['Index', 'Scan', 'Peptide', 'IsDecoy', 'score']]
    ids.rename(columns={'IsDecoy': 'Decoy', 'score': 'Score'}, inplace=True)

    def format_modification(row):
        seq = row['Peptide']
        seq = seq.replace('C(Carbamidomethylation)', 'C(+57.02)')
        seq = seq.replace('M(Oxidation (M))', 'M(+15.99)')
        seq = seq.replace('N(Deamidation (NQ))', 'N(+.98)')
        seq = seq.replace('Q(Deamidation (NQ))', 'Q(+.98)')
        return seq

    ids['Peptide'] = ids.apply(format_modification, axis=1)

    return ids

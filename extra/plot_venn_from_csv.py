import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import os
import re


def process_protein_id(protein):
    protein = re.sub(r'\[|\]|\'| ', '', protein)
    protein = re.sub(';', ',', protein)
    return list(set(protein.split(',')))


def process_peptide_id(peptide):
    return peptide.split(';')


def get_results_from_csv(filename, column_names):
    ids = pd.read_csv(filename)
    items = []
    if 'Peptide' in column_names:
        for _, row in ids.iterrows():
            peps = list(set(row['Peptide'].split(';')))
            items += zip([row['Index']]*len(peps), peps)
        return pd.DataFrame(items, columns=['Index', 'Peptide'])
    elif 'Protein' in column_names:
        for _, row in ids.iterrows():
            proteins = process_protein_id(row['Protein'])
            items.extend(zip([row['Index']]*len(proteins), proteins))
        return pd.DataFrame(items, columns=['Index', 'Protein'])
    else:
        raise KeyError


def get_results(filename, column_names):
    file_extension = os.path.splitext(filename)[1]
    if file_extension == '.csv':
        return get_results_from_csv(filename, column_names)


def plot_venn(files, column_names):
    filenames = list(files.values())
    sets = [set([tuple(item) for item in get_results(filename, column_names)[column_names].values.tolist()])
            for filename in filenames]

    if len(sets) == 2:
        venn2(sets, set_labels=list(files.keys()))
    elif len(sets) == 3:
        venn3(sets, set_labels=list(files.keys()))
    plt.show()

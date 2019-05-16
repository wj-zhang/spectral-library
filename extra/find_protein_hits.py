import pandas as pd
from pyteomics import fasta
import re as re


def find_protein_hits(id_file, db_file, hit_file):
    ids = pd.read_csv(id_file)
    db = fasta.read(db_file)

    def remove_modification(row):
        return re.sub(r'\(.*?\)', '', row['Peptide'])

    ids['Pep'] = ids.apply(remove_modification, axis=1)
    ids.drop_duplicates(subset=['Pep'], keep='first', inplace=True)

    patterns = [re.escape(s) for s in ids['Pep']]
    patterns = re.compile('|'.join(patterns))

    hits = []
    for description, seq in db:
        if patterns.search(seq):
            hits.append(description)

    print('{0} protein hits found.'.format(len(hits)))

    with open(hit_file, 'w') as file:
        for hit in hits:
            file.write(hit + '\n')

import pandas as pd


def compare_peptide_lists(a, b):
    a, b = pd.read_csv(a), pd.read_csv(b)

    a.index = a['Index'].astype(str) + a['Scan'].astype(str) + a['Peptide']
    b.index = b['Index'].astype(str) + b['Scan'].astype(str) + b['Peptide']

    a_only = a[~a.index.isin(b.index)]
    b_only = b[~b.index.isin(a.index)]

    ab_common = pd.merge(a, b, how='inner', on=['Index', 'Scan', 'Peptide'])

    return a_only, b_only, ab_common

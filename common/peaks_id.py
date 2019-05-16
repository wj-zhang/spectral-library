"""
peaks_id - read peptide identifications from PEAKS pep.xml file into python data frame
"""

from pyteomics import pepxml


def read(source, keep=None, remove_duplicates=False):
    psms = pepxml.DataFrame(source)
    psms = psms[['index', '-10lgP', 'assumed_charge', 'start_scan', 'modified_peptide', 'peptide', 'spectrum', 'protein']]

    if remove_duplicates:
        psms.drop_duplicates(subset=['spectrum', 'start_scan'], keep=False, inplace=True)

    psms.sort_values(['-10lgP'], ascending=False, inplace=True)
    if keep == 'first':
        psms.drop_duplicates(subset=['modified_peptide', 'assumed_charge'], keep='first', inplace=True)

    return psms

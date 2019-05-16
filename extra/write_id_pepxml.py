import jinja2
import pandas as pd
from _datetime import datetime


class Parameters:
    def __init__(self):
        self.date = datetime.now().strftime("%Y-%m-%d %H:%M")


class Base:
    def __init__(self, source_file, df):
        self.source_file = source_file
        self.psms = (Psm(index, row) for index, row in df.iterrows())


class Psm:
    def __init__(self, index, row):
        self.spectrum = row['Source File']
        scan = row['Scan'].split(":")[1]
        self.start_scan = scan
        self.end_scan = scan
        self.precursor_neutral_mass = row['Mass']
        self.assumed_charge = int(round(row['Mass']/row['m/z']))
        self.index = index
        self.hit_rank = 1
        self.sequence = row['Peptide']
        try:
            proteins = row['Accession'].split(":")
        except AttributeError:
            proteins = ['']
        self.protein = proteins[0]
        self.alternative_proteins = proteins[1:]
        self.num_total_proteins = len(proteins)
        self.mass_diff = row['ppm']
        self.scores = {'-10lgP': str(row['-10lgP'])}


def write(id_file, pepxml_file, template_file = 'template.jinja'):
    template_loader = jinja2.FileSystemLoader('../extra/')
    template_env = jinja2.Environment(loader=template_loader)
    template = template_env.get_template(template_file)

    id_df_grouped = pd.read_csv(id_file).groupby(['Source File'])

    template_vars = {
        'path_to_output': pepxml_file,
        'parameters': Parameters(),
        'bases': (Base(source_file, id_df) for source_file, id_df in id_df_grouped)
    }

    with open(template_vars['path_to_output'], 'w') as output:
        output.write(template.render(template_vars))

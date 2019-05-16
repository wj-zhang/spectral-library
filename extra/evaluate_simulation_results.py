import pandas as pd


def evaluate(ture_id_file, search_id_file):
    # 0 mixture, 1 single, 2 true positive, 3 false positive, 4 false negative, 5 precision, 6 recall, 7 F-measure

    true_ids = pd.read_csv(ture_id_file)
    search_ids = pd.read_csv(search_id_file)
    ids = true_ids.join(search_ids.set_index('Index'), on='Index', rsuffix='_s')

    stat = {}
    stat['total'] = [0] * 8
    for idx, id in ids.iterrows():

            t_pep_pairs = list(zip(id['Peptide'].split(';'), id['Charge'].split(';')))
            if not pd.isnull(id['Peptide_s']):
                if not isinstance(id['Charge_s'], str):
                    charges = str(int(id['Charge_s']))
                else:
                    charges = id['Charge_s']
                s_pep_pairs = list(zip(id['Peptide_s'].split(';'), charges.split(';')))
            else:
                s_pep_pairs = []

            mode = id['Mode']
            if not mode in stat:

                stat[mode] = [0] * 8

            stat[mode][0] += 1
            stat[mode][1] += len(t_pep_pairs)

            stat[mode][2] += len(set(t_pep_pairs) & set(s_pep_pairs))
            stat[mode][3] += len(set(s_pep_pairs) - set(t_pep_pairs))
            stat[mode][4] += len(set(t_pep_pairs) - set(s_pep_pairs))

    for key in stat:
        if not key == 'total':
            compute_stat(stat, key)
            for i in range(5):
                stat['total'][i] += stat[key][i]
    key = 'total'
    compute_stat(stat, key)
    return stat


def compute_stat(stat, key):
    stat[key][5] = stat[key][2] / (stat[key][2] + stat[key][3])
    stat[key][6] = stat[key][2] / (stat[key][2] + stat[key][4])
    stat[key][7] = 2 * stat[key][5]* stat[key][6] / (stat[key][5] + stat[key][6])




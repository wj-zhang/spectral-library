import ms2
from sklearn.utils import shuffle
import math
from pyteomics import mgf, auxiliary
import peaks_modifications
import pandas as pd
import splib_modifications
import numpy as np
import random
from bisect import bisect_left


allowed_max_range_mz = 0
actual_max_range_mz = 2
removed_num = 0


def combine_identical_peaks(mz, intensities, tol=1.5e-5):
    idx = 0
    while idx < len(mz) - 1:
        idx_next = idx + 1
        if mz[idx_next] - mz[idx] < mz[idx] * tol:
            mz[idx] = (mz[idx]*intensities[idx] + mz[idx_next]*intensities[idx_next]) \
                      / (intensities[idx] + intensities[idx_next])
            intensities[idx] = intensities[idx] + intensities[idx_next]
            mz = np.delete(mz, idx_next)
            intensities = np.delete(intensities, idx_next)
        else:
            idx += 1

    return mz, intensities


def combine_multi_spectra(spectra, weight, mode=0):
    assert(len(spectra) == len(weight))

    ms = {
        'params': {
            'title': [], # string
            'pepmass': [], # float tuple
            'charge': auxiliary.structures.ChargeList([]), # int list
            'scans': [], # string
            'rtinseconds': [] # string
        },
        'm/z array': np.array([]),
        'intensity array': np.array([])
    }

    median_idx = (len(spectra) - 1) // 2
    spectra.insert(0, spectra.pop(median_idx))
    for i in range(len(spectra)):
        spectrum = spectra[i]
        ms['m/z array'] = np.append(ms['m/z array'], spectrum['m/z array'])
        ms['intensity array'] = np.append(ms['intensity array'], spectrum['intensity array']*weight[i])

        if not spectrum['params']['modified seq'] == '':
            ms['params']['title'].append(spectrum['params']['modified seq'])
        ms['params']['pepmass'].append(float(spectrum['params']['scan'][2]))  # string
        ms['params']['charge'].append(int(spectrum['params']['charge'][0])) # string
        ms['params']['scans'].append(spectrum['params']['scan'][0]) # string
        ms['params']['rtinseconds'].append(str(float(spectrum['params']['RTime'])*60)) # float

    sorted_idx = np.argsort(ms['m/z array'])
    ms['m/z array'] = ms['m/z array'][sorted_idx]
    ms['intensity array'] = ms['intensity array'][sorted_idx]
    ms['m/z array'], ms['intensity array'] = combine_identical_peaks(ms['m/z array'], ms['intensity array'])

    meta_info = ';' + str(len(weight)) + ' ' + str(mode)
    ms['params']['title'] = 'mod_seq=' + ' '.join(ms['params']['title']) + meta_info

    global actual_max_range_mz, removed_num
    # range_mz = (max(ms['params']['pepmass']) - min(ms['params']['pepmass'])) / 2
    mz_median = ms['params']['pepmass'][0]
    range_mz = max(max(ms['params']['pepmass']) - mz_median, mz_median - min(ms['params']['pepmass']))
    if range_mz >= allowed_max_range_mz:
        removed_num += 1
        return None
    if range_mz > actual_max_range_mz:
        actual_max_range_mz = range_mz

    # if len(spectra) > 1:
    #     ms['params']['pepmass'].insert(0, (max(ms['params']['pepmass']) + min(ms['params']['pepmass'])) / 2)
    ms['params']['pepmass'] = tuple(ms['params']['pepmass'])
    ms['params']['scans'] = ' '.join(ms['params']['scans'])
    ms['params']['rtinseconds'] = ' '.join(ms['params']['rtinseconds'])

    return ms


def construct_library_and_query(num_mixture, weights_list, base_spectra):
    library = []
    query = []

    if num_mixture == 0:
        for spectrum in base_spectra:
            query.append(combine_multi_spectra([spectrum], [1]))
        return library, query

    k, m = divmod(len(base_spectra) // num_mixture, len(weights_list))

    positions = np.array(list(range(len(base_spectra) // num_mixture))) * num_mixture
    positions = shuffle(positions, random_state=0)

    pos_i = 0
    for i in range(len(weights_list)):
        for j in range(k + min(i+1, m) - min(i, m)):
            pos = positions[pos_i]
            mixture_spectrum = combine_multi_spectra(base_spectra[pos:pos + num_mixture], weights_list[i], mode=i)
            if mixture_spectrum is not None:
                query.append(mixture_spectrum)
                library += [ms for ms in base_spectra[pos:pos + num_mixture] if not ms['params']['modified seq'] == '']
            pos_i += 1
    return library, query


def random_split(a_list, ratio, random_state=0):
    a_list = shuffle(a_list, random_state=random_state)
    a_list = shuffle(a_list, random_state=random_state)
    len_list = len(a_list)
    sum_pre_probs = 0
    splits = {}
    for key in ratio:
        idx_begin = math.floor(sum_pre_probs * len_list)
        sum_pre_probs += ratio[key]
        idx_end = math.floor(sum_pre_probs * len_list)
        splits[key] = a_list[idx_begin:idx_end]
    return splits


def preprocess_option(option):
    if 1 in option['ratio']:
        option['weights'][1] = [[1]]

    sum_ = sum([key*value for key, value in option['ratio'].items()])
    for key, value in option['ratio'].items():
        option['ratio'][key] = key * value / sum_

    # for key, value in option['weights'].items():
    #     for idx in range(len(value)):
    #         sum_ = sum(value[idx])
    #         option['weights'][key][idx] = [v / sum_ for v in value[idx]]

    global allowed_max_range_mz, actual_max_range_mz, removed_num
    allowed_max_range_mz = option['mz range']
    actual_max_range_mz = 0
    removed_num = 0

    return option


def generate(library_file, simulation_our_library_file, simulation_st_library_file, simulation_query_file, id_file, infor_file, option, random_state=0, is_shuffle=False, write=True, remove_peaks_ratio=0.3):
    random.seed(random_state)
    np.random.seed(random_state)

    option = preprocess_option(option)

    base_spectra_list = get_library_spectra(library_file)
    # base_spectra_list = randomly_remove_spectra_from_library(base_spectra_list)

    splits_library = random_split(base_spectra_list, option['ratio'], random_state=random_state)

    library_spectra = []
    query_spectra = []
    for num_mixture in option['ratio']:
        spectrum_list = splits_library[num_mixture]
        spectrum_list.sort(key=lambda x: float(x['params']['scan'][2]), reverse=False)
        library, query = construct_library_and_query(num_mixture, option['weights'][num_mixture], spectrum_list)
        library_spectra += library
        query_spectra += query

    if is_shuffle:
        query_spectra = shuffle(query_spectra, random_state=random_state)

    infor = infor_simulation(query_spectra, library_spectra, option)

    header = ms2.read_header(library_file)

    if write:
        write_ms2_file(simulation_our_library_file, header, library_spectra)
        library_spectra = convert_to_SectraST_library(library_spectra)
        write_ms2_file(simulation_st_library_file, header, library_spectra)

    # query_spectra = randomly_remove_peaks(query_spectra, ratio=remove_peaks_ratio)
    query_spectra = add_false_precursors_from_library_to_query(query_spectra, option, library_spectra, random_state, do=True)

    query_spectra, ids = post_process_query(query_spectra)

    if write:
        write_mgf_file(simulation_query_file, query_spectra)

        ids.to_csv(id_file, index=False)

        query_spectra = post_process_query_for_respect(query_spectra)
        write_mgf_file(simulation_query_file.replace('.mgf', '_respect.mgf'), query_spectra)

        with open(infor_file, 'w') as infor_write:
            infor_write.write(str(infor))

    return infor


def randomly_remove_spectra_from_library(library_spectra, ratio=0.05):
    idxs = shuffle(list(range(len(library_spectra))))
    num = math.floor(len(library_spectra) * ratio)
    idxs = idxs[:num]
    for idx in idxs:
        library_spectra[idx]['params']['modified seq'] = ''
    return library_spectra


def randomly_remove_peaks(query_spectra, ratio=0.3):

    for idx in range(len(query_spectra)):
        masses, intensities = query_spectra[idx]['m/z array'], query_spectra[idx]['intensity array']
        idxs = shuffle(list(range(len(masses))))
        num = math.floor(len(masses) * (1-ratio))
        idxs = idxs[:num]
        idxs.sort()
        query_spectra[idx]['m/z array'], query_spectra[idx]['intensity array'] = masses[idxs], intensities[idxs]

    return query_spectra


def write_ms2_file(filename, header, spectra):
    with open(filename, 'w') as destination:
        ms2.write_header(destination, header)
        for psm in spectra:
            ms2.write(destination, psm)


def write_mgf_file(filename, spectra):
    with open(filename, 'w') as destination:
        mgf.write(spectra, destination)


def convert_to_SectraST_library(spectra):
    for idx in range(len(spectra)):
        spectra[idx]['params']['modified seq'] = splib_modifications.convert_modification_backward(spectra[idx]['params']['modified seq'])
    return spectra


def get_library_spectra(library_file):
    spectra = list(ms2.read(library_file))
    for idx in range(len(spectra)):
        spectra[idx]['params']['scan'] = (str(idx + 1), str(idx + 1), spectra[idx]['params']['scan'][2])

    return spectra


def infor_simulation(query_spectra, library_spectra, option):
    infor = {}
    infor['library'] = len(library_spectra)
    infor['query'] = [len(query_spectra), 0]
    for key in option['weights']:
        infor[key] = [0]*len(option['weights'][key])

    for spectrum in query_spectra:
        num, mod = spectrum['params']['title'].split(';')[1].split()
        num, mod = int(num), int(mod)
        infor[num][mod] += 1
        infor['query'][1] += num

    infor['max mz range'] = actual_max_range_mz
    infor['removed spectra'] = removed_num
    return infor


def post_process_query_for_respect(query_spectra):
    for idx in range(len(query_spectra)):
        query_spectra[idx]['params']['rtinseconds'] = idx + 1
    return query_spectra


def post_process_query(query_spectra):
    ids = []
    for idx in range(len(query_spectra)):
        spectrum = query_spectra[idx]
        peptile_str, mode = spectrum['params']['title'][8:].split(';')
        peptides_list = [peaks_modifications.convert_modification_backward(peptide) for peptide in peptile_str.split(' ')]
        peptides_list = [peptide for peptide in peptides_list if not peptide == '']
        peptiles = ';'.join(peptides_list)
        mode = mode.replace(' ', ';')
        from_scans = spectrum['params']['scans'].replace(' ', ';')
        charges = ';'.join([str(charge) for charge in spectrum['params']['charge']][:len(peptides_list)])

        # if len(peptides_list) > 0:
        ids.append([idx+1, idx+1, peptiles, charges, from_scans, mode])

        spectrum['params']['title'] = 'index=' + str(idx)
        spectrum['params']['scans'] = str(idx + 1)
        spectrum['params']['rtinseconds'] = spectrum['params']['rtinseconds'].split(' ')[0]
        spectrum['intensity array'] = spectrum['intensity array'] * 10000

        query_spectra[idx] = spectrum

    ids = pd.DataFrame(ids, columns=['Index', 'Scan', 'Peptide', 'Charge', 'Library', 'Mode'])
    return query_spectra, ids


def binary_search(a, x):
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    else:
        return -1


def add_false_precursors_from_library_to_query(query_spectra, option, library_spectra, random_state, ratio=(0, 0, 0, 1), do=True):
    masses_list = []
    charges_list = []
    for spectrum in library_spectra:
        mz, z = float(spectrum['params']['scan'][2]), int(spectrum['params']['charge'][0])
        masses_list.append(mz)
        charges_list.append(z)

    sort_idx = np.argsort(masses_list)
    masses_list = [masses_list[idx] for idx in sort_idx]
    charges_list = [charges_list[idx] for idx in sort_idx]

    spectra_processed = []
    has_visited = [False]*len(query_spectra)
    for i in list(option['ratio'].keys()):
        spectra_i = []
        for idx in range(len(query_spectra)):
            if not has_visited[idx]:
                if len(query_spectra[idx]['params']['charge']) == i:
                    spectra_i.append(query_spectra[idx])
                    has_visited[idx] = True

        spectra_i = shuffle(spectra_i, random_state=random_state)

        spectra_i_processed = []
        points = [0] + [math.floor(sum(ratio[:i+1]) * len(spectra_i)) for i in range(len(ratio)-1)] + [len(spectra_i)]
        for num in range(len(ratio)):
            spectra_i_processed += add_false_precursors_from_library(spectra_i[points[num]:points[num+1]], num, library_spectra, sort_idx, masses_list, charges_list, do=do)

        spectra_processed += spectra_i_processed

    return spectra_processed


def add_false_precursors_from_library(spectra, num, library_spectra, sort_idx, masses_list, charges_list, weight=0.05, do=True):
    if num == 0 or not do:
        return spectra

    for idx in range(len(spectra)):
        mzs, zs = spectra[idx]['params']['pepmass'], spectra[idx]['params']['charge']
        mzs, zs = [mz for mz in mzs if mz is not None], [z.real for z in zs]
        for _ in range(num):
            while True:
                if random.uniform(0, 1) < 0.5:
                    mass = min(mzs)
                    pos = binary_search(masses_list, mass)
                    if pos == 0:
                        continue

                    mzs.append(masses_list[pos - 1])
                    zs.append(charges_list[pos - 1])
                    spectra[idx]['m/z array'] = np.append(spectra[idx]['m/z array'], library_spectra[sort_idx[pos - 1]]['m/z array'])
                    spectra[idx]['intensity array'] = np.append(spectra[idx]['intensity array'], library_spectra[sort_idx[pos - 1]]['intensity array'] * weight)
                    break
                else:
                    mass = max(mzs)
                    pos = binary_search(masses_list, mass)
                    if pos == len(masses_list) - 1:
                        continue

                    mzs.append(masses_list[pos + 1])
                    zs.append(charges_list[pos + 1])

                    spectra[idx]['m/z array'] = np.append(spectra[idx]['m/z array'], library_spectra[sort_idx[pos + 1]]['m/z array'])
                    spectra[idx]['intensity array'] = np.append(spectra[idx]['intensity array'], library_spectra[sort_idx[pos + 1]]['intensity array'] * weight)
                    break

        spectra[idx]['params']['pepmass'] = tuple(mzs)
        spectra[idx]['params']['charge'] = auxiliary.structures.ChargeList([])
        for z in zs:
            spectra[idx]['params']['charge'].append(z)

        sorted_idx_peaks = np.argsort(spectra[idx]['m/z array'])
        spectra[idx]['m/z array'] = spectra[idx]['m/z array'][sorted_idx_peaks]
        spectra[idx]['intensity array'] = spectra[idx]['intensity array'][sorted_idx_peaks]
        spectra[idx]['m/z array'], spectra[idx]['intensity array'] = combine_identical_peaks(spectra[idx]['m/z array'], spectra[idx]['intensity array'])

    return spectra


def add_false_precursors(spectra, num, masses_list, charges_list):
    if num == 0:
        return spectra

    for idx in range(len(spectra)):
        mzs, zs = spectra[idx]['params']['pepmass'], spectra[idx]['params']['charge']
        mzs, zs = [mz for mz in mzs if mz is not None], [z.real for z in zs]
        for _ in range(num):
            while True:
                if random.uniform(0, 1) < 0.5:
                    mass = min(mzs)
                    pos = binary_search(masses_list, mass)
                    if pos == 0:
                        continue

                    mzs.append(masses_list[pos - 1])
                    zs.append(charges_list[pos - 1])
                    break
                else:
                    mass = max(mzs)
                    pos = binary_search(masses_list, mass)
                    if pos == len(masses_list) - 1:
                        continue

                    mzs.append(masses_list[pos + 1])
                    zs.append(charges_list[pos + 1])
                    break

        spectra[idx]['params']['pepmass'] = tuple(mzs)
        spectra[idx]['params']['charge'] = auxiliary.structures.ChargeList([])
        for z in zs:
            spectra[idx]['params']['charge'].append(z)
    return spectra


def add_false_precursors_to_query(query_spectra, option, random_state, ratio=(0.4, 0.3, 0.2, 0.1)):
    masses_list = []
    charges_list = []
    for spectrum in query_spectra:
        mzs, zs = spectrum['params']['pepmass'], spectrum['params']['charge']
        mzs, zs = [mz for mz in mzs if mz is not None], [z.real for z in zs]
        masses_list += mzs
        charges_list += zs

    sort_idx = np.argsort(masses_list)
    masses_list = [masses_list[idx] for idx in sort_idx]
    charges_list = [charges_list[idx] for idx in charges_list]

    spectra_processed = []
    has_visited = [False]*len(query_spectra)
    for i in list(option['ratio'].keys()):
        spectra_i = []
        for idx in range(len(query_spectra)):
            if not has_visited[idx]:
                if len(query_spectra[idx]['params']['charge']) == i:
                    spectra_i.append(query_spectra[idx])
                    has_visited[idx] = True

        spectra_i = shuffle(spectra_i, random_state=random_state)

        spectra_i_processed = []
        points = [0] + [math.floor(sum(ratio[:i+1]) * len(spectra_i)) for i in range(len(ratio)-1)] + [len(spectra_i)]
        for num in range(len(ratio)):
            spectra_i_processed += add_false_precursors(spectra_i[points[num]:points[num+1]], num, masses_list, charges_list)

        spectra_processed += spectra_i_processed

    return spectra_processed


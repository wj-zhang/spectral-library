import numpy as np
from spectrum_snr import spectrum_snr
from dot_product import dot_product
from copy import deepcopy
from spectrum_alignment import spectrum_alignment


def create(ms2_list, ion_ppm):
    if len(ms2_list) == 1:
        ms2_list[0]['params']['NReplicates'] = 1
        return ms2_list[0]

    snr_list = sort_spectra_by_quality(ms2_list)
    remove_dissimilar_replicates(ms2_list, ion_ppm, snr_list)
    spectrum_aligned = align_spectra(ms2_list, snr_list, ion_ppm)
    remove_noise_peaks(spectrum_aligned, len(ms2_list))
    return spectrum_aligned


def sort_spectra_by_quality(ms2_list, inplace=True):
    snr_list = [spectrum_snr(ms2['intensity array']) for ms2 in ms2_list]
    indices = sorted(range(len(snr_list)), key=snr_list.__getitem__, reverse=True)
    snr_list = [snr_list[i] for i in indices]
    ms2_list[:] = [ms2_list[i] for i in indices]
    return snr_list


def remove_dissimilar_replicates(ms2_list, ion_ppm, snr_list, inplace=True):
    index_pair, first_ms2 = (0, 0), None
    for idx, ms2 in enumerate(ms2_list):
        if first_ms2:
            if dot_product(first_ms2, ms2, ion_ppm=ion_ppm) < 0.6:
                if idx - i > index_pair[1] - index_pair[0]:
                    index_pair = (i, idx)
                i, first_ms2 = idx, ms2
        else:
            i, first_ms2 = idx, ms2

    idx += 1
    if idx - i > index_pair[1] - index_pair[0]:
        index_pair = (i, idx)

    ms2_list[:] = ms2_list[index_pair[0]:index_pair[1]]
    snr_list[:] = snr_list[index_pair[0]:index_pair[1]]


def align_spectra(ms2_list, snr_list, ion_ppm):
    ms_merged = None
    for idx, ms2 in enumerate(ms2_list):
        if ms_merged is None:
            ms_merged = deepcopy(ms2)
            ms_merged['count array'] = np.ones(ms2['m/z array'].shape)
            w_merged = snr_list[idx]
            continue

        w = snr_list[idx]
        w_sum = w_merged + w
        v_merged, v = w_merged/w_sum, w/w_sum
        w_merged = w_sum/2

        alignment = spectrum_alignment(ms_merged, ms2, ion_ppm=ion_ppm)
        for i, j in alignment.items():
            ms_merged['m/z array'][i] = v_merged * ms_merged['m/z array'][i] + v * ms2['m/z array'][j]
            ms_merged['intensity array'][i] = v_merged * ms_merged['intensity array'][i] + v * ms2['intensity array'][j]
            ms_merged['count array'][i] += 1

        indices = list(alignment.values())
        ms2_mz = np.delete(ms2['m/z array'], indices)
        ms2_intensity = np.delete(ms2['intensity array'], indices)
        ms2_count = np.ones(ms2_mz.shape)

        ms_merged['m/z array'] = np.concatenate((ms_merged['m/z array'], ms2_mz))
        ms_merged['intensity array'] = np.concatenate((ms_merged['intensity array'], ms2_intensity))
        ms_merged['count array'] = np.concatenate((ms_merged['count array'], ms2_count))

        indices = np.argsort(ms_merged['m/z array'])
        ms_merged['m/z array'] = ms_merged['m/z array'][indices]
        ms_merged['intensity array'] = ms_merged['intensity array'][indices]
        ms_merged['count array'] = ms_merged['count array'][indices]

    ms_merged['params']['NReplicates'] = len(ms2_list)
    return ms_merged


def remove_noise_peaks(spectrum_aligned, spectrum_count, inplace=True):
    count = int((1+spectrum_count)*0.6) -1
    larger_than_count = spectrum_aligned['count array'] > count

    spectrum_aligned['m/z array'] = spectrum_aligned['m/z array'][larger_than_count]
    spectrum_aligned['intensity array'] = spectrum_aligned['intensity array'][larger_than_count]

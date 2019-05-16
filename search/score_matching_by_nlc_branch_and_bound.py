import numpy as np
import nonnegative_linear_combination
from spectrum_alignment import spectrum_alignment
import copy


spectra_unique_global = []
vectors_unique_global = []
scores_unique_global = []
ids_list_global = []
ids_list_copy_global = []
score_global = 0
weights_global = []
hits_global = []
query_global = []
remaining_bounds_global = []


def scoring(query, spectra_list, ion_ppm):
    ids_list = []
    spectra_unique = []
    vectors_unique = []
    scores_no_dups = []
    for spectra in spectra_list:
        ids = []
        for spectrum in spectra:
            try:
                id = spectra_unique.index(spectrum)
            except ValueError:
                alignment = spectrum_alignment(query, spectrum, ion_ppm=ion_ppm)
                if len(alignment) < 5:
                    continue
                vector = construct_vector(query['intensity array'], spectrum['intensity array'], list(alignment.keys()), list(alignment.values()))
                score = np.dot(query['intensity array'], vector[:-1])

                id = len(spectra_unique)
                spectra_unique.append(spectrum)
                vectors_unique.append(vector)
                scores_no_dups.append(score)
            ids.append(id)
        if ids:
            ids_list.append(ids)

    if len(ids_list) == 0:
        return None, None, None, None

    global spectra_unique_global, vectors_unique_global, scores_unique_global, score_global, weights_global, hits_global, ids_list_global, query_global, remaining_bounds_global, ids_list_copy_global, _MAX_NUM

    for i in range(len(ids_list)):
        id_scores = zip(ids_list[i], [scores_no_dups[id] for id in ids_list[i]])
        id_scores = sorted(id_scores, key=lambda x: x[1], reverse=True)
        ids_list[i] = list(list(zip(*id_scores))[0])

    ids_list = sorted(ids_list, key=lambda x: scores_no_dups[x[0]], reverse=True)
    ids_list_copy_global = copy.deepcopy(ids_list)

    max_num = 3

    for i in range(len(ids_list)):
        ids_list[i] = ids_list[i][:min(max_num, len(ids_list[i]))]

    spectra_unique_global = spectra_unique
    vectors_unique_global = vectors_unique
    scores_unique_global = scores_no_dups
    score_global = 0
    weights_global = []
    hits_global = []
    ids_list_global = ids_list
    query_global = query
    remaining_bounds_global = [0]

    for i in reversed(range(len(ids_list))):
        remaining_bounds_global.append(remaining_bounds_global[-1] + scores_no_dups[ids_list[i][0]] ** 2)
    remaining_bounds_global.reverse()
    del remaining_bounds_global[0]

    pre_bound = 0
    pre_hits = []
    now_cand_list_idx = 0
    branch_and_bound(pre_hits, pre_bound, now_cand_list_idx)

    hit_ids = [id for id in hits_global if id is not None]
    weights = [weight for weight in weights_global if weight is not None]
    if hit_ids:
        return [spectra_unique[id] for id in hit_ids], np.sqrt(score_global), weights, [np.linalg.norm(vectors_unique[id][:-1]) for id in hit_ids]
    else:
        return None, None, None, None


def branch_and_bound(pre_hits, pre_bound, now_cand_list_idx):
    if now_cand_list_idx == len(ids_list_global):
        C, d = get_C_ids(pre_hits, vectors_unique_global, query_global)
        w = nonnegative_linear_combination.nlc(C, d)
        score = 1 - np.linalg.norm(C.dot(w).squeeze() - d)**2
        global score_global, hits_global, weights_global
        if score > score_global:
            score_global = score
            hits_global = copy.copy(pre_hits)
            weights_global = copy.copy(w)
        return

    hit = False
    increas_num = 0
    for id in ids_list_global[now_cand_list_idx]:
        if id in pre_hits:
            if len(ids_list_global[now_cand_list_idx]) + 1 <= len(ids_list_copy_global[now_cand_list_idx]):
                increas_num += 1
                ids_list_global[now_cand_list_idx].append(ids_list_copy_global[now_cand_list_idx][len(ids_list_global[now_cand_list_idx])])
            continue
        hit = True
        pre_bound_now = pre_bound + scores_unique_global[id] ** 2
        bound_now = pre_bound_now + remaining_bounds_global[now_cand_list_idx]
        if bound_now <= score_global:
            return
        pre_hits.append(id)
        branch_and_bound(pre_hits, pre_bound_now, now_cand_list_idx + 1)
        del pre_hits[-1]

    for i in range(increas_num):
        del ids_list_global[now_cand_list_idx][-1]

    if not hit:
        pre_bound_now = pre_bound
        bound_now = pre_bound_now + remaining_bounds_global[now_cand_list_idx]
        if bound_now <= score_global:
            return
        pre_hits.append(None)
        branch_and_bound(pre_hits, pre_bound_now, now_cand_list_idx + 1)
        del pre_hits[-1]


def construct_vector(a, b, idx_a, idx_b):
    vector = np.zeros(len(a) + 1)
    vector[idx_a] = b[idx_b]
    vector[-1] = np.sqrt(max(np.linalg.norm(b) ** 2 - np.linalg.norm(vector) ** 2, 0))
    return vector


def get_C_ids(hit_ids, vectors_no_dups, query):
    ids_non_none = [id for id in hit_ids if id is not None]
    C = []
    num = 0
    for id in ids_non_none:
        a = np.append(vectors_no_dups[id][:-1], [0]*len(ids_non_none))
        a[len(vectors_no_dups[id]) - 1 + num] = vectors_no_dups[id][-1]
        num += 1
        C.append(a)
    return np.array(C).transpose(), np.append(query['intensity array'], [0]*len(ids_non_none))

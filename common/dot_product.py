"""
dot_product - calculate the dot product between two arrays
"""

import numpy as np
from sklearn.preprocessing import normalize


def dot_product(*args, **kwds):
    try:
        ion_ppm = kwds['ion_ppm']
    except KeyError:
        ion_ppm = 0.00002

    if len(args) == 2 and isinstance(args[0], dict):
        masses_a, intensities_a = list(args[0]['m/z array']), list(args[0]['intensity array'])
        masses_b, intensities_b = list(args[1]['m/z array']), list(args[1]['intensity array'])
    elif len(args) == 4 and isinstance(args[0], (list, np.ndarray)):
        masses_a, intensities_a = list(args[0]), list(args[1])
        masses_b, intensities_b = list(args[2]), list(args[3])
    else:
        raise ValueError

    intensities_a = normalize_list(intensities_a)
    intensities_b = normalize_list(intensities_b)

    return aligned_dot_product(masses_a, intensities_a, masses_b, intensities_b, ion_ppm)


def aligned_dot_product(masses_a, intensities_a, masses_b, intensities_b, ion_ppm):
    a_iter = iter(zip(masses_a, intensities_a))
    b_iter = iter(zip(masses_b, intensities_b))

    def find_first_valid_pair(a, b):
        try:
            while True:
                error = a[0] * ion_ppm
                if b[0]-a[0] > error:
                    a = next(a_iter)
                elif a[0]-b[0] > error:
                    b = next(b_iter)
                else:
                    a_array.append(a)
                    b_array.append(b)
                    return a, b
        except StopIteration:
            return None, None

    def extend_to_right(a, b):
        a_next, b_next = None, None
        a_flag, b_flag = True, True

        while a_flag or b_flag:
            a_flag, b_flag = False, False

            try:
                while True:
                    error = a[0] * ion_ppm
                    b_next = next(b_iter) if b_next is None else b_next
                    if b_next[0]-a[0] > error:
                        break
                    else:
                        b = b_next
                        b_array.append(b)
                        b_next = None
                        b_flag = True
            except StopIteration:
                pass

            try:
                while True:
                    error = b[0] * ion_ppm
                    a_next = next(a_iter) if a_next is None else a_next
                    if a_next[0]-b[0] > error:
                        break
                    else:
                        a = a_next
                        a_array.append(a)
                        a_next = None
                        a_flag = True
            except StopIteration:
                pass

        return a_next, b_next

    def dynamic_programming():
        len_a, len_b = len(a_array), len(b_array)
        matrix = [[0] * (len_b + 1) for i in [0, 1]]

        k = 0
        for i in range(len_a):
            error = a_array[i][0] * ion_ppm
            for j in range(len_b):
                intensity = a_array[i][1] * b_array[j][1] if abs(a_array[i][0] - b_array[j][0]) < error else 0
                matrix[k][j + 1] = max(matrix[k][j], matrix[1 - k][j + 1], matrix[1 - k][j] + intensity)
            k = 1 - k

        return matrix[1 - k][len_b]

    def cal_product():
        if len(a_array) == 1 and len(b_array) == 1:
            return a_array[0][1] * b_array[0][1]
        if len(a_array) > 1 or len(b_array) > 1:
            return dynamic_programming()

    try:
        prod = 0
        a = next(a_iter)
        b = next(b_iter)

        while a is not None and b is not None:
            a_array, b_array = [], []
            a, b = find_first_valid_pair(a, b)
            if a is None:
                break
            a, b = extend_to_right(a, b)
            prod += cal_product()
    except StopIteration:
        print("Iter stopped with error.")

    return prod


def normalize_list(x):
    x = np.asarray(x).reshape(-1, 1)
    x = normalize(x, axis=0)
    return x.reshape(-1, )

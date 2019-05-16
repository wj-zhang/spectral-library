import numpy as np


def spectrum_alignment(*args, **kwds):
    try:
        ion_ppm = kwds['ion_ppm']
    except KeyError:
        ion_ppm = 0.00002

    if len(args) == 2 and isinstance(args[0], dict):
        return align(args[0]['m/z array'], args[0]['intensity array'],
                                   args[1]['m/z array'], args[1]['intensity array'], ion_ppm)
    elif len(args) == 4 and isinstance(args[0], (list, np.ndarray)):
        return align(args[0], args[1], args[2], args[3], ion_ppm)
    else:
        raise ValueError


def align(masses_a, intensities_a, masses_b, intensities_b, ion_ppm):
    alignment = {}
    la, ra, lb, rb = 0, 0, 0, 0

    def find_first_valid_pair(la, lb):
        try:
            while True:
                error = masses_a[la] * ion_ppm
                if masses_b[lb]-masses_a[la] > error:
                    la += 1
                elif masses_a[la]-masses_b[lb] > error:
                    lb += 1
                else:
                    return la, lb
        except IndexError:
            return None, None

    def extend_to_right(ra, rb):
        a_flag, b_flag = True, True

        while a_flag or b_flag:
            a_flag, b_flag = False, False

            try:
                error = masses_a[ra] * ion_ppm
                while True:
                    rb += 1
                    if masses_b[rb]-masses_a[ra] > error:
                        rb -= 1
                        break
                    else:
                        b_flag = True
            except IndexError:
                rb -= 1

            try:
                error = masses_b[rb] * ion_ppm
                while True:
                    ra += 1
                    if masses_a[ra]-masses_b[rb] > error:
                        ra -= 1
                        break
                    else:
                        a_flag = True
            except IndexError:
                ra -= 1

        return ra, rb

    def dynamic_programming(la, ra, lb, rb):
        len_a, len_b = ra-la+1, rb-lb+1
        traces, products = [[0] * (len_b + 1) for i in range(len_a+1)], [[0] * (len_b + 1) for i in range(len_a+1)]

        for i in range(len_a):
            error = masses_a[la + i] * ion_ppm
            for j in range(len_b):
                if abs(masses_a[la + i] - masses_b[lb + j]) < error:
                    intensity = intensities_a[la + i] * intensities_b[lb + j]
                else:
                    intensity = 0
                u = [products[i+1][j], products[i][j + 1], products[i][j] + intensity]
                traces[i+1][j+1] = np.argmax(u)
                products[i+1][j+1] = u[traces[i+1][j+1]]

        i, j = len_a, len_b
        v = [(0, -1), (-1, 0), (-1, -1)]

        while i > 0 and j > 0:
            w = traces[i][j]
            if w == 2:
                alignment[la+i-1] = lb+j-1
            i, j = i+v[w][0], j+v[w][1]

    def cal_product(la, ra, lb, rb):
        if la==ra and lb==rb:
            alignment[la] = lb
        else:
            dynamic_programming(la, ra, lb, rb)

    len_a, len_b = len(masses_a), len(masses_b)
    while la < len_a and lb < len_b:
        la, lb = find_first_valid_pair(la, lb)
        if la is None:
            break
        ra, rb = la, lb
        ra, rb = extend_to_right(ra, rb)
        cal_product(la, ra, lb, rb)
        la, lb = ra+1, rb+1

    return alignment

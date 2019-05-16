"""
ms2 - read and write MS/MS data in MS2 format
"""

import numpy as np


def read(source):
    params, masses, intensities = {}, [], []

    with open(source, 'r') as source:
        for line in source:
            s = line.strip().split('\t')
            if not s or s[0] == 'H':
                pass
            elif s[0] == 'S':
                if masses:
                    yield {'params': params, 'm/z array': np.array(masses), 'intensity array': np.array(intensities)}
                    params, masses, intensities = {}, [], []
                params['scan'] = tuple(s[1:])
            elif s[0] == 'I':
                params[s[1]] = s[2]
            elif s[0] == 'Z':
                params['charge'] = tuple(s[1:])
            elif s[0] == 'D':
                params[s[1]] = s[2]
            elif s[0] == 'N':
                params[s[1]] = int(s[2])
            elif s[0] == 'M':
                params[s[1]] = float(s[2])
            elif s[0] == 'L':
                params[s[1]] = s[2]
            else:
                masses.append(float(s[0]))
                intensities.append(float(s[1]))

    yield {'params': params, 'm/z array': np.array(masses), 'intensity array': np.array(intensities)}


def read_header(source):
    with open(source, 'r') as source:
        header = {}
        for line in source:
            if line[0] != 'H':
                break
            l = line.split('\t')
            header[l[1]] = l[2].strip()

        return header


def write(destination, psm):
    params = psm['params']
    masses = ['{:.4f}'.format(x) for x in psm['m/z array']]
    intensities = ['{:.4f}'.format(x) for x in psm['intensity array']]

    params_str = 'S\t' + '\t'.join(map(str, params['scan'])) + '\n' \
                 + 'I\tRTime\t' + str(params['RTime']) + '\n' \
                 + 'Z\t' + '\t'.join(map(str, params['charge'])) + '\n' \
                 + 'D\tseq\t' + params['seq'] + '\n' \
                 + 'D\tmodified seq\t' + params['modified seq'] + '\n'

    variable_dict = {
        'NReplicates': 'N',
        'Protein': 'L',
        'iRT': 'M',
        'CCS': 'M'
    }

    for k, v in variable_dict.items():
        try:
            value = params[k]
            params_str += v + '\t' + k + '\t' + str(value) + '\n'
        except KeyError:
            pass

    destination.write(params_str)

    for mass, intensity in zip(masses, intensities):
        destination.write('\t'.join([mass, intensity]) + '\n')


def write_header(destination, header):
    for key, value in header.items():
        destination.write('H\t' + '\t'.join([key, value]) + '\n')

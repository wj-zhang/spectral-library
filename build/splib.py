"""
splib - read spectral library in splib format
"""

import re
import mmap
import struct
import numpy as np


def read(splib_file):
    file = open(splib_file, 'rb')
    mm = mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ)

    _read_header(mm)

    try:
        while True:
            yield _read_spectrum(mm)
    except StopIteration:
        pass

    mm.close()
    file.close()


_annotation_ion_types = frozenset(b'abcxyz')
_ignore_annotations = False


def _parse_annotation(raw):
    if raw[0] not in _annotation_ion_types:
        return None

    first_annotation = raw.split(b',', 1)[0]
    if b'i' in first_annotation:
        return None

    ion_sep = first_annotation.find(b'/')
    if ion_sep == -1:
        ion_sep = len(first_annotation)

    first_annotation_substring = first_annotation[:ion_sep]

    if b'-' in first_annotation_substring or b'+' in first_annotation_substring:
        return None

    charge_sep = first_annotation.find(b'^')
    if charge_sep != -1:
        ion_type = first_annotation[:charge_sep].decode('UTF-8')
        charge = int(first_annotation[charge_sep + 1: ion_sep])
    else:
        ion_type = first_annotation[:ion_sep].decode('UTF-8')
        charge = 1

    return ion_type, charge


def _parse_retention_time(comment):
    comments = comment.decode(encoding='UTF-8').split(' ')
    rts = [x[14:] for x in comments if 'RetentionTime' in x]
    if len(rts) != 1:
        return 0
    rts = [float(rt) for rt in rts[0].split(',')]

    return sum(rts) / len(rts)


def _read_header(mm):
    version = struct.unpack('i', mm.read(4))[0]
    sub_version = struct.unpack('i', mm.read(4))[0]
    file_name = mm.readline()
    num_lines = struct.unpack('i', mm.read(4))[0]
    comments = [mm.readline() for i in range(num_lines)]

    return version, sub_version, file_name, num_lines, comments


def _read_spectrum(mm):
    read_bytes = mm.read(4)
    if not read_bytes:
        raise StopIteration

    identifier = struct.unpack('i', read_bytes)[0]
    name = mm.readline().strip()
    modified_seq = name[name.find(b'.')+1 : name.rfind(b'.')].decode(encoding='UTF-8')
    seq = re.sub(r'\[.*?\]', '', modified_seq)

    slash_sep = name.rfind(b'/')
    z = int(name[slash_sep+1 : slash_sep+2])
    mz = struct.unpack('d', mm.read(8))[0]
    status = mm.readline().strip()

    num_peaks = struct.unpack('i', mm.read(4))[0]
    masses = np.empty((num_peaks,), np.float64)
    intensities = np.empty((num_peaks,), np.float64)
    annotations = np.empty((num_peaks,), object)
    for i in range(num_peaks):
        masses[i] = struct.unpack('d', mm.read(8))[0]
        intensities[i] = struct.unpack('d', mm.read(8))[0]
        annotation = mm.readline().strip()
        if not _ignore_annotations:
            annotations[i] = _parse_annotation(annotation)
        info = mm.readline()

    comment = mm.readline()
    rt = _parse_retention_time(comment)
    if b' Remark=DECOY_' in comment:
        seq = '#' + seq

    params = {'scan': tuple([identifier, identifier, mz]),
              'charge': tuple([z, mz]),
              'RTime': rt,
              'seq': seq,
              'modified seq': modified_seq
              }

    return {'params': params, 'm/z array': masses, 'intensity array': intensities}

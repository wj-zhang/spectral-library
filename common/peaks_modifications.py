import re
from modifications import modifications_env

peaks_modifications = {
    '(+57.02)': 4,
    '(+15.99)': 35,
    '(+.98)': 7,
    '(+42.01)': 1,
    '(+40.02)': 26,
    '(-18.01)': 27,
    '(-17.02)': 28
}

mod_env = modifications_env(peaks_modifications)

convert_modification_list_forward = ['.{}'.format(re.escape(k)) for k in mod_env['originals_to_symbols'].keys()]
convert_modification_regex_forward = re.compile('|'.join(convert_modification_list_forward))


def _convert_modification_forward(m):
    modification = mod_env['originals_to_symbols'][m.string[m.start() + 1:m.end()]]
    if m.string[m.start()].isupper():
        modification += m.string[m.start()]
    return modification


def convert_modification_forward(peptide):
    return convert_modification_regex_forward.sub(_convert_modification_forward, peptide)


convert_modification_list_backward = ['{}.'.format(re.escape(k)) for k in mod_env['symbols_to_originals'].keys()]
convert_modification_regex_backward = re.compile('|'.join(convert_modification_list_backward))


def _convert_modification_backward(m):
    return m.string[m.end()-1] + mod_env['symbols_to_originals'][m.string[m.start():m.end()-1]]


def convert_modification_backward(peptide):
    return convert_modification_regex_backward.sub(_convert_modification_backward, peptide)

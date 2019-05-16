import re
from modifications import modifications_env


splib_modifications = {
    'C[160]': 4,
    'M[147]': 35,
    'Q[111]': 28,
    'E[111]': 27,
    'C[143]': 26,
    'n[43]': 1
}

mod_env = modifications_env(splib_modifications)

convert_modification_list_forward = [re.escape(k) for k in mod_env['originals_to_symbols'].keys()]
convert_modification_regex_forward = re.compile('|'.join(convert_modification_list_forward))


def _convert_modification_forward(m):
    modification = mod_env['originals_to_symbols'][m.string[m.start():m.end()]]
    if m.string[m.start()].isupper():
        modification += m.string[m.start()]
    return modification


def convert_modification_forward(peptide):
    return convert_modification_regex_forward.sub(_convert_modification_forward, peptide)


convert_modification_list_backward = ['{}[A-Z]'.format(re.escape(k)) for k in mod_env['symbols_to_originals'].keys()]
convert_modification_regex_backward = re.compile('|'.join(convert_modification_list_backward))


def _convert_modification_backward(m):
    modification = mod_env['symbols_to_originals'][m.string[m.start():m.end() - 1]]
    if modification.islower():
        modification += m.string[m.end()-1]
    return modification


def convert_modification_backward(peptide):
    return convert_modification_regex_backward.sub(_convert_modification_backward, peptide)

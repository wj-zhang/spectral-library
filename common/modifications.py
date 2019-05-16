import re
from pyteomics import mass
from pyteomics.parser import std_labels


unimod_db = mass.Unimod(source=open('../data/unimod.xml', 'r'))


def convert_between_int_and_str(x):
    if isinstance(x, int):
        y = [chr(ord('a') + int(i)) for i in str(x)]
        y = ''.join(y)
    else:
        y = [str(ord(i)-ord('a')) for i in x]
        y = int(''.join(y))
    return y


def get_modification_by_record_id(id):
    for mod in unimod_db.mods:
        if mod['record_id'] == id:
            return mod


def modifications_env(mods):
    if isinstance(mods, dict):
        env = modifications_env(list(mods.values()))
        originals_to_symbols = {o: env['digits_to_symbols'][d] for o, d in mods.items()}
        symbols_to_originals = {s: o for o, s in originals_to_symbols.items()}
        env['originals_to_symbols'] = originals_to_symbols
        env['symbols_to_originals'] = symbols_to_originals
        return env

    if all(isinstance(m, int) for m in mods):
        digits = mods
        symbols = [convert_between_int_and_str(m) for m in mods]
    elif all(isinstance(m, str) for m in mods):
        symbols = mods
        digits = [convert_between_int_and_str(m) for m in mods]

    digits_to_symbols = {d:s for (d,s) in zip(digits, symbols)}
    symbols_to_digits = {s:d for d,s in digits_to_symbols.items()}

    digits_to_symbols_str = {str(d):s for d,s in digits_to_symbols.items()}
    symbols_to_digits_str = {s:d for d,s in digits_to_symbols_str.items()}

    aa_mass = mass.std_aa_mass
    for d, s in digits_to_symbols.items():
        mod = get_modification_by_record_id(d)
        aa_mass[s] = mod['mono_mass']

    aa_labels = std_labels + symbols

    env = {
        'digits': digits,
        'symbols': symbols,
        'digits_to_symbols': digits_to_symbols,
        'symbols_to_digits': symbols_to_digits,
        'digits_to_symbols_str': digits_to_symbols_str,
        'symbols_to_digits_str': symbols_to_digits_str,
        'aa_mass': aa_mass,
        'aa_labels': aa_labels
    }

    return env

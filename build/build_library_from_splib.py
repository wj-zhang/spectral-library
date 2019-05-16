import splib
import datetime
from preprocessing import preprocess
from splib_modifications import mod_env, convert_modification_forward


def build(splib_file):
    reader = splib.read(splib_file)

    for psm in reader:
        modified_seq = convert_modification_forward(psm['params']['modified seq'])
        psm['params']['modified seq'] = modified_seq
        if '[' in modified_seq:
            print('Unknown modifications: %s' % modified_seq)
            raise KeyError

        masses, intensities = preprocess(psm, peptide=modified_seq, mod_env=mod_env)

        if masses is None:
            continue

        yield {'params': psm['params'], 'm/z array': masses, 'intensity array': intensities}


def build_header(splib_file):
    header = {
        'CreationDate': datetime.datetime.now().strftime("%c"),
        'CreationFile': splib_file,
        'Modifications': str(mod_env['digits'])
    }

    return header

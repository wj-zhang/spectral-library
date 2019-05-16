from pyteomics import mgf


def correct(spectra, precursor_mass_shift=None, ion_mass_shift=None, mode='library'):
    precursor_mass_shift, ion_mass_shift = precursor_mass_shift/1000000.0, ion_mass_shift/1000000.0

    if isinstance(spectra, str):
        if mode is 'query':
            with mgf.read(spectra) as reader:
                spectra = [correcting(spectrum, precursor_mass_shift, ion_mass_shift) for spectrum in reader]
        elif mode is 'library':
            with mgf.read(spectra) as reader:
                spectra = [correcting(spectrum, precursor_mass_shift, ion_mass_shift) for spectrum in reader if spectrum['params']['sequence']]
    else:
        if mode is 'query':
            spectra = [correcting(spectrum, precursor_mass_shift, ion_mass_shift) for spectrum in spectra]
        elif mode is 'library':
            spectra = [correcting(spectrum, precursor_mass_shift, ion_mass_shift) for spectrum in spectra if spectrum['params']['sequence']]

    return spectra


def correcting(spectrum, precursor_mass_shift, ion_mass_shift):
    if precursor_mass_shift:
        mzs = []
        for mz in spectrum['params']['pepmass']:
            if mz:
                mz -= mz * precursor_mass_shift
                mzs.append(mz)
        spectrum['params']['pepmass'] = tuple(mzs)

    if ion_mass_shift:
        mzs = spectrum['m/z array']
        spectrum['m/z array'] = [mz-mz*ion_mass_shift for mz in mzs]

    return spectrum


def correct_from_pepxml(spectra, index, precursor_mass_shift=None, ion_mass_shift=None):
    precursor_mass_shift, ion_mass_shift = precursor_mass_shift / 1000000.0, ion_mass_shift / 1000000.0
    return [correcting(spectra[idx], precursor_mass_shift, ion_mass_shift) for idx in index]

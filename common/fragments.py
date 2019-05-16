from pyteomics import mass


def ion_masses(peptide, aa_mass=mass.std_aa_mass, types=('b', 'y'), charge=1):
    ion_pos = [i+1 for i, c in enumerate(peptide) if c.isupper()]
    ion_pos = ion_pos[:-1]
    b_ion_masses = [mass.fast_mass2(peptide[:i], ion_type=types[0], charge=charge, aa_mass=aa_mass) for i in ion_pos]
    y_ion_masses = [mass.fast_mass2(peptide[i:], ion_type=types[1], charge=charge, aa_mass=aa_mass) for i in ion_pos]

    return b_ion_masses, y_ion_masses[::-1]


def ion_annotation(masses, intensities, peptide, aa_mass=mass.std_aa_mass, ion_ppm=0.00002, types=('b', 'y'), charge=1):
    if 'b' not in types or 'y' not in types:
        raise TypeError

    b_ion_masses, y_ion_masses = ion_masses(peptide, aa_mass=aa_mass, types=types, charge=charge)

    b_ion_flags = single_type_ion_annotation(masses, intensities, b_ion_masses, ion_ppm)
    y_ion_flags = single_type_ion_annotation(masses, intensities, y_ion_masses, ion_ppm)

    ion_flags = merge(b_ion_flags, y_ion_flags)

    return ion_flags


def single_type_ion_annotation(masses, intensities, ion_masses, ion_ppm):
    ion_flags = {}

    def add_ion_flag(ion_idx, idx, diff, intensity):
        try:
            existing_ion_flag = ion_flags[ion_idx]
            if existing_ion_flag[2] < intensity:
                ion_flags[ion_idx] = (idx, diff, intensity)
        except KeyError:
            ion_flags[ion_idx] = (idx, diff, intensity)

    ion_masses_iter = iter(ion_masses)
    ion_mass_prev, ion_idx, ion_mass = 0, 0, next(ion_masses_iter)

    for idx, mass in enumerate(masses):
        try:
            while mass > ion_mass:
                ion_mass_prev, ion_idx, ion_mass = ion_mass, ion_idx+1, next(ion_masses_iter)
        except StopIteration:
            ion_mass = float('inf')

        error = mass * ion_ppm
        diff_prev, diff_next = mass-ion_mass_prev, ion_mass-mass
        if diff_prev < diff_next and diff_prev < error:
            add_ion_flag(ion_idx, idx, diff_prev, intensities[idx])
        if diff_prev > diff_next and diff_next < error:
            add_ion_flag(ion_idx+1, idx, diff_next, intensities[idx])

    return {value[0]:(key, value[1]) for key, value in ion_flags.items()}


def merge(b_ion_flags, y_ion_flags):
    ion_flags = {}

    for key, value in b_ion_flags.items():
        ion_flags['b' + str(value[0])] = key

    for key, value in y_ion_flags.items():
        try:
            b_value = b_ion_flags[key]
            if b_value[1] < value[1]:
                ion_flags['y' + str(value[0])] = key
        except KeyError:
            ion_flags['y' + str(value[0])] = key

    return ion_flags

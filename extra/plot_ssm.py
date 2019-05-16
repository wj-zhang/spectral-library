import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import ms2
from pyteomics import mgf
from pyteomics import mass
from fragments import ion_annotation
from spectrum_alignment import spectrum_alignment
import matplotlib


matplotlib.use('Agg')

colors = {'a': 'green', 'b': 'blue', 'y': 'red', 'u': 'black', None: 'darkgrey'}
zorders = {'a': 2, 'b': 3, 'y': 3, 'u': 1, None: 0}


def plot(query_spectrum, library_spectrum, peptide=None, score=0, aa_mass=mass.std_aa_mass, ion_ppm=0.00002, out_file=None):
    query_masses = query_spectrum["m/z array"]
    query_intensities = query_spectrum["intensity array"]
    query_ions = [None] * len(query_masses)
    if peptide:
        ion_flags = ion_annotation(query_masses, query_intensities, peptide, aa_mass=aa_mass, ion_ppm=ion_ppm)
        for key, value in ion_flags.items():
            query_ions[value] = key

    library_masses = library_spectrum["m/z array"]
    library_intensities = library_spectrum["intensity array"]
    library_ions = [None] * len(library_masses)
    if peptide:
        ion_flags = ion_annotation(library_masses, library_intensities, peptide, aa_mass=aa_mass, ion_ppm=ion_ppm)
        for key, value in ion_flags.items():
            library_ions[value] = key

    alignment = spectrum_alignment(query_spectrum, library_spectrum, ion_ppm=ion_ppm)
    for i, j in alignment.items():
        query_ions[i] = query_ions[i] or 'u'
        library_ions[j] = library_ions[j] or 'u'

    _plot(query_masses, query_intensities, query_ions, library_masses, library_intensities, library_ions, peptide,
          score, out_file)


def _plot(query_masses, query_intensities, query_ions, library_masses, library_intensities, library_ions, peptide,
          score, out_file):
    _, ax = plt.subplots(figsize=(20, 10))

    intensities = query_intensities / np.max(query_intensities)
    for i, (mass, intensity, ion) in enumerate(zip(query_masses, intensities, query_ions)):
        ion_type = ion[0] if ion else None
        color, zorder = colors[ion_type], zorders[ion_type]
        ax.plot([mass, mass], [0, intensity], color=color, zorder=zorder)
        if ion and ion is not 'u':
            ax.text(mass-5, intensity+0.05, '{}'.format(ion), color=color, rotation=270)

    intensities = library_intensities / np.max(library_intensities)
    for i, (mass, intensity, ion) in enumerate(zip(library_masses, intensities, library_ions)):
        ion_type = ion[0] if ion else None
        color, zorder = colors[ion_type], zorders[ion_type]
        ax.plot([mass, mass], [0, -intensity], color=color, zorder=zorder)
        if ion and ion is not 'u':
            ax.text(mass-5, -intensity-0.05, '{}'.format(ion), color=color, rotation=270)

    ax.axhline(0, color='black')

    max_mz = max(np.max(query_masses), np.max(library_masses)) + 10
    min_mz = min(np.min(query_masses), np.min(library_masses)) - 10

    ax.set_xticks(np.arange(0, max_mz, 200))
    ax.set_xlim(min_mz, max_mz)
    y_ticks = np.arange(-1, 1.05, 0.25)
    y_tick_labels = np.arange(-1, 1.05, 0.25)
    y_tick_labels[y_tick_labels < 0] = -y_tick_labels[y_tick_labels < 0]
    y_tick_labels = ['{:.0%}'.format(l) for l in y_tick_labels]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels)
    ax.set_ylim(-1.15, 1.05)

    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.grid(b=True, which='major', color='lightgrey', linestyle='--', linewidth=1.0)
    ax.grid(b=True, which='minor', color='lightgrey', linestyle='--', linewidth=0.5)

    ax.tick_params(axis='both', which='both', labelsize='small')

    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')

    if peptide:
        ax.text(0.5, 1.06, '{}, Score: {:.2f}'.format(peptide, score), horizontalalignment='center',
                verticalalignment='bottom', fontsize='x-large', fontweight='bold', transform=plt.gca().transAxes)

    if out_file:
        plt.savefig(out_file, dpi=300, bbox_inches='tight')


def plot_ssm(query_file, library_file, psm, aa_mass=mass.std_aa_mass, ion_ppm=0.00002, out_file=None):
    query_spectrum, library_spectrum, peptide = None, None, None

    with mgf.read(query_file) as reader:
        for spectrum in reader:
            if int(spectrum['params']['title'][6:]) + 1 == psm['index']:
                query_spectrum = spectrum
                break

    reader = ms2.read(library_file)
    for spectrum in reader:
        if spectrum["params"]["scan"][0] == str(psm['library']):
            library_spectrum = spectrum
            peptide = spectrum['params']['modified seq']
            break

    if query_spectrum is None or library_spectrum is None:
        print('No query spectrum or library spectrum found')

    score = psm.get('score', 0.0)

    plot(query_spectrum, library_spectrum, peptide=peptide, score=score, aa_mass=aa_mass, ion_ppm=ion_ppm, out_file=out_file)


def plot_pickle(query, hit, aa_mass=mass.std_aa_mass, ion_ppm=0.00002, out_file=None):
    query_spectrum, library_spectrum, peptide = query, hit, hit['params']['modified seq']

    score = 0.0

    plot(query_spectrum, library_spectrum, peptide=peptide, aa_mass=aa_mass, ion_ppm=ion_ppm, out_file=out_file)

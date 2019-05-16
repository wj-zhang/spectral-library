import numpy as np


def spectrum_snr(intensities):
    if len(intensities) < 6:
        return 0.01

    largest_six = np.sort(intensities[np.argpartition(intensities, -6)[-6:]])
    average_from_two_to_six = sum(largest_six[:-1]) / 5
    median = np.median(intensities)

    return average_from_two_to_six / median

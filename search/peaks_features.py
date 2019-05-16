import pandas as pd


def read(source):
    features = pd.read_csv(source)
    features = features[['Feature Id', 'm/z', 'z', 'RT', 'RT Begin', 'RT End', 'Quality']]
    features = features[features['Quality'] > 15]
    return features


import numpy as np
from scipy.spatial.distance import cdist
from astropy.coordinates import angular_separation
from difflib import SequenceMatcher
import pandas as pd

df = pd.read_csv("UCC_cat_20230507.csv")
x, y = df['GLON'], df['GLAT']
pmRA, pmDE, plx = df['pmRA'], df['pmDE'], df['plx']
coords = np.array([x, y]).T
# Find the distances to all clusters, for all clusters
dist = cdist(coords, coords)
# Change distance to itself from 0 to inf
msk = dist == 0.
dist[msk] = np.inf

dups_fnames = []
for i, cl in enumerate(dist):
    idx = np.argsort(cl)[:10]
    for j in idx:
        d = round(angular_separation(x[i], y[i], x[j], y[j]) * 60, 2)
        if d < 5:
            for fnamei in df['fnames'][i].split(';'):
                for fnamej in df['fnames'][j].split(';'):
                    sm_ratio = SequenceMatcher(None, fnamei, fnamej).ratio()
                    if sm_ratio > 0.5:
                        print(fnamei, fnamej, d, sm_ratio)

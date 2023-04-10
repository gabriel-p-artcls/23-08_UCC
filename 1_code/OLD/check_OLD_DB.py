
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from astropy.coordinates import angular_separation
import pandas as pd


df_old = pd.read_csv("../1_code/OLD_DBs.csv")

"""
Obtain the distance to the closest cluster for each cluster in the DB. This is
helpful in finding duplicates with different naming conventions
"""
x, y = df_old['GLON'], df_old['GLAT']
pmRA, pmDE, plx = df_old['pmRA'], df_old['pmDE'], df_old['plx']
coords = np.array([x, y]).T

# Find the distances to all clusters, for all clusters
dist = cdist(coords, coords)
msk = dist == 0.
dist[msk] = np.inf
dist_idx = np.argmin(dist, 0)

cl_dists, clusts = [], []
for i, j in enumerate(dist_idx):
    # d = round(dist[i][j] * 60, 3)
    d = round(angular_separation(x[i], y[i], x[j], y[j]) * 60, 3)
    pm_d = np.sqrt((pmRA[i]-pmRA[j])**2 + (pmDE[i]-pmDE[j])**2)
    plx_d = abs(plx[i] - plx[j])
    if d < 5 and pm_d < 2 and plx_d < .05:
        cl_dists.append(dist[i][j])
        clusts.append(
            [df_old['DB'][i], df_old['ID'][i], df_old['DB'][j],
             df_old['ID'][j], d, round(pm_d, 2), round(plx_d, 4)])
idx = np.argsort(cl_dists)
clusts = np.array(clusts)

for cl in clusts[idx]:
    print(cl)

"""
Count the number of clusters with no data for a given parameter
"""
# no_plx, no_pmRA, no_pmDE, no_plxpm, just_hao, msk = 0, 0, 0, 0, 0, []
# for i, cl in df_old.iterrows():
#     if np.isnan(cl['plx']):
#         no_plx += 1
#     if np.isnan(cl['pmRA']):
#         no_pmRA += 1
#     if np.isnan(cl['pmDE']):
#         no_pmDE += 1
#     if np.isnan(cl['plx']) and np.isnan(cl['pmRA']) and np.isnan(cl['pmDE']):
#         no_plxpm += 1
#         msk.append(True)
#     else:
#         msk.append(False)
#     if cl['DB'] == 5:
#         just_hao += 1
# print(no_plx, no_pmRA, no_pmDE, no_plxpm, just_hao)
# msk = np.array(msk)

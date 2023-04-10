
import numpy as np
from difflib import SequenceMatcher
from scipy.spatial.distance import cdist
import pandas as pd

"""
Script to clean the HAO21 database from clusters listed multiple times
and to bring the naming to a common convention
"""

dbs_folder = '../0_data/databases/old/'

data_CG20 = pd.read_csv(
    dbs_folder + "CG20.csv", sep=',', comment='#', index_col=False)
# names_CG20 = list(data_CG20['Name'])
names_CG20 = [_.lower().replace('_', '').replace(' ', '').replace('-', '')
              for _ in data_CG20['Name']]

data = pd.read_csv(
    dbs_folder + "HAO21_original.csv", sep=',', comment='#', index_col=False)
names = list(data['Cluster'])

names_f = []
for cl in names:
    # Remove trailing spaces
    cl = cl.strip()
    # # Cleaned and renamed manually using CG20 as template
    # if 'vdBergh-Hagen' in cl:
    #     cl = cl.replace('vdBergh-Hagen_', 'BH')
    # Fix FSR clusters
    if cl.startswith("FSR"):
        if cl[4] == "0":
            cl_id = int(cl[4:])
            cl = "FSR_" + str(cl_id)
    # Fix ESO clusters
    if cl.startswith("ESO"):
        cl_id = cl[4:].split('-')
        if len(cl_id) == 1:
            cl_id = cl[4:].split('_')
        cl_id = str(int(cl_id[0])) + '_' + str(int(cl_id[1]))
        cl = "ESO_" + str(cl_id)

    names_f.append(cl.lower().replace('_', '').replace(' ', '').replace('-', ''))

dups = []
for i, cl in enumerate(names_f):
    try:
        # Only store if duplicate is found
        i_dup = i + 1 + names_f[i + 1:].index(cl)
    except:
        continue

    xy_d = np.sqrt(
        (data['GLON'][i] - data['GLON'][i_dup])**2
        + (data['GLAT'][i] - data['GLAT'][i_dup])**2)
    dups.append([cl, i, i_dup, names[i], round(xy_d * 60, 2)])

dups = np.array(dups).T
dups_d = dups[-1].astype(float)
d_sort = np.argsort(dups_d)[::-1]
drop_r = []
for j in d_sort:
    try:
        k = names_CG20.index(dups[0][j])
    except:
        k = None

    if k is not None:
        # This duplicated cluster is present in CG20
        print("In CG20   ", dups[3][j], 2 + int(dups[1][j]),
              2 + int(dups[2][j]), k)
        # Drop it from dataframe
        if int(dups[1][j]) in drop_r:
            print(int(dups[1][j]), "already stored")
        if int(dups[2][j]) in drop_r:
            print(int(dups[2][j]), "already stored")
        drop_r += [int(dups[1][j]), int(dups[2][j])]

    if k is None:
        # This duplicated cluster is not present in CG20
        print("Not in CG20", dups[3][j], 2 + int(dups[1][j]),
              2 + int(dups[2][j]), dups[4][j])
        # Keep the one with the largest number of members
        N1 = data['N'][int(dups[1][j])]
        N2 = data['N'][int(dups[2][j])]
        if N1 >= N2:
            if int(dups[2][j]) in drop_r:
                print(int(dups[2][j]), "already stored")
            drop_r.append(int(dups[2][j]))
        else:
            if int(dups[1][j]) in drop_r:
                print(int(dups[1][j]), "already stored")
            drop_r.append(int(dups[1][j]))

data_f = data.drop(index=drop_r)

names_f = []
for cl in data_f['Cluster']:
    # Remove trailing spaces
    cl = cl.strip()
    # Fix FSR clusters
    if cl.startswith("FSR"):
        if cl[4] == "0":
            cl_id = int(cl[4:])
            cl = "FSR_" + str(cl_id)
    # Fix ESO clusters
    if cl.startswith("ESO"):
        cl_id = cl[4:].split('-')
        if len(cl_id) == 1:
            cl_id = cl[4:].split('_')
        cl_id = str(int(cl_id[0])) + '_' + str(int(cl_id[1]))
        cl = "ESO_" + str(cl_id)

    names_f.append(cl)

data_f['Cluster'] = names_f
data_f.to_csv("HAO21_clean.csv", index=False)
breakpoint()

"""
"""


def print_dists(data):
    """
    """
    x, y = data['GLON'], data['GLAT']
    coords = np.array([x, y]).T

    # Find the distances to all clusters, for all clusters
    dist = cdist(coords, coords)
    msk = dist == 0.
    dist[msk] = np.inf
    dist_idx = np.argmin(dist, 0)

    cl_dists, clusts = [], []
    for i, j in enumerate(dist_idx):
        cl_dists.append(dist[i][j])
        clusts.append(
            [data['Cluster'][i], data['Cluster'][j], round(dist[i][j] * 60, 3)])
    idx = np.argsort(cl_dists)
    clusts = np.array(clusts)
    for cl in clusts[idx]:
        # if 'vdBergh-Hagen' in cl[0] or 'BH' in cl[0]:
        sim_index = similar(cl[0], cl[1])
        if sim_index > 0.5 and float(cl[2]) < 10:
            print(cl, sim_index)


def similar(a, b):
    """
    """
    a, b = a.strip(), b.strip()
    a = a.lower().replace('_', '').replace(' ', '').replace('-', '')
    b = b.lower().replace('_', '').replace(' ', '').replace('-', '')
    sim_index = SequenceMatcher(None, a, b).ratio()
    return round(sim_index, 3)


print_dists(data)
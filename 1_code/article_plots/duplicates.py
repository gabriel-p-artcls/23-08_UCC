
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("/home/gabriel/Github/UCC/add_New_DB/UCC_cat_20230626_in.csv")

dups = []
for i, row in df.iterrows():
    if str(row['dups_fnames']) != 'nan':
        dups.append(row['fnames'].split(';')[0])

clusters_groups = (
    'ocsn', 'cwwdl', 'ubc', 'ufmg', 'ryu', 'hsc', 'cwnu', 'mwsc', 'oc', 'fof',
    'upk', 'hxhwl', 'lisc', 'hxwhb', 'xdocc')
Ndups = {}
N_dup_all = 0
for cl_ini in clusters_groups:
    Nini = 0
    for clust in dups:
        if clust.startswith(cl_ini):
            Nini += 1
    Ndups[cl_ini] = Nini
    N_dup_all += Nini

print(N_dup_all)
print(Ndups)
breakpoint()


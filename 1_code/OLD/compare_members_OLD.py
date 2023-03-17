
from os import listdir
import pandas as pd
import numpy as np

prob_min = 0.5
max_mag = 18

"""
Compare our outputs with CG2020
"""

path0 = "../0_data/cantat_gaudin_et_al_2020/CG_2020_members.csv"
# fastMP
path1 = "/home/gabriel/Documentos/fastMP_out/all/"

# CG2020
data_CG20 = pd.read_csv(path0)
names_CG = np.array([_.lower() for _ in data_CG20['Cluster']])

# Raw data
folders = ('S1', 'S2', 'M1', 'L1', 'L2', 'L3')

files = listdir(path1)
for file in files:

    # CG2020
    msk_cl = names_CG == file.split('.')[0]
    data_CG_ids = [str(_) for _ in data_CG20[msk_cl]['GaiaDR2']]
    N_CG = len(data_CG_ids)

    # fastMP
    data1 = pd.read_csv(path1 + file)
    for fold in folders:
        try:
            data1 = pd.read_csv("../0_data/" + fold + '/' + file)
        except:
            pass
    data_p1 = pd.read_csv(path1 + file)
    msk1 = (data_p1['probs_final'] > prob_min)
    data1 = data1[msk1]
    msk2 = data1['Gmag'] < max_mag
    data1 = data1[msk2]
    data1_ids = [_.replace('Gaia EDR3 ', '') for _ in data1['EDR3Name']]
    N_fastMP = len(data1)

    ids_match = list(set(data_CG_ids) & set(data1_ids))
    N_ids_inters = len(ids_match)

    # Percent recovered from CG20
    perc_diff = 100 * (1 - (N_CG - N_ids_inters) / N_CG)

    print("{}: N_CG: {}, N_fastMP: {} --> {:.0f}".format(file, N_CG, N_fastMP, perc_diff))

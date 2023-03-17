
import pandas as pd
import numpy as np
from os import listdir
import matplotlib.pyplot as plt


# CG2020
path0 = "../0_data/cantat_gaudin_et_al_2020/CG_2020_members.csv"
data_CG = pd.read_csv(path0)
names_CG = np.array([_.lower() for _ in data_CG['Cluster']])

# fastMP
path_p1 = "/media/gabriel/rest/Dropbox_nosync/Papers/2023/GDR3_members/1_code/members_fastMP/out_1/"

# Nmembs = {
# 'ngc_6475.csv': 1251,
# }

files = listdir(path_p1)
# np.random.shuffle(files)

prob_min1 = .5

for file in files:
    print(file)

    if "_full" in file:
        continue

    # if file[:-4] not in ('Pismis_19',):
    #     continue

    data = pd.read_csv("../0_data/CG_2017/" + file.lower())

    # CG2020
    msk0 = names_CG == file[:-4].lower()
    data_CGm = data_CG[msk0]

    data1 = pd.read_csv(path_p1 + file)
    msk1 = (data1['probs'] > prob_min1)
    data1 = data1[msk1]

    plt.figure(figsize=(15, 9))
    plt.subplot(231)

    plt.scatter(data1['GLON'], data1['GLAT'], alpha=.5, c=data1['probs'],
                vmin=prob_min1,
                label=f"fastMP_{path_p1.split('/')[-2]}: {msk1.sum()}")
    plt.colorbar()
    plt.scatter(data_CGm['GLON'], data_CGm['GLAT'], c='k', marker='x', lw=.7,
                label=f"CG20: {len(data_CGm)}")
    plt.legend()

    plt.subplot(232)
    plt.scatter(data1['pmRA'], data1['pmDE'], alpha=.5, c=data1['probs'],
                vmin=prob_min1)
    plt.colorbar()
    plt.scatter(
        data_CGm['pmRA*'], data_CGm['pmDE'], c='k', marker='x', lw=.7)

    plt.subplot(233)
    plt.scatter(data1['BP-RP'], data1['Gmag'], c=data1['probs'],
                alpha=.5, vmin=prob_min1)
    plt.colorbar()
    plt.scatter(
        data_CGm['BP-RP'], data_CGm['Gmag'], c='k', marker='x', lw=.7,
        alpha=.5, zorder=3)
    plt.gca().invert_yaxis()

    plt.subplot(234)
    plt.hist(data_CGm['Plx'], color='k', alpha=.5, density=True)
    plt.hist(data1['Plx'], color='b', alpha=.5, density=True)
    plt.axvline(np.median(data1['Plx']), c='b', ls=':')

    plt.subplot(235)
    plt.hist(data1['probs'], color='b', alpha=.5)

    plt.suptitle(file)
    plt.show()


import pandas as pd
import numpy as np
from os import listdir
import matplotlib.pyplot as plt


path0 = "../0_data/cantat_gaudin_et_al_2020/CG_2020_members.csv"

# # pyUPMASK results
# path2 = "../2_pipeline/pyUPMASK/78_large_fr_flags/"
# sep2 = " "
# New method results
path2 = "../2_pipeline/new_method/78_large_fr_flags/M65/"
# path2 = "../2_pipeline/new_method/303_no_flags/256_M51/"
sep2 = ","

files = listdir(path2)
np.random.shuffle(files)

prob_min = 0.5

for file in files:
    print(file)
    # if file != "LP_5.csv":
    #     continue

    # CG2020
    data_Pmemb = pd.read_csv(path0)
    msk0 = data_Pmemb['Cluster'] == file.split('.')[0]
    data0 = data_Pmemb[msk0]

    data2 = pd.read_csv(path2 + file, sep=sep2)
    lon_min, lon_max = data2['GLON'].min(), data2['GLON'].max()
    lat_min, lat_max = data2['GLAT'].min(), data2['GLAT'].max()
    msk2 = data2['probs_final'] > prob_min
    data2_18 = data2[msk2 & (data2['Gmag'] <= 18)]
    data2_inf = data2[msk2 & (data2['Gmag'] > 18)]

    plt.subplot(221)
    plt.scatter(data2_inf['GLON'], data2_inf['GLAT'], c='k', marker='^')
    plt.scatter(data2_18['GLON'], data2_18['GLAT'], c=data2_18['Gmag'],
                label=len(data2_18), cmap='Reds')
    plt.colorbar()
    plt.scatter(
        data0['GLON'], data0['GLAT'], c='b', marker='x', alpha=.5,
        label=len(data0))
    plt.legend()
    plt.xlim(lon_min, lon_max)
    plt.ylim(lat_min, lat_max)

    plt.subplot(222)
    plt.scatter(data2_inf['pmRA'], data2_inf['pmDE'], marker='^', c='k',
                alpha=.5)
    plt.scatter(data0['pmRA*'], data0['pmDE'], c='b', alpha=.5)
    plt.scatter(data2_18['pmRA'], data2_18['pmDE'], marker='x', c='r',
                alpha=.5)

    plt.subplot(223)
    plt.hist(data2_inf['Plx'], color='k', alpha=.5)
    plt.hist(data0['Plx'], color='b', alpha=.5)
    plt.hist(data2_18['Plx'], color='r', alpha=.5)

    plt.subplot(224)
    plt.scatter(data0['BP-RP'], data0['Gmag'], c='b', alpha=.5)
    plt.scatter(data2_18['BP-RP'], data2_18['Gmag'], marker='x', c='r',
                alpha=.5)
    plt.scatter(data2_inf['BP-RP'], data2_inf['Gmag'], marker='^', c='k',
                alpha=.5)
    plt.gca().invert_yaxis()

    plt.suptitle(file)
    plt.show()

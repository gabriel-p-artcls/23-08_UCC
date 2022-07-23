
import pandas as pd
import numpy as np
from os import listdir
import matplotlib.pyplot as plt


# # pyUPMASK results
# path2 = "../2_pipeline/pyUPMASK/78_large_fr_flags/"
# sep2 = " "

# Method 1
path1 = "../2_pipeline/new_method/78_large_fr_flags/M51/"

# Method 2
path2 = "../2_pipeline/new_method/78_large_fr_flags/M59/"
# path2 = "../2_pipeline/new_method/303_no_flags/256_M51/"

sep = ","
files = listdir(path2)

prob_min = 0.5

for file in files:
    print(file)
    # if file != "LP_5.csv":
    #     continue

    # Method 1
    data1 = pd.read_csv(path1 + file, sep=sep)
    msk1 = data1['probs_final'] > prob_min
    data1 = data1[msk1]

    data2 = pd.read_csv(path2 + file, sep=sep)
    msk2 = data2['probs_final'] > prob_min
    data2 = data2[msk2]

    plt.subplot(231)
    plt.scatter(data1['GLON'], data1['GLAT'], c='b', alpha=.5,
                label=len(data1))
    plt.scatter(data2['GLON'], data2['GLAT'], c='r', alpha=.5,
                label=len(data2))
    plt.legend()

    plt.subplot(232)
    plt.scatter(data1['pmRA'], data1['pmDE'], c='b', alpha=.5)
    plt.scatter(data2['pmRA'], data2['pmDE'], c='r', alpha=.5)

    plt.subplot(233)
    plt.hist(data1['Plx'], color='b', alpha=.5)
    plt.hist(data2['Plx'], color='r', alpha=.5)

    plt.subplot(234)
    plt.scatter(data1['BP-RP'], data1['Gmag'], c='b', alpha=.5)
    plt.scatter(data2['BP-RP'], data2['Gmag'], c='r', alpha=.5)
    plt.gca().invert_yaxis()

    plt.subplot(235)
    plt.hist(data1['probs_final'], color='b', alpha=.5)
    plt.hist(data2['probs_final'], color='r', alpha=.5)

    plt.suptitle(file)
    plt.show()

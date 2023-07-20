
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
Print/plot data on the CG20 clusters.
"""

# date = "0702"
date = "0712"
clpath = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/"
final_dbs_path = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/UCC_cat_20230702_out.csv"

# CG2020
path0 = "../0_data/CG_2020_members.csv.gz"

print("Reading fastMP output file...")
fastMP_db = pd.read_csv(final_dbs_path)
fnames = [_.split(';')[0] for _ in fastMP_db['fnames']]


def main():
    plot('NGC_3144')
    breakpoint()


def plot(clname):

    fname0 = clname.lower().replace('_', '').replace('-', '')
    i = fnames.index(fname0)
    row = fastMP_db.iloc[i]
    file_name = clpath + row['quad'] + '/datafiles/' + fname0 + '.parquet'
    cl2 = pd.read_parquet(file_name)

    probs = cl2['probs'].values
    msk = probs > .5
    Nmembs = 25 #int(row['N_membs'])
    if msk.sum() < Nmembs:
        # Select the 'N_membs' stars with the largest probabilities
        idx = np.argsort(probs)[::-1][:Nmembs]
        msk = np.full(len(probs), False)
        msk[idx] = True
    cl2 = cl2[msk]

    data_CG = pd.read_csv(path0)

    msk = data_CG['Cluster'] == clname
    cl1 = data_CG[msk]

    plt.suptitle(f"{clname}: N_C20={len(cl1)}, N_fMP={len(cl2)}")
    plt.subplot(221)
    plt.scatter(cl1['GLON'], cl1['GLAT'], alpha=.5, label='CG20')
    plt.scatter(cl2['GLON'], cl2['GLAT'], alpha=.5, label='fastMP')
    plt.legend()
    plt.subplot(222)
    plt.scatter(cl1['pmRA*'], cl1['pmDE'], alpha=.5)
    plt.scatter(cl2['pmRA'], cl2['pmDE'], alpha=.5)
    plt.subplot(223)
    plt.hist(cl1['Plx'], alpha=.5)
    plt.hist(cl2['Plx'], alpha=.5)
    plt.subplot(224)
    plt.scatter(cl1['BP-RP'], cl1['Gmag'], alpha=.5)
    plt.scatter(cl2['BP-RP'], cl2['Gmag'], alpha=.5)
    plt.gca().invert_yaxis()
    plt.show()


if __name__ == '__main__':
    # plt.style.use('science')
    main()

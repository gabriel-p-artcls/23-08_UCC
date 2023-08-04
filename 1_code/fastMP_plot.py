
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
Print/plot data on the CG20 clusters.
"""

clpath = "/home/gabriel/Github/UCC/"
final_dbs_path = clpath + "add_New_DB/UCC_cat_20230702.csv"

print("Reading fastMP output file...")
fastMP_db = pd.read_csv(final_dbs_path)
fnames = [_.split(';')[0] for _ in fastMP_db['fnames']]


def main():
    plot('ngc2516')
    breakpoint()


def plot(fname0):

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
    
    membs = cl2[msk]
    field = cl2[~msk]

    plt.suptitle(f"{fname0}: N_fMP={len(membs)}")
    plt.subplot(221)
    plt.scatter(field['GLON'], field['GLAT'], alpha=.35, c='grey', s=2)
    plt.scatter(membs['GLON'], membs['GLAT'], facecolors='none', ec='r', lw=2, label='fastMP')
    plt.subplot(222)
    plt.scatter(field['pmRA'], field['pmDE'], alpha=.35, c='grey', s=2)
    plt.scatter(membs['pmRA'], membs['pmDE'], facecolors='none', ec='r', lw=2)
    plt.subplot(223)
    plx_min, plx_max = min(membs['Plx']), max(membs['Plx'])
    msk = (field['Plx'] > plx_min) & (field['Plx'] < plx_max)
    plt.hist(field['Plx'][msk], alpha=.35, color='grey', density=True)
    plt.hist(membs['Plx'], alpha=.5, color='r', density=True)
    plt.xlim(plx_min, plx_max)
    plt.subplot(224)
    plt.scatter(field['BP-RP'], field['Gmag'], alpha=.35, c='grey', s=2)
    plt.scatter(membs['BP-RP'], membs['Gmag'], facecolors='none', ec='r', lw=2)
    plt.xlim(min(membs['BP-RP']) - .2, max(membs['BP-RP']) + .2)
    plt.ylim(max(membs['Gmag']) + .5, min(membs['Gmag']) - .5)
    # plt.gca().invert_yaxis()
    plt.show()


if __name__ == '__main__':
    # plt.style.use('science')
    main()

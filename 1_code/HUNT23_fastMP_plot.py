
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord


clpath = "/home/gabriel/Descargas/out/"
final_dbs_path = "/home/gabriel/Github/UCC/add_New_DB/UCC_cat_20230520.csv"
hunt23_membs_path = "../0_data/hunt23_members.csv.gz"
print("Reading HUNT23 members...")
hunt23_membs = pd.read_csv(hunt23_membs_path)
print("Reading fastMP output file...")
fastMP_db = pd.read_csv(final_dbs_path)
fnames = [_.split(';') for _ in fastMP_db['fnames']]


def main():
    hunt23_name = "Berkeley_29"
    make_plot(hunt23_name)
    breakpoint()


def make_plot(hunt23_name):

    fastMP_name = hunt23_name.lower().replace('_', '').replace(' ', '')

    if hunt23_name.startswith('VDBH_'):
        hunt23_name = 'BH_' + hunt23_name.split('_')[1]
    if hunt23_name.startswith('VDB_'):
        hunt23_name = 'vdBergh_' + hunt23_name.split('_')[1]

    # Read HUNT23 members data
    msk1 = hunt23_membs['name'] == hunt23_name
    hunt23_cl = hunt23_membs[msk1]

    # Read fastMP output catalogue
    i = np.nan
    for i, fname in enumerate(fnames):
        if fastMP_name in fname:
            break
    if np.isnan(i):
        print("Cluster not found")
        return

    row = fastMP_db.iloc[i]
    fname0 = row['fnames'].split(';')[0]
    # Search for the cluster csv file
    file_name = clpath + row['quad'] + '/datafiles/' + fname0 + '.csv.gz'
    try:
        cl_d = pd.read_csv(file_name)
    except:
        print(f"Not found: {file_name}")

    # probs = cl_d['probs'].values
    # msk = probs > 0.5
    # if msk.sum() < 25:
    #     # Select the 'N_membs_min' stars with the largest probabilities
    #     idx = np.argsort(probs)[::-1][:25]
    #     msk = np.full(len(probs), False)
    #     msk[idx] = True
    # fastMP_cl = cl_d[msk]

    probs = cl_d['probs'].values
    msk = probs > 0.5
    Nmembs = int(row['Nmembs'])
    if msk.sum() < Nmembs:
        # Select the 'N_membs' stars with the largest probabilities
        idx = np.argsort(probs)[::-1][:Nmembs]
        msk = np.full(len(probs), False)
        msk[idx] = True
    fastMP_cl = cl_d[msk]

    plot(hunt23_cl, fastMP_cl)


def plot(clust1, clust2):
    """
    """
    plt.suptitle(f"N_H23={len(clust1)}, N_fMP={len(clust2)}")
    plt.subplot(231)
    plt.scatter(clust1['l'], clust1['b'], alpha=.5, label='HUNT23')
    plt.scatter(clust2['GLON'], clust2['GLAT'], alpha=.5, label='fastMP')
    # plt.colorbar()
    plt.legend()
    plt.subplot(232)
    plt.scatter(clust1['pmra'], clust1['pmdec'], alpha=.5)
    plt.scatter(clust2['pmRA'], clust2['pmDE'], alpha=.5)
    plt.subplot(233)
    plt.hist(clust1['parallax'], alpha=.5, density=True)
    plt.hist(clust2['Plx'], alpha=.5, density=True)
    plt.subplot(234)
    plt.scatter(clust1['bp_rp'], clust1['phot_g_mean_mag'], alpha=.5)
    plt.scatter(clust2['BP-RP'], clust2['Gmag'], alpha=.5)
    plt.gca().invert_yaxis()
    plt.subplot(235)
    plt.hist(clust1['probability'], alpha=.5) #, density=True)
    plt.hist(clust2['probs'], alpha=.5) #, density=True)
    plt.xlabel('Probabilities')

    plt.subplot(236)
    d_pc = 1000 / clust1['parallax'].values
    d_pc = np.clip(d_pc, 20000, 1)
    cc = SkyCoord(
        ra=clust1['ra'].values*u.degree, dec=clust1['dec'].values*u.degree,
        distance=d_pc*u.pc)
    x, y, z = cc.cartesian.x, cc.cartesian.y, cc.cartesian.z
    plt.scatter(x, y, alpha=.5)

    d_pc = 1000 / clust2['Plx'].values
    d_pc = np.clip(d_pc, 20000, 1)
    cc = SkyCoord(
        ra=clust2['RA_ICRS'].values*u.degree,
        dec=clust2['DE_ICRS'].values*u.degree, distance=d_pc*u.pc)
    x, y, z = cc.cartesian.x, cc.cartesian.y, cc.cartesian.z
    plt.scatter(x, y, alpha=.5)

    plt.xlabel('X')
    plt.ylabel('Y')

    plt.show()


if __name__ == '__main__':
    # plt.style.use('science')
    main()

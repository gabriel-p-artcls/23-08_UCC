
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord


hunt23_membs_path = "../0_data/hunt23_members.parquet"
print("Reading HUNT23 members...")
hunt23_membs = pd.read_parquet(hunt23_membs_path)

# date = "0712"
# clpath = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/"
# final_dbs_path = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/UCC_cat_20230702_out.csv"
clpath = "/home/gabriel/Github/UCC/"
final_dbs_path = clpath + "add_New_DB/UCC_cat_20230702.csv"
print("Reading fastMP output file...")
fastMP_db = pd.read_csv(final_dbs_path)
fnames = [_.split(';') for _ in fastMP_db['fnames']]


def main():
    hunt23_name = "Haffner_9"
    make_plot(hunt23_name)
    breakpoint()


def make_plot(hunt23_name, Nmembs_manual=None):
    """
    """
    fastMP_name = hunt23_name.lower().replace('_', '').replace(
        ' ', '').replace('-', '').replace('.', '').replace('+', 'p')

    if hunt23_name.startswith('VDBH_'):
        hunt23_name = 'BH_' + hunt23_name.split('_')[1]
    if hunt23_name.startswith('VDB_'):
        hunt23_name = 'vdBergh_' + hunt23_name.split('_')[1]

    # Read HUNT23 members data
    msk1 = hunt23_membs['Name'] == hunt23_name
    hunt23_cl = hunt23_membs[msk1]
    msk = hunt23_cl['Prob'] > .5
    print(f"H23: {len(hunt23_cl) - msk.sum()} stars removed with P<0.5")
    hunt23_cl = hunt23_cl[msk]

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
    # Search for the cluster file
    file_name = clpath + row['quad'] + '/datafiles/' + fname0 + '.parquet'
    try:
        cl_d = pd.read_parquet(file_name)
    except:
        print(f"Not found: {file_name}")

    probs = cl_d['probs'].values
    if Nmembs_manual is None:
        msk = probs > .5
        Nmembs = 25 #int(row['N_membs'])
    else:
        Nmembs = Nmembs_manual
        msk = probs > 10
    if msk.sum() < Nmembs:
        # Select the 'N_membs' stars with the largest probabilities
        idx = np.argsort(probs)[::-1][:Nmembs]
        msk = np.full(len(probs), False)
        msk[idx] = True
    fastMP_cl = cl_d[msk]

    plot(hunt23_name, hunt23_cl, fastMP_cl)


def plot(name, clust1, clust2):
    """
    """
    plt.suptitle(f"{name}, N_H23={len(clust1)}, N_fMP={len(clust2)}")
    plt.subplot(231)
    plt.scatter(clust1['GLON'], clust1['GLAT'], alpha=.5, label='HUNT23')
    plt.scatter(clust2['GLON'], clust2['GLAT'], alpha=.5, label='fastMP')
    # plt.colorbar()
    plt.legend()

    plt.subplot(232)
    plt.scatter(clust1['pmRA'], clust1['pmDE'], alpha=.5)
    plt.scatter(clust2['pmRA'], clust2['pmDE'], alpha=.5)
    # pmra_t, pmdec_t = 4.74*clust1['pmRA']/clust1['Plx'],\
    #     4.74*clust1['pmDE']/clust1['Plx']
    # plt.scatter(pmra_t, pmdec_t, alpha=.5)
    # pmra_t, pmdec_t = 4.74*clust2['pmRA']/clust2['Plx'],\
    #     4.74*clust2['pmDE']/clust2['Plx']
    # plt.scatter(pmra_t, pmdec_t, alpha=.5)

    plt.subplot(233)
    plt.hist(clust1['Plx'], alpha=.5, density=True)
    plt.hist(clust2['Plx'], alpha=.5, density=True)
    plt.subplot(234)
    plt.scatter(clust1['BP-RP'], clust1['Gmag'], alpha=.5)
    plt.scatter(clust2['BP-RP'], clust2['Gmag'], alpha=.5)
    plt.gca().invert_yaxis()
    plt.subplot(235)
    plt.hist(clust1['Prob'], alpha=.5) #, density=True)
    plt.hist(clust2['probs'], alpha=.5) #, density=True)
    plt.xlabel('Probabilities')

    plt.subplot(236)
    d_pc = 1000 / clust1['Plx'].values
    d_pc = np.clip(d_pc, 1, 20000)
    # cc = SkyCoord(
    #     ra=clust1['RA_ICRS'].values*u.degree, dec=clust1['DE_ICRS'].values*u.degree,
    #     distance=d_pc*u.pc)
    # x, y, z = cc.cartesian.x, cc.cartesian.y, cc.cartesian.z
    # cc = SkyCoord(
    #     l=clust1['GLON'].values*u.degree, b=clust1['GLAT'].values*u.degree,
    #     distance=d_pc*u.pc, frame='galactic')
    # x, y, _ = cc.cartesian.x, cc.cartesian.y, cc.cartesian.z

    import astropy.coordinates as coord
    gc_frame = coord.Galactocentric()
    # Galactic coordinates.
    eq = SkyCoord(ra=clust1['RA_ICRS'].values*u.degree, dec=clust1['DE_ICRS'].values*u.degree, frame='icrs')
    lb = eq.transform_to('galactic')
    lon = lb.l.wrap_at(180 * u.deg).radian * u.radian
    lat = lb.b.radian * u.radian
    coords = SkyCoord(l=lon, b=lat, distance=d_pc*u.pc, frame='galactic')
    # Galactocentric coordinates.
    c_glct = coords.transform_to(gc_frame)
    x, y, z = c_glct.x, c_glct.y, c_glct.z
    # x_kpc, y_kpc, z_kpc = x_pc.value/1000, y_pc.value/1000, z_pc.value/1000
    plt.scatter(x, y, alpha=.5)
    xmin1, xmax1 = min(x.value), max(x.value)
    ymin1, ymax1 = min(y.value), max(y.value)

    d_pc = 1000 / clust2['Plx'].values
    d_pc = np.clip(d_pc, 1, 20000)
    # cc = SkyCoord(
    #     ra=clust2['RA_ICRS'].values*u.degree,
    #     dec=clust2['DE_ICRS'].values*u.degree, distance=d_pc*u.pc)
    # x, y, z = cc.cartesian.x, cc.cartesian.y, cc.cartesian.z
    # cc = SkyCoord(
    #     l=clust2['GLON'].values*u.degree, b=clust2['GLAT'].values*u.degree,
    #     distance=d_pc*u.pc, frame='galactic')
    # x, y, _ = cc.cartesian.x, cc.cartesian.y, cc.cartesian.z

    eq = SkyCoord(ra=clust2['RA_ICRS'].values*u.degree, dec=clust2['DE_ICRS'].values*u.degree, frame='icrs')
    lb = eq.transform_to('galactic')
    lon = lb.l.wrap_at(180 * u.deg).radian * u.radian
    lat = lb.b.radian * u.radian
    coords = SkyCoord(l=lon, b=lat, distance=d_pc*u.pc, frame='galactic')
    # Galactocentric coordinates.
    c_glct = coords.transform_to(gc_frame)
    x, y, z = c_glct.x, c_glct.y, c_glct.z
    # x_kpc, y_kpc, z_kpc = x_pc.value/1000, y_pc.value/1000, z_pc.value/1000
    plt.scatter(x, y, alpha=.5)
    xmin2, xmax2 = min(x.value), max(x.value)
    ymin2, ymax2 = min(y.value), max(y.value)

    xmin, xmax = min(xmin1, xmin2), max(xmax1, xmax2)
    ymin, ymax = min(ymin1, ymin2), max(ymax1, ymax2)
    xr = xmax - xmin
    yr = ymax - ymin
    xyr = max(xr, yr) * .5
    xc = .5 * (xmin + xmax)
    yc = .5 * (ymin + ymax)

    plt.xlim(xc - xyr, xc + xyr)
    plt.ylim(yc - xyr, yc + xyr)

    plt.xlabel('X')
    plt.ylabel('Y')

    plt.show()


if __name__ == '__main__':
    # plt.style.use('science')
    main()

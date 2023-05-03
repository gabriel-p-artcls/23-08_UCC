
import os
from astropy.coordinates import angular_separation
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import spatial
import main_process_GDR3_query as G3Q
import time as t

# LOCAL RUN
# insert at 1, 0 is the script path (or '' in REPL)
import sys
sys.path.insert(1, '/home/gabriel/Github/fastmp/')  # Path to fastMP
GAIADR3_path = '/media/gabriel/backup/gabriel/GaiaDR3/'
frames_path = GAIADR3_path + 'datafiles_G20/'
frames_ranges = GAIADR3_path + 'files_G20/frame_ranges.txt'
comb_DBs = "/home/gabriel/Github/web_sites/UCC/datafiles/UCC_cat_20230503.csv"

# # CLUSTER RUN
# # Path to the database with Gaia DR3 data
# GAIADR3_path = '/home/gperren/GaiaDR3/'
# frames_path = GAIADR3_path + 'datafiles_G20/'
# # File that contains the regions delimited by each frame (in this folder)
# frames_ranges = GAIADR3_path + 'files_G20/frame_ranges.txt'
# # Full database of clusters (in this folder)
# comb_DBs = "UCC_cat_20230503.csv"

from fastmp import fastMP

# Maximum magnitude to retrieve
max_mag = 20

# Run ID
run_ID = os.path.basename(__file__).split('.')[0][-1:]
# Create output folder in not present
Path(f"out_{run_ID}").mkdir(parents=True, exist_ok=True)
# Output folder
out_path = f"out_{run_ID}/"


def main(N_cl_extra=10):
    """
    N_cl_extra: number of extra clusters in frame to detect
    """
    # Read data
    frames_data, database_full = read_input()

    # Parameters used to search for close-by clusters
    xys = np.array([
        database_full['GLON'].values, database_full['GLAT'].values]).T
    tree = spatial.cKDTree(xys)
    close_cl_idx = tree.query(xys, k=N_cl_extra + 1)

    # Split into 5 jobs
    N_r = int(len(database_full) / 5.)
    if run_ID == '1':
        clusters_list = database_full[:N_r]
    elif run_ID == '2':
        clusters_list = database_full[N_r:2*N_r]
    elif run_ID == '3':
        clusters_list = database_full[2*N_r:3*N_r]
    elif run_ID == '4':
        clusters_list = database_full[3*N_r:4*N_r]
    elif run_ID == '5':
        clusters_list = database_full[4*N_r:]

    # # Full list
    # clusters_list = database_full
    # # Shuffle
    # clusters_list = clusters_list.sample(frac=1).reset_index(drop=True)

    for index, cl in clusters_list.iterrows():

        # if 'teutschj084364511' not in cl['fnames']:
        #     continue

        fname0 = cl['fnames'].split(';')[0]

        print("\n" + str(index), fname0, cl['GLON'], cl['GLAT'],
              cl['pmRA'], cl['pmDE'], cl['plx'])

        # Get close clusters coords
        centers_ex = get_close_cls(
            index, database_full, close_cl_idx, cl['dups_fnames'])

        # Generate frame
        data = get_frame(frames_path, frames_data, cl)
        # Store full file
        # data.to_csv(out_path + fname0 + "_full.csv", index=False)
        # Read from file
        # data = pd.read_csv(out_path + fname0 + "_full.csv")

        # Extract center coordinates
        xy_c, vpd_c, plx_c = (cl['GLON'], cl['GLAT']), None, None
        if not np.isnan(cl['pmRA']):
            vpd_c = (cl['pmRA'], cl['pmDE'])
        if not np.isnan(cl['plx']):
            plx_c = cl['plx']

        fixed_centers = False
        if vpd_c is None and plx_c is None:
            fixed_centers = True

        # Generate input data array for fastMP
        X = np.array([
            data['GLON'].values, data['GLAT'].values, data['pmRA'].values,
            data['pmDE'].values, data['Plx'].values, data['e_pmRA'].values,
            data['e_pmDE'].values, data['e_Plx'].values])

        # Process with fastMP
        start = t.time()
        while True:
            print("Fixed centers?:", fixed_centers)
            probs_all = fastMP(
                xy_c=xy_c, vpd_c=vpd_c, plx_c=plx_c, centers_ex=centers_ex,
                fixed_centers=fixed_centers).fit(X)

            bad_center = check_centers(
                *X[:5, :], xy_c, vpd_c, plx_c, probs_all)

            if bad_center == '000' or fixed_centers is True:
                break
            else:
                # print("Re-run with fixed_centers = True")
                fixed_centers = True

        bad_center = check_centers(*X[:5, :], xy_c, vpd_c, plx_c, probs_all)
        print("*** {}: Nmembs(P>0.5)={}, t={}, cents={}".format(
              fname0, (probs_all > 0.5).sum(),
              round(t.time() - start, 1), bad_center))

        # Write member stars for cluster
        df = filter_stars_not_used(data, probs_all)
        df.to_csv(out_path + fname0 + ".csv.gz", index=False,
                  compression='gzip')


def read_input():
    """
    Read input file with the list of clusters to process
    """
    frames_data = pd.read_csv(frames_ranges)
    database_full = pd.read_csv(comb_DBs)
    return frames_data, database_full


def get_frame(frames_path, frames_data, cl):
    """
    """
    if not np.isnan(cl['plx']):
        c_plx = cl['plx']
    else:
        c_plx = None

    if c_plx is None:
        box_s_eq = 1
    else:
        if c_plx > 15:
            box_s_eq = 30
        elif c_plx > 4:
            box_s_eq = 20
        elif c_plx > 2:
            box_s_eq = 6
        elif c_plx > 1.5:
            box_s_eq = 5
        elif c_plx > 1:
            box_s_eq = 4
        elif c_plx > .75:
            box_s_eq = 3
        elif c_plx > .5:
            box_s_eq = 2
        elif c_plx > .25:
            box_s_eq = 1.5
        elif c_plx > .1:
            box_s_eq = 1
        else:
            box_s_eq = .5

    if 'Ryu' in cl['ID']:
        box_s_eq = 10 / 60

    # Filter by parallax if possible
    plx_min = -2
    if c_plx is not None:
        if c_plx > 15:
            plx_p = 5
        elif c_plx > 4:
            plx_p = 2
        elif c_plx > 2:
            plx_p = 1
        elif c_plx > 1:
            plx_p = .7
        else:
            plx_p = .6
        plx_min = c_plx - plx_p

    data = G3Q.run(frames_path, frames_data, cl['RA_ICRS'], cl['DE_ICRS'],
                   box_s_eq, plx_min, max_mag)

    return data


def get_close_cls(idx, database_full, close_cl_idx, dups):
    """
    Get data on the closest clusters to the one being processed

    idx: Index to the cluster in the full list
    """
    # Indexes to the closest clusters in XY
    ex_cls_idx = close_cl_idx[1][idx][1:]

    duplicate_cls = []
    if str(dups) != 'nan':
        duplicate_cls = dups.split(';')

    centers_ex = []
    for i in ex_cls_idx:

        # Check if this close cluster is identified as a probable duplicate
        # of this cluster. If it is, do not add it to the list of extra
        # clusters in the frame
        skip_cl = False
        if duplicate_cls:
            for dup_fname_i in database_full['fnames'][i].split(';'):
                if dup_fname_i in duplicate_cls:
                    # print("skip", database_full['fnames'][i])
                    skip_cl = True
                    break
            if skip_cl:
                continue

        ex_cl_dict = {
            'xy': [database_full['GLON'][i], database_full['GLAT'][i]]}

        # Only use clusters with defined centers in PMs and Plx, otherwise
        # non-bonafide clusters disrupt the process for established clusters
        # like NGC 2516
        if np.isnan(database_full['pmRA'][i])\
                or np.isnan(database_full['plx'][i]):
            continue

        if not np.isnan(database_full['pmRA'][i]):
            ex_cl_dict['pms'] = [
                database_full['pmRA'][i], database_full['pmDE'][i]]
        if not np.isnan(database_full['plx'][i]):
            ex_cl_dict['plx'] = [database_full['plx'][i]]

        # print(database_full['ID'][i], ex_cl_dict)

        centers_ex.append(ex_cl_dict)

    return centers_ex


def check_centers(
    lon, lat, pmRA, pmDE, plx, xy_c, vpd_c, plx_c, probs_all, N_membs_min=25
):
    """
    """
    # Select high-quality members
    msk = probs_all > 0.5
    if msk.sum() < N_membs_min:
        idx = np.argsort(probs_all)[::-1][:N_membs_min]
        msk = np.full(len(probs_all), False)
        msk[idx] = True

    # Centers of selected members
    xy_c_f = np.nanmedian([lon[msk], lat[msk]], 1)
    vpd_c_f = np.nanmedian([pmRA[msk], pmDE[msk]], 1)
    plx_c_f = np.nanmedian(plx[msk])

    bad_center_xy, bad_center_pm, bad_center_plx = '0', '0', '0'

    # 5 arcmin maximum
    d_arcmin = angular_separation(xy_c_f[0], xy_c_f[1], xy_c[0], xy_c[1]) * 60
    if d_arcmin > 5:
        # print("d_arcmin: {:.1f}".format(d_arcmin))
        # print(xy_c, xy_c_f)
        bad_center_xy = '1'

    # Relative difference
    if vpd_c is not None:
        pm_max = []
        for vpd_c_i in abs(np.array(vpd_c)):
            if vpd_c_i > 10:
                pm_max.append(20)
            elif vpd_c_i > 1:
                pm_max.append(25)
            elif vpd_c_i > 0.1:
                pm_max.append(35)
            elif vpd_c_i > 0.01:
                pm_max.append(50)
            else:
                pm_max.append(70)
        pmra_p = 100 * abs(vpd_c_f[0] - vpd_c[0]) / (vpd_c[0] + 0.001)
        pmde_p = 100 * abs(vpd_c_f[1] - vpd_c[1]) / (vpd_c[1] + 0.001)
        if pmra_p > pm_max[0] or pmde_p > pm_max[1]:
            # print("pm: {:.2f} {:.2f}".format(pmra_p, pmde_p))
            # print(vpd_c, vpd_c_f)
            bad_center_pm = '1'

    # Relative difference
    if plx_c is not None:
        if plx_c > 0.2:
            plx_max = 25
        elif plx_c > 0.1:
            plx_max = 30
        elif plx_c > 0.05:
            plx_max = 35
        elif plx_c > 0.01:
            plx_max = 50
        else:
            plx_max = 70
        plx_p = 100 * abs(plx_c_f - plx_c) / (plx_c + 0.001)
        if abs(plx_p) > plx_max:
            # print("plx: {:.2f}".format(plx_p))
            # print(plx_c, plx_c_f)
            bad_center_plx = '1'

    bad_center = bad_center_xy + bad_center_pm + bad_center_plx

    return bad_center


def filter_stars_not_used(data, probs_all, prob_min=0):
    """
    """
    data['probs'] = np.round(probs_all, 5)
    msk = probs_all >= prob_min

    df = data[msk]
    # Order by probabilities
    df = df.sort_values('probs', ascending=False)

    return df


if __name__ == '__main__':
    main()

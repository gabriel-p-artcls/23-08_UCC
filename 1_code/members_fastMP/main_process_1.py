
import os
from pathlib import Path
import numpy as np
import pandas as pd
import csv
from scipy import spatial
import main_process_GDR3_query as G3Q
import time as t

# # LOCAL RUN
# # insert at 1, 0 is the script path (or '' in REPL)
# import sys
# sys.path.insert(1, '/home/gabriel/Github/fastmp/')  # Path to fastMP
# frames_path = '/media/gabriel/backup/gabriel/GaiaDR3/datafiles_G20/'
# comb_DBs = "../../2_pipeline/3_final_DB.csv"

# CLUSTER RUN
# Path to the database with Gaia DR3 data
frames_path = '/home/gperren/GaiaDR3/datafiles_G20/'
# Full database of clusters (in this folder)
comb_DBs = "3_final_DB.csv"

from fastmp import fastMP

# File that contains the regions delimited by each frame (in this folder)
frames_ranges = 'frame_ranges.txt'

# Maximum magnitude to retrieve
max_mag = 20

# Minimum probability to store
prob_min = 0.5

# Run ID
run_ID = os.path.basename(__file__).split('.')[0][-1:]
# Create output folder in not present
Path(f"out_{run_ID}").mkdir(parents=True, exist_ok=True)
# Output folder
out_path = f"out_{run_ID}/"


def main(N_cl_extra=5):
    """
    N_cl_extra: number of extra clusters in frame to detect
    """
    # Generate final table
    with open(out_path + f"table_{run_ID}.csv", "w") as f:
        f.write(
            "DB,DB_i,ID,fname,dups_fname,dups_names,UCC_ID,"
            + "Class,Class_v,GLON,GLAT,RA_ICRS,DE_ICRS,plx,pmRA,pmDE,RV,"
            + "N_membs\n")

    # Read data
    frames_data, database_full = read_input()

    # Parameters used to search for close-by clusters
    database_fnames = database_full['fname']
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

        # if cl['fname'] not in ('ngc2451a',):
        #     continue
        print(index, cl['fname'], cl['GLON'], cl['GLAT'], cl['pmRA'],
              cl['pmDE'], cl['plx'])

        # Get close clusters coords
        centers_ex = get_close_cls(
            database_full, database_fnames, close_cl_idx, cl['fname'],
            cl['dups_fnames'])

        # Generate frame
        data = get_frame(frames_path, frames_data, cl)
        # # Store full file
        # data.to_csv(out_path + cl['fname'] + "_full.csv", index=False)
        # # Read from file
        # data = pd.read_csv(out_path + cl['fname'] + "_full.csv")

        # Extract center coordinates
        xy_c, vpd_c, plx_c, fixed_centers = (cl['GLON'], cl['GLAT']), None,\
            None, False
        if not np.isnan(cl['pmRA']):
            vpd_c = (cl['pmRA'], cl['pmDE'])
        if not np.isnan(cl['plx']):
            plx_c = cl['plx']
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
            probs_all, N_membs = fastMP(
                xy_c=xy_c, vpd_c=vpd_c, plx_c=plx_c, centers_ex=centers_ex,
                fixed_centers=fixed_centers).fit(X)

            break_flaf = check_centers(X, xy_c, vpd_c, plx_c, probs_all)

            if break_flaf is True or fixed_centers is True:
                break
            else:
                # print("Re-run with fixed_centers = True")
                fixed_centers = True

        try:
            classif = get_classif(X[0], X[1], X[2], X[3], X[4], probs_all)
        except:
            classif = 'FFF'

        print("{}: Nmembs {}, P>0.5 {}, {}; t={}, {}".format(
              cl['fname'], N_membs, (probs_all > 0.5).sum(), classif,
              round(t.time() - start, 1), fixed_centers))

        # Write output
        write_out(run_ID, out_path, cl, data, probs_all, N_membs, classif)
        # write_out(run_ID, out_path, cl)


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


def get_close_cls(database_full, database_fnames, close_cl_idx, fname, dups):
    """
    Get data on the closest clusters to the one being processed
    """
    # Index to the cluster in the full list
    idx = database_full.index[database_fnames == fname][0]
    # Indexes to the closest clusters in XY
    ex_cls_idx = close_cl_idx[1][idx][1:]

    duplicate_cls = []
    if str(dups) != 'nan':
        duplicate_cls = dups.split(',')

    centers_ex = []
    for i in ex_cls_idx:

        # Check if this close cluster is identified as a probable duplicate
        # of this cluster. If it is, do not add it to the list of extra
        # clusters in the frame
        if duplicate_cls:
            if database_fnames[i] in duplicate_cls:
                continue

        ex_cl_dict = {
            'xy': [database_full['GLON'][i], database_full['GLAT'][i]]}
        if not np.isnan(database_full['pmRA'][i]):
            ex_cl_dict['pms'] = [
                database_full['pmRA'][i], database_full['pmDE'][i]]
        if not np.isnan(database_full['plx'][i]):
            ex_cl_dict['plx'] = [database_full['plx'][i]]

        if np.isnan(database_full['pmRA'][i]) or\
                np.isnan(database_full['plx'][i]):
            continue
        ex_cl_dict = {
            'xy': [database_full['GLON'][i], database_full['GLAT'][i]],
            'pms': [database_full['pmRA'][i], database_full['pmDE'][i]],
            'plx': [database_full['plx'][i]]
        }
        # print(database_full['ID'][i], ex_cl_dict)
        centers_ex.append(ex_cl_dict)

    return centers_ex


def get_classif(
    lon, lat, pmRA, pmDE, plx, probs_final, rad_max=2, N_membs_min=25
):
    """
    """
    # Filter stars beyond 2 times the 95th percentile distance from the center
    msk_membs = probs_final > 0.5
    if msk_membs.sum() < N_membs_min:
        msk_membs = probs_final > 0.25
    if msk_membs.sum() < N_membs_min:
        return 'DDD'

    xy_c = np.nanmedian([lon[msk_membs], lat[msk_membs]], 1)
    xy = np.array([lon[msk_membs], lat[msk_membs]]).T
    xy_rads = spatial.distance.cdist(xy, np.array([xy_c])).T[0]
    xy_rad = np.percentile(xy_rads, 95)
    xy = np.array([lon, lat]).T
    xy_rads = spatial.distance.cdist(xy, np.array([xy_c])).T[0]
    msk_xy = xy_rads < rad_max * xy_rad
    lon, lat, pmRA, pmDE, plx, probs_final = lon[msk_xy], lat[msk_xy],\
        pmRA[msk_xy], pmDE[msk_xy], plx[msk_xy], probs_final[msk_xy]

    # Select most probable members
    msk_membs = probs_final > 0.5
    if msk_membs.sum() < N_membs_min:
        msk_membs = probs_final > 0.25
    if msk_membs.sum() < N_membs_min:
        return 'DDD'

    lon_m, lat_m, pmRA_m, pmDE_m, plx_m = lon[msk_membs], lat[msk_membs],\
        pmRA[msk_membs], pmDE[msk_membs], plx[msk_membs]

    xy_c = np.nanmedian([lon_m, lat_m], 1)
    vpd_c = np.nanmedian([pmRA_m, pmDE_m], 1)
    plx_c = np.nanmedian(plx_m)

    # Median distances to centers for members
    xy = np.array([lon_m, lat_m]).T
    xy_dists = spatial.distance.cdist(xy, np.array([xy_c])).T[0]
    xy_50 = np.nanmedian(xy_dists)
    pm = np.array([pmRA_m, pmDE_m]).T
    pm_dists = spatial.distance.cdist(pm, np.array([vpd_c])).T[0]
    pm_50 = np.nanmedian(pm_dists)
    plx_dists = abs(plx_m - plx_c)
    plx_50 = np.nanmedian(plx_dists)
    # Count member stars within median distances
    N_memb_xy = (xy_dists < xy_50).sum()
    N_memb_pm = (pm_dists < pm_50).sum()
    N_memb_plx = (plx_dists < plx_50).sum()

    # Median distances to centers for field stars
    xy = np.array([lon[~msk_membs], lat[~msk_membs]]).T
    xy_dists_f = spatial.distance.cdist(xy, np.array([xy_c])).T[0]
    pm = np.array([pmRA[~msk_membs], pmDE[~msk_membs]]).T
    pm_dists_f = spatial.distance.cdist(pm, np.array([vpd_c])).T[0]
    plx_dists_f = abs(plx[~msk_membs] - plx_c)
    # Count field stars within median distances
    N_field_xy = (xy_dists_f < xy_50).sum()
    N_field_pm = (pm_dists_f < pm_50).sum()
    N_field_plx = (plx_dists_f < plx_50).sum()

    def ABCD_classif(Nm, Nf):
        """Obtain 'ABCD' classification"""
        if Nm == 0:
            return "D"
        if Nf == 0:
            return "A"
        N_ratio = Nm / Nf

        if N_ratio >= 1:
            cl = "A"
        elif N_ratio < 1 and N_ratio >= 0.5:
            cl = "B"
        elif N_ratio < 0.5 and N_ratio > 0.1:
            cl = "C"
        else:
            cl = "D"
        return cl

    c_xy = ABCD_classif(N_memb_xy, N_field_xy)
    c_pm = ABCD_classif(N_memb_pm, N_field_pm)
    c_plx = ABCD_classif(N_memb_plx, N_field_plx)
    classif = c_xy + c_pm + c_plx

    return classif


def check_centers(X, xy_c, vpd_c, plx_c, probs_all):
    """
    """
    msk = probs_all > 0.5
    lon, lat, pmRA, pmDE, plx = X[:5, msk]

    xy_c_f = np.median([lon, lat], 1)
    vpd_c_f = np.median([pmRA, pmDE], 1)
    plx_c_f = np.median(plx)

    break_flag = True

    # 5 arcmin maximum
    x_d, y_d = abs(xy_c_f - xy_c) / 60
    if x_d > 5 or y_d > 5:
        break_flag = False
        # print("xy: {:.1f} {:.1f}".format(x_d, y_d))
        # print(xy_c, xy_c_f)
        return break_flag

    # 5% relative difference maximum
    pm_max = 5
    if vpd_c is not None:
        pmra_p, pmde_p = 100 * (vpd_c_f - vpd_c) / vpd_c
        if abs(pmra_p) > pm_max or abs(pmde_p) > pm_max:
            # print("pm: {:.2f} {:.2f}".format(pmra_p, pmde_p))
            # print(vpd_c, vpd_c_f)
            break_flag = False
            return break_flag

    if plx_c is not None:
        if plx_c >= 4:
            plx_max = 0.5
        elif plx_c < 4 and plx_c > 2:
            plx_max = 0.2
        elif plx_c > 1 and plx_c < 2:
            plx_max = 0.1
        elif plx_c > 0.5 and plx_c < 1:
            plx_max = 0.05
        elif plx_c > 0.1 and plx_c < 0.5:
            plx_max = 0.01
        else:
            plx_max = 0.005
        plx_p = abs(plx_c_f - plx_c)
        if abs(plx_p) > plx_max:
            # print("plx: {:.2f}".format(plx_p))
            # print(plx_c, plx_c_f)
            break_flag = False
            return break_flag

    return break_flag

# def write_out(run_ID, out_path, cl):
def write_out(run_ID, out_path, cl, data, probs_all, N_membs, classif,
              prob_min=0.5):
    """
    """
    data['probs'] = np.round(probs_all, 2)
    msk = probs_all > prob_min

    if N_membs > msk.sum():
        idx = np.argsort(probs_all)[::-1][:N_membs]
        msk = np.full(len(probs_all), False)
        msk[idx] = True

    df = data[msk]
    N_membs = len(df)

    # Write member stars for cluster
    df.to_csv(out_path + cl['fname'] + ".csv.gz", index=False,
              compression='gzip')

    lon, lat = np.nanmedian(df['GLON']), np.nanmedian(df['GLAT'])
    ra, dec = np.nanmedian(df['RA_ICRS']), np.nanmedian(df['DE_ICRS'])
    plx = np.nanmedian(df['Plx'])
    pmRA, pmDE = np.nanmedian(df['pmRA']), np.nanmedian(df['pmDE'])
    RV = np.nanmedian(df['RV'])

    lon, lat = round(lon, 4), round(lat, 4)
    ra, dec = round(ra, 4), round(dec, 4)
    plx = round(plx, 4)
    pmRA, pmDE = round(pmRA, 4), round(pmDE, 4)
    RV = round(RV, 4)

    abcd_v = UCC_value(classif)

    df_row = pd.DataFrame(data={
        "DB": [cl['DB']], "DB_i": [cl['DB_i']], "ID": [cl['ID']],
        'fname': [cl['fname']], "dups_fnames": [cl['dups_fnames']],
        "UCC_ID": [cl['UCC_ID']], 'Class': classif, "Class_v": [abcd_v], 
        "GLON": [lon], "GLAT": [lat], "RA_ICRS": [ra], "DE_ICRS": [dec],
        "plx": [plx], "pmRA": [pmRA], "pmDE": [pmDE], "RV": [RV],
        "N_membs": [N_membs]})
    df_row.to_csv(out_path + f"table_{run_ID}.csv", mode='a', header=False,
                  index=False, na_rep='nan', quoting=csv.QUOTE_NONNUMERIC)


def UCC_value(abcd_c):
    """
    """
    class_2_num = {'A': 1., 'B': 0.5, 'C': 0.25, 'D': 0.1}
    c1, c2, c3 = abcd_c
    abcd_v = (class_2_num[c1] + class_2_num[c2] + class_2_num[c3]) / 3
    return round(abcd_v, 2)


if __name__ == '__main__':
    main()

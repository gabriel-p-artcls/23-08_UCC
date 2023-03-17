
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import spatial
# import time as t
import sys
import main_process_GDR3_query as G3Q
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/home/gabriel/Github/fastmp/')  # Path to fastMP
from fastmp import fastMP

# Path to the database with Gaia DR3 data
frames_path = '/media/gabriel/backup/gabriel/GaiaDR3/datafiles/'

# Full database of clusters (in this folder)
gaiadr3_path_full = "database_full.csv"
# File that contains the regions delimited by each frame (in this folder)
frames_ranges = 'frame_ranges.txt'

# Maximum magnitude to retrieve
max_mag = 19

# Minimum probability to store
prob_min = 0.5

# Run ID
run_ID = os.path.basename(__file__).split('.')[0][-1:]
# Create output folder
Path(f"out_{run_ID}").mkdir(parents=True, exist_ok=True)


def main(N_cl_extra=5):
    """
    N_cl_extra: number of extra clusters in frame to detect
    """

    # Path to the database to process
    gaiadr3_path_partial = f"process_partial_{run_ID}.csv"
    # Output folder
    out_path = f"out_{run_ID}/"

    # Read data
    clusters_list, frames_data, database_full = read_input(
        gaiadr3_path_partial, frames_ranges)

    # Parameters used to search for close-by clusters
    database_full_names = np.array([_.lower() for _ in database_full['ID']])
    xys = np.array([
        database_full['GLON'].values, database_full['GLAT'].values]).T
    tree = spatial.cKDTree(xys)
    close_cl_idx = tree.query(xys, k=N_cl_extra + 1)

    # aa = os.listdir("/media/gabriel/rest/Dropbox_nosync/Papers/2023/GDR3_members/1_code/members_fastMP/out_1/")
    # aa = [_[:-4] for _ in aa]
    # clusters_list = clusters_list.sample(frac=0.1)

    for index, cl in clusters_list.iterrows():
        # t0 = t.time()
        name = cl['ID']

        # if name in aa:
        # if name not in ('HE1564',):
        # # if 'LISC36' not in name:
            # continue
        print(index, name)

        # Get close clusters coords
        centers_ex = get_close_cls(
            database_full, database_full_names, close_cl_idx, name)

        # Generate frame
        # t1 = t.time()
        data = get_frame(frames_path, frames_data, cl)
        data.to_csv(out_path + name + "_full.csv", index=False)
        # data = pd.read_csv(out_path + name + "_full.csv")
        # t1f = t.time() - t1

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

        classif = get_classif(X[0], X[1], X[2], X[3], X[4], probs_all)

        print("{}: Nmembs {}, P>0.5 {}, {}".format(
              name, N_membs, (probs_all > 0.5).sum(), classif))

        # Write output
        write_out(out_path, name, data, probs_all, N_membs)

        # print("Times: {:.1f}+{:.1f} = {:.1f}".format(t1f, t2f, t.time() - t0))

        # Remove full file
        os.remove(out_path + name + "_full.csv")


def read_input(gaiadr3_path_partial, frame_ranges):
    """
    Read input file with the list of clusters to process
    """
    clusters_list = pd.read_csv(gaiadr3_path_partial)
    # Shuffle
    clusters_list = clusters_list.sample(frac=1).reset_index(drop=True)
    frames_data = pd.read_csv(frame_ranges)
    database_full = pd.read_csv(gaiadr3_path_full)
    return clusters_list, frames_data, database_full


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


def get_close_cls(GDR3_f, GDR3_f_names, close_cl_idx, name):
    """
    Get data on the closest clusters to the one being processed
    """
    # Index to the cluster in the full list
    idx = GDR3_f.index[GDR3_f_names == name.lower()][0]
    # Indexes to the closest clusters in XY
    ex_cls_idx = close_cl_idx[1][idx][1:]
    centers_ex = []
    for i in ex_cls_idx:
        ex_cl_dict = {'xy': [GDR3_f['GLON'][i], GDR3_f['GLAT'][i]]}
        if not np.isnan(GDR3_f['pmRA'][i]):
            ex_cl_dict['pms'] = [GDR3_f['pmRA'][i], GDR3_f['pmDE'][i]]
        if not np.isnan(GDR3_f['plx'][i]):
            ex_cl_dict['plx'] = [GDR3_f['plx'][i]]

        # print(GDR3_f['Name'][i], GDR3_f['GLON'][i], GDR3_f['GLAT'][i])

        centers_ex.append(ex_cl_dict)

    # ex_cl_dict = {'xy': [300.4, -16]}
    # centers_ex.append(ex_cl_dict)

    return centers_ex


def get_classif(lon, lat, pmRA, pmDE, plx, probs_final):
    """
    """
    # Select most probable members
    msk_membs = probs_final > 0.5

    xy_c = np.median([lon[msk_membs], lat[msk_membs]], 1)
    vpd_c = np.median([pmRA[msk_membs], pmDE[msk_membs]], 1)
    plx_c = np.median(plx[msk_membs])

    # Median distances to centers for members
    xy = np.array([lon[msk_membs], lat[msk_membs]]).T
    xy_rads = spatial.distance.cdist(xy, np.array([xy_c])).T[0]
    xy_05 = np.median(xy_rads)
    pm = np.array([pmRA[msk_membs], pmDE[msk_membs]]).T
    pm_rads = spatial.distance.cdist(pm, np.array([vpd_c])).T[0]
    pm_05 = np.median(pm_rads)
    plx_rad = abs(plx[msk_membs] - plx_c)
    plx_05 = np.median(plx_rad)
    # Count member stars within median distances
    N_memb_xy = (xy_rads < xy_05).sum()
    N_memb_pm = (pm_rads < pm_05).sum()
    N_memb_plx = (plx_rad < plx_05).sum()

    # Median distances to centers for field stars
    xy = np.array([lon[~msk_membs], lat[~msk_membs]]).T
    xy_rads_f = spatial.distance.cdist(xy, np.array([xy_c])).T[0]
    pm = np.array([pmRA[~msk_membs], pmDE[~msk_membs]]).T
    pm_rads_f = spatial.distance.cdist(pm, np.array([vpd_c])).T[0]
    plx_rad_f = abs(plx[~msk_membs] - plx_c)
    # Count field stars within median distances
    N_field_xy = (xy_rads_f < xy_05).sum()
    N_field_pm = (pm_rads_f < pm_05).sum()
    N_field_plx = (plx_rad_f < plx_05).sum()

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
        if plx_c > 2:
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


def write_out(out_path, cl_name, data, probs_all, N_membs, prob_min=0.5):
    """
    """
    data['probs'] = np.round(probs_all, 2)
    msk = probs_all > prob_min

    if N_membs > msk.sum():
        p_sort = np.sort(probs_all)[::-1]
        prob_min = p_sort[N_membs]
        msk = probs_all > prob_min

    df = data[msk]

    df.to_csv(out_path + cl_name + ".csv", index=False)


if __name__ == '__main__':
    main()

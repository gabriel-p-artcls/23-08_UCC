
import numpy as np
import pandas as pd
from scipy import spatial
import main_process_GDR3_query as G3Q

# LOCAL RUN
# insert at 1, 0 is the script path (or '' in REPL)
import sys
sys.path.insert(1, '/home/gabriel/Github/fastmp/')  # Path to fastMP
from fastmp import fastMP
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
# from fastmp import fastMP

# Maximum magnitude to retrieve
max_mag = 20


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

    for index, cl in database_full.iterrows():

        if 'KHARCHENKO12' in cl['DB'] or 'BICA19' in cl['DB'] or\
                'CANTAT20' in cl['DB'] or 'DIAS21' in cl['DB'] or\
                'LOKTIN17' in cl['DB']:
            pass
        else:
            continue

        fname0 = cl['fnames'].split(';')[0]

        if fname0 != 'ngc2516':
            continue
        if 'ubc' in fname0 or 'upk' in fname0 or 'fof' in fname0\
                or 'ufmg' in fname0:
            continue

        # Extract center coordinates
        xy_c, vpd_c, plx_c = (cl['GLON'], cl['GLAT']), None, None
        if not np.isnan(cl['pmRA']):
            vpd_c = (cl['pmRA'], cl['pmDE'])
        if not np.isnan(cl['plx']):
            plx_c = cl['plx']

        if vpd_c is None and plx_c is None:
            continue

        txt = 'nan'
        if 'LOKTIN17' in cl['DB']:
            txt = 'L17'

        # Get close clusters coords
        centers_ex = get_close_cls(
            index, database_full, close_cl_idx, cl['dups_fnames'])

        # Generate frame
        data = get_frame(frames_path, frames_data, cl)
        data.to_csv(fname0 + "_full.csv", index=False)
        breakpoint()

        # Generate input data array for fastMP
        X = np.array([
            data['GLON'].values, data['GLAT'].values, data['pmRA'].values,
            data['pmDE'].values, data['Plx'].values, data['e_pmRA'].values,
            data['e_pmDE'].values, data['e_Plx'].values])

        # Process with fastMP
        fixed_centers = True
        probs_all1 = fastMP(
            xy_c=xy_c, centers_ex=centers_ex,
            fixed_centers=fixed_centers).fit(X)
        probs_all2 = fastMP(
            xy_c=xy_c, vpd_c=vpd_c, plx_c=plx_c, centers_ex=centers_ex,
            fixed_centers=fixed_centers).fit(X)

        msk1 = probs_all1 > 0.5
        if msk1.sum() > 25:
            pmRA1 = np.nanmedian(data[msk1]['pmRA'])
            pmDE1 = np.nanmedian(data[msk1]['pmDE'])
        else:
            is1 = np.argsort(probs_all1)[::-1]
            pmRA1 = np.nanmedian(data['pmRA'].values[is1][:25])
            pmDE1 = np.nanmedian(data['pmDE'].values[is1][:25])
        d1 = np.sqrt((vpd_c[0]-pmRA1)**2 + (vpd_c[1]-pmDE1)**2)
        msk2 = probs_all2 > 0.5
        if msk2.sum() > 25:
            pmRA2 = np.nanmedian(data[msk2]['pmRA'])
            pmDE2 = np.nanmedian(data[msk2]['pmDE'])
        else:
            is2 = np.argsort(probs_all2)[::-1]
            pmRA2 = np.nanmedian(data['pmRA'].values[is2][:25])
            pmDE2 = np.nanmedian(data['pmDE'].values[is2][:25])
        d2 = np.sqrt((vpd_c[0]-pmRA2)**2 + (vpd_c[1]-pmDE2)**2)

        d3 = np.sqrt((pmRA1-pmRA2)**2 + (pmDE1-pmDE2)**2)

        print("*** {}, {}, {}, {}, {:.3f}, {:.3f}, {:.3f}".format(
              txt, fname0, msk1.sum(), msk2.sum(), d1, d2, d3))


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


if __name__ == '__main__':
    main()

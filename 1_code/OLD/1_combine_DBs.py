
import datetime
import numpy as np
import json
import csv
import pandas as pd
from scipy.spatial.distance import cdist
from astropy.coordinates import SkyCoord
from astropy.coordinates import angular_separation
import astropy.units as u

"""
Script used to combine all the  databases using the clusters' names
"""

dbs_folder = '/home/gabriel/Github/web_sites/UCC/datafiles/'
DBs_json = "all_dbs.json"
out_path = '../2_pipeline/'


def main(N_dups=10):
    """
    """
    dbs_used, DB_data, DB_names_match, DB_names_orig = get_data_and_names()

    unique_names, unique_names_orig = get_unique_names(
        DB_names_match, DB_names_orig)

    cl_dict = get_matches(
        dbs_used, DB_data, DB_names_match, unique_names, unique_names_orig)

    comb_dbs = combine_DBs(cl_dict)

    comb_dbs['ID'] = rm_dup_names(comb_dbs['ID'])

    comb_dbs['ID'], comb_dbs['fnames'] = assign_fname(comb_dbs['ID'])

    comb_dbs['UCC_ID'] = assign_UCC_ids(comb_dbs['GLON'], comb_dbs['GLAT'])
    # dup_check(ucc_ids, comb_dbs)

    comb_dbs['dups_fnames'] = dups_identify(comb_dbs, N_dups)

    print("Writing to file...")
    pd.DataFrame(comb_dbs).to_csv(
        out_path + '1_combined_DBs.csv', na_rep='nan', index=False,
        quoting=csv.QUOTE_NONNUMERIC)
    
    print("Combined database written to file")


def get_data_and_names():
    """
    1. For each DB extract and store all its data --> DB_data
    2. For each cluster in each DB extract and standardize its
       name(s) --> DB_names_match
    3. For each cluster in each DB extract its original unedited
       name(s) --> DB_names_orig
    """
    with open(dbs_folder + DBs_json) as f:
        dbs_used = json.load(f)

    DB_data, DB_names_match, DB_names_orig, N_all = {}, {}, {}, 0
    for DB, _ in dbs_used.items():
        # Load DB data
        data = pd.read_csv(
            dbs_folder + "databases/" + DB + ".csv", index_col=False)
        print(DB, len(data))
        N_all += len(data)
        # Store in dictionary
        DB_data[DB] = data

        # Extract and standardize all names
        names_all = data[dbs_used[DB]['names']]
        names_final, names_orig = [], []
        for i, names in enumerate(names_all):
            names_l, names_l_orig = [], []
            names_s = names.split(',')
            for name in names_s:
                name = cluster_rename(name)

                names_l.append(name.lower().replace('_', '').replace(
                               ' ', '').replace('-', '').replace("'", ''))
                names_l_orig.append(name)

            names_final.append(names_l)
            # Store original names
            names_orig.append(names_l_orig)

        DB_names_match[DB] = names_final
        DB_names_orig[DB] = names_orig

    print(f"\n{N_all} clusters in all DBs")

    return dbs_used, DB_data, DB_names_match, DB_names_orig


def cluster_rename(name):
    """
    Standardize the naming of these clusters watching for 

    FSR XXX w leading zeros
    FSR XXX w/o leading zeros
    FSR_XXX w leading zeros
    FSR_XXX w/o leading zeros

    ESO XXX-YY w leading zeros
    ESO XXX-YY w/o leading zeros
    ESO_XXX_YY w leading zeros
    ESO_XXX_YY w/o leading zeros
    ESO_XXX-YY w leading zeros
    ESO_XXX-YY w/o leading zeros
    """
    name = name.strip()

    if name.startswith("FSR"):
        if ' ' in name or '_' in name:
            if '_' in name:
                n2 = name.split('_')[1]
            else:
                n2 = name.split(' ')[1]
            n2 = int(n2)
            if n2 < 10:
                n2 = '000' + str(n2)
            elif n2 < 100:
                n2 = '00' + str(n2)
            elif n2 < 1000:
                n2 = '0' + str(n2)
            else:
                n2 = str(n2)
            name = "FSR_" + n2

    if name.startswith("ESO"):
        if ' ' in name[4:]:
            n1, n2 = name[4:].split(' ')
        elif '_' in name[4:]:
            n1, n2 = name[4:].split('_')
        elif '' in name[4:]:
            n1, n2 = name[4:].split('-')

        n1 = int(n1)
        if n1 < 10:
            n1 = '00' + str(n1)
        elif n1 < 100:
            n1 = '0' + str(n1)
        else:
            n1 = str(n1)
        n2 = int(n2)
        if n2 < 10:
            n2 = '0' + str(n2)
        else:
            n2 = str(n2)
        name = "ESO_" + n1 + '_' + n2

    return name


def get_unique_names(DB_names_match, DB_names_orig):
    """
    Identify unique names for all the DBs. An entry can store more than
    one unique name if the DB lists several names as belonging to the same
    cluster.

    If any cluster's name is identified with another name by any of the DBs
    (eg: Pismis 6 & NGC 2645),these two clusters *that could be shown as
    separate clusters by other DBs* will be merged into a single cluster
    in the final list.
    """
    all_names, all_names_orig = [], []
    for k, names_match in DB_names_match.items():
        all_names += names_match
        all_names_orig += DB_names_orig[k]

    print("Extracting unique names...")

    match_dict = {}
    for i, names_l in enumerate(all_names):
        uid = "cl" + str(i)

        clid_m = []
        for name in names_l:
            # Check if any of these names is already in dictionary
            for clid, d_names in match_dict.items():
                if name in d_names[1]:
                    clid_m += [clid]
            # Remove possible duplicate id
            clid_m = list(set(clid_m))

        if len(clid_m) == 0:
            match_dict[uid] = [[], []]
            match_dict[uid][0] = list(dict.fromkeys(all_names_orig[i]))
            match_dict[uid][1] = list(dict.fromkeys(names_l))
        elif len(clid_m) == 1:
            clid_m = clid_m[0]
            match_dict[clid_m][0] += all_names_orig[i]
            match_dict[clid_m][1] += names_l
            # Remove duplicates
            match_dict[clid_m][0] = list(dict.fromkeys(match_dict[clid_m][0]))
            match_dict[clid_m][1] = list(dict.fromkeys(match_dict[clid_m][1]))
        else:
            # Create new entry
            match_dict[uid] = [[], []]
            match_dict[uid][0] += all_names_orig[i]
            match_dict[uid][1] += names_l
            # Merge old entries into new entry
            for clid in clid_m:
                match_dict[uid][0] += list(match_dict[clid][0])
                match_dict[uid][1] += list(match_dict[clid][1])
                # Remove old entries from dictionary
                del match_dict[clid]
            # Remove duplicates
            match_dict[uid][0] = list(dict.fromkeys(match_dict[uid][0]))
            match_dict[uid][1] = list(dict.fromkeys(match_dict[uid][1]))

    unique_names, unique_names_orig = [], []
    for k, v in match_dict.items():
        unique_names_orig.append(v[0])
        unique_names.append(v[1])

    # # Check for duplicates
    # for i, cl in enumerate(unique_names):
    #     flag_match = False
    #     for cl0 in cl:
    #         for j, cl_rest in enumerate(unique_names[i + 1:]):
    #             for cl1 in cl_rest:
    #                 if cl0 == cl1:
    #                     flag_match = True
    #                     break
    #             if flag_match:
    #                 break
    #         if flag_match:
    #             break
    #     if flag_match:
    #         print(i, i + 1 + j, cl, unique_names[i + 1 + j])
    # print("end check")

    print(f"N={len(unique_names)} unique names identified")

    return unique_names, unique_names_orig


def get_matches(
    dbs_used, DB_data, DB_names_match, unique_names, unique_names_orig
):
    """
    If any name in the name lists stored in 'unique_names' matches any name
    in the lists of names in any DB, all those names are assumed to belong to
    the same unique cluster
    """
    print("Matching databases...")
    cl_dict = {}
    # For each list of unique names
    for q, unique_n in enumerate(unique_names):

        # For each name in list
        cl_str = ';'.join(unique_names_orig[q])
        cl_dict[cl_str] = {
            'DB': [], 'DB_i': [], 'RA': [], 'DE': [], 'plx': [], 'pmra': [],
            'pmde': []}

        # For each DB
        for DB_ID, names_db in DB_names_match.items():
            df = DB_data[DB_ID]
            cols = []
            for v in dbs_used[DB_ID]['pos'].split(','):
                if str(v) == 'None':
                    v = None
                cols.append(v)
            # Remove Rv column
            ra, de, plx, pmra, pmde = cols[:-1]

            # For each name in this list of unique names
            for name in unique_n:
                db_match = False

                # For each name list in this DB
                for i, name_l in enumerate(names_db):

                    # Unique name is in this list of names for this DB
                    if name in name_l:
                        db_match = True
                        cl_dict[cl_str]['DB'].append(DB_ID)
                        cl_dict[cl_str]['DB_i'].append(i)
                        # Extract row from this DB
                        row = df.iloc[i]
                        # if DB_ID != 'DIAS21':
                        #     # This DB has bad data for RA
                        cl_dict[cl_str]['RA'].append(row[ra])
                        cl_dict[cl_str]['DE'].append(row[de])
                        if DB_ID == 'KHARCHENKO12':
                            # This DB has bad data for these params
                            continue
                        if plx is not None:
                            cl_dict[cl_str]['plx'].append(row[plx])
                        if pmra is not None:
                            cl_dict[cl_str]['pmra'].append(row[pmra])
                        if pmde is not None:
                            cl_dict[cl_str]['pmde'].append(row[pmde])

                    if db_match:
                        break
                if db_match:
                    break

    return cl_dict


def combine_DBs(cl_dict):
    """
    Store unique values for each cluster
    """
    print("Generating final data...")

    db_l, db_i_l, names_l, ra_l, dec_l, glon_l, glat_l, plx_l, pmRA_l,\
        pmDE_l = [[] for _ in range(10)]
    for names, v in cl_dict.items():

        DBs, DBs_i, ra, dec, plx, pmRA, pmDE = v['DB'], v['DB_i'], v['RA'],\
            v['DE'], v['plx'], v['pmra'], v['pmde']

        # Store DBs and the indexes in them where the cluster is located
        db_l.append(";".join(str(_) for _ in DBs))
        db_i_l.append(";".join(str(_) for _ in DBs_i))

        names_l.append(names)

        ra_m = np.nanmedian(ra)
        dec_m = np.nanmedian(dec)
        ra_l.append(round(ra_m, 5))
        dec_l.append(round(dec_m, 5))

        lon, lat = radec2lonlat(ra_m, dec_m)
        glon_l.append(lon)
        glat_l.append(lat)

        if np.isnan(plx).all():
            plx_m = np.nan
        else:
            plx_m = np.nanmedian(plx)
        plx_l.append(round(plx_m, 3))

        if np.isnan(pmRA).all():
            pmRA_m = np.nan
        else:
            pmRA_m = np.nanmedian(pmRA)
        pmRA_l.append(round(pmRA_m, 3))

        if np.isnan(pmDE).all():
            pmDE_m = np.nan
        else:
            pmDE_m = np.nanmedian(pmDE)
        pmDE_l.append(round(pmDE_m, 3))

    # Store combined databases
    final_DB = {
        'DB': db_l, 'DB_i': db_i_l, 'ID': names_l, 'RA_ICRS': ra_l,
        'DE_ICRS': dec_l, 'GLON': glon_l, 'GLAT': glat_l, 'plx': plx_l,
        'pmRA': pmRA_l, 'pmDE': pmDE_l
    }

    return final_DB


def radec2lonlat(ra, dec):
    gc = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    lb = gc.transform_to('galactic')
    lon, lat = lb.l.value, lb.b.value
    return np.round(lon, 5), np.round(lat, 5)


def rm_dup_names(all_names):
    """
    This removes duplicates such as "NGC_2516" and "NGC 2516", "UBC 123" and
    "UBC123" and "OC-4567" and "OC 4567"
    """
    no_dup_names = []
    for names in all_names:
        names = names.split(';')
        names_temp = []
        for name in names:
            name = name.strip()
            if 'UBC' in name and 'UBC ' not in name and 'UBC_' not in name:
                name = name.replace('UBC', 'UBC ')
            if 'UBC_' in name:
                name = name.replace('UBC_', 'UBC ')

            if 'UFMG' in name and 'UFMG ' not in name and 'UFMG_' not in name:
                name = name.replace('UFMG', 'UFMG ')

            if 'LISC' in name and 'LISC ' not in name and 'LISC_' not in name\
                    and 'LISC-' not in name:
                name = name.replace('LISC', 'LISC ')

            if 'OC-' in name:
                name = name.replace('OC-', 'OC ')
            # Removes duplicates such as "NGC_2516" and "NGC 2516" 
            name = name.replace('_', ' ')
            names_temp.append(name)

        # Equivalent to set() but maintains order
        no_dup_names.append(';'.join(list(dict.fromkeys(names_temp))))

    return no_dup_names


def assign_fname(all_names):
    """
    Assign names used for files and urls
    """
    all_names_reorder, all_fnames_reorder = [], []
    for names in all_names:
        names = names.split(';')
        fnames_temp = []
        for i, name in enumerate(names):
            name = name.strip()
            # We replace '+' with 'p' to avoid duplicating names for clusters
            # like 'Juchert J0644.8-0925' and 'Juchert_J0644.8+0925'
            name = name.lower().replace('_', '').replace(' ', '').replace(
                '-', '').replace('.', '').replace("'", '').replace('+', 'p')
            fnames_temp.append(name)

        names_reorder, fnames_reorder = preferred_names(names, fnames_temp)

        all_names_reorder.append(";".join(names_reorder))
        # all_fnames_reorder.append(";".join(fnames_reorder))
        all_fnames_reorder.append(';'.join(list(dict.fromkeys(fnames_reorder))))

    return all_names_reorder, all_fnames_reorder


def preferred_names(names, fnames):
    """
    Use naming conventions according to this list of preferred names
    """
    names_lst = (
        'blanco', 'westerlund', 'ngc', 'melotte', 'trumpler', 'ruprecht',
        'berkeley', 'pismis', 'vdbh', 'loden', 'kronberger', 'collinder',
        'haffner', 'tombaugh', 'dolidze', 'auner', 'waterloo', 'basel',
        'bochum', 'hogg', 'carraro', 'lynga', 'johansson', 'mamajek',
        'platais', 'harvard', 'czernik', 'koposov', 'eso', 'ascc', 'teutsch',
        'alessi', 'king', 'saurer', 'fsr', 'juchert', 'antalova', 'stephenson')

    # Replace with another name according to the preference list
    if len(names) == 1:
        return names, fnames

    # Always move 'MWSC to the last position'
    if "mwsc" in fnames[0]:
        names = names[1:] + [names[0]]
        fnames = fnames[1:] + [fnames[0]]

    # # Select the first name listed
    def find_preferred_name(fnames):
        """Replace with another name according to the preference list"""
        for name_prefer in names_lst:
            for i, name in enumerate(fnames):
                if name_prefer in name:
                    return i
        return None

    def reorder_names(nms_lst, i):
        nms_reorder = list(nms_lst)
        name0 = nms_reorder[i]
        del nms_reorder[i]
        nms_reorder = [name0] + nms_reorder
        return nms_reorder

    i = find_preferred_name(fnames)
    if i is not None:
        # Reorder
        names_reorder = reorder_names(names, i)
        fnames_reorder = reorder_names(fnames, i)
    else:
        names_reorder, fnames_reorder = names, fnames

    return names_reorder, fnames_reorder


def assign_UCC_ids(glon, glat):
    """
    Format: UCC GXXX.X+YY.Y
    """
    lonlat = np.array([glon, glat]).T
    lonlat = trunc(lonlat)
    
    ucc_ids = []
    for idx, ll in enumerate(lonlat):
        lon, lat = str(ll[0]), str(ll[1])

        if ll[0] < 10:
            lon = '00' + lon
        elif ll[0] < 100:
            lon = '0' + lon

        if ll[1] >= 10:
            lat = '+' + lat
        elif ll[1] < 10 and ll[1] > 0:
            lat = '+0' + lat
        elif ll[1] == 0:
            lat = '+0' + lat.replace('-', '')
        elif ll[1] < 0 and ll[1] >= -10:
            lat = '-0' + lat[1:]
        elif ll[1] < -10:
            pass

        ucc_id = 'UCC G' + lon + lat

        i = 0
        while True:
            if i > 25:
                ucc_id += "ERROR"
                print("ERROR NAMING")
                break
            if ucc_id in ucc_ids:
                if i == 0:
                    # Add a letter to the end
                    ucc_id += ascii_lowercase[i]
                else:
                    # Replace last letter
                    ucc_id = ucc_id[:-1] + ascii_lowercase[i]
                i += 1
            else:
                break
        # if i > 0:
        #     print(ucc_id, dbs['ID'][idx], dbs['GLON'][idx], dbs['GLAT'][idx])

        ucc_ids.append(ucc_id)

    return ucc_ids


def trunc(values, decs=1):
    return np.trunc(values*10**decs)/(10**decs)


def dup_check(names, db):
    """
    Check for duplicates in 'names' list
    """
    for i, cl0 in enumerate(names):
        for j, cl1 in enumerate(names[i + 1:]):
            if cl0 == cl1:
                print(i, i + 1 + j, cl0, dbs['ID'][i], dbs['ID'][i + 1 + j])
                break
    print("End duplicates check")


def dups_identify(df, N_dups):
    """
    Find the closest clusters to all clusters
    """
    print("Estimating 2D (GLON, GLAT) distances...")
    x, y = df['GLON'], df['GLAT']
    pmRA, pmDE, plx = df['pmRA'], df['pmDE'], df['plx']
    coords = np.array([x, y]).T
    # Find the distances to all clusters, for all clusters
    dist = cdist(coords, coords)
    # Change distance to itself from 0 to inf
    msk = dist == 0.
    dist[msk] = np.inf

    print(f"Finding duplicates (max={N_dups})...")
    dups_fnames = []
    for i, cl in enumerate(dist):
        idx = np.argsort(cl)[:N_dups]

        dups_fname = []
        for j in idx:
            # Angular distance in arcmin (rounded)
            d = round(angular_separation(x[i], y[i], x[j], y[j]) * 60, 2)
            # PMs distance
            pm_d = np.sqrt((pmRA[i]-pmRA[j])**2 + (pmDE[i]-pmDE[j])**2)
            # Parallax distance
            plx_d = abs(plx[i] - plx[j])

            dup_flag = duplicate_find(d, pm_d, plx_d, plx[i])

            if dup_flag:
                fname = df['fnames'][j].split(';')[0]
                dups_fname.append(fname)

        if dups_fname:
            print(i, df['fnames'][i], len(dups_fname), dups_fname)
            dups_fname = ";".join(dups_fname)
        else:
            dups_fname = 'nan'

        dups_fnames.append(dups_fname)

    return dups_fnames


def duplicate_find(d, pm_d, plx_d, plx):
    """
    Identify a cluster as a duplicate following an arbitrary definition
    that depends on the parallax
    """
    if plx >= 4:
        rad, plx_r, pm_r = 15, 0.5, 1
    elif 3 <= plx and plx < 4:
        rad, plx_r, pm_r = 10, 0.25, 0.5
    elif 2 <= plx and plx < 3:
        rad, plx_r, pm_r = 5, 0.15, 0.25
    elif 1 <= plx and plx < 2:
        rad, plx_r, pm_r = 2.5, 0.1, 0.15
    else:
        rad, plx_r, pm_r = 1, 0.05, 0.1

    if pm_d < pm_r and plx_d < plx_r and d < rad:
        return True

    return False


if __name__ == '__main__':
    # plt.style.use('science')
    main()


import numpy as np
import json
import csv
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

"""
Script used to combine all the  databases using the clusters' names
"""

dbs_folder = '/home/gabriel/Github/web_sites/UCC/datafiles/'
DBs_json = "all_dbs.json"
out_path = '../2_pipeline/'


def main():
    """
    """
    dbs_used, DB_data, DB_names_match, DB_names_orig = get_data_and_names()

    unique_names, unique_names_orig = get_unique_names(
        DB_names_match, DB_names_orig)

    cl_dict = get_matches(
        dbs_used, DB_data, DB_names_match, unique_names, unique_names_orig)

    final_DB = combine_DBs(cl_dict)

    print("Writing to file...")
    pd.DataFrame(final_DB).to_csv(
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
                               ' ', '').replace('-', ''))
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


def lonlat2radec(lon, lat):
    gc = SkyCoord(l=lon * u.degree, b=lat * u.degree, frame='galactic')
    rd = gc.transform_to('fk5')
    ra, dec = rd.ra.value, rd.dec.value
    return np.round(ra, 5), np.round(dec, 5)


def radec2lonlat(ra, dec):
    gc = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    lb = gc.transform_to('galactic')
    lon, lat = lb.l.value, lb.b.value
    return np.round(lon, 5), np.round(lat, 5)


if __name__ == '__main__':
    # plt.style.use('science')
    main()

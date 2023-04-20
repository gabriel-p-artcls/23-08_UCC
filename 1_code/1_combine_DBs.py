
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

"""
Script used to combine all the  databases using the clusters' names

Output (row example):

UCC_ID: UCC GXXX.X-YY.Y
DBs: "1469"
name: cluster
all_names: "name1, name2, etc"
ra, dec, glon, glat: "val1, val2, ..."
plx, pmRa, pmDE: "val1, val2, ..."

"""

dbs_folder = '/home/gabriel/Github/web_sites/UCC/datafiles/databases/'
out_path = '../2_pipeline/'

# Data stored per cluster (per database):
# Name(s), RA, DEC, plx, pmRA, pmDE

dbs_used = {
    'KHARCHENKO12': ['name', 'ra', 'dec', None, 'pm_ra', 'pm_dec'],
    'CASTRO18': ["Cluster", "_RA.icrs", "_DE.icrs", "Plx", "pmRA*", "pmDE"],
    'BICA19': ['Name', 'RA_ICRS', 'DE_ICRS', None, None, None],
    'CASTRO19': ["Cluster", "RA_ICRS", "DE_ICRS", "Plx", "pmRA*", "pmDE"],
    'SIM19': ["ID", "RA_ICRS", "DE_ICRS", "plx", "pmRA*", "pmDE"],
    'LIUPANG19': ["ID", "_RA.icrs", "_DE.icrs", "plx", "pmRA", "pmDE"],
    'FERREIRA19': ["Name", "RA_ICRS", "DE_ICRS", None, "pmRA*", "pmDE"],
    'CANTAT20': ["Cluster", "RA_ICRS", "DE_ICRS", "plx", "pmRA*", "pmDE"],
    'CASTRO20': ["Cluster", "RA_ICRS", "DE_ICRS", "plx", "pmRA", "pmDE"],
    'FERREIRA20': ["Name", "RA_ICRS", "DE_ICRS", "plxcl", "pmRAcl", "pmDEcl"],
    'HAO20': ["Name", "ra", "dec", 'plx', "pmra", "pmde"],
    'DIAS21': ['Cluster', 'RA_ICRS', 'DE_ICRS', 'Plx', 'pmRA', 'pmDE'],
    'CASADO21': ["Name", "RA_ICRS", "DE_ICRS", "plx", "pmra", "pmde"],
    'FERREIRA21': ["Name", "RAJ2000", "DEJ2000", "plx", "pmRA", "pmDE"],
    'HUNT21': ['Name', 'RA_ICRS', 'DE_ICRS', 'plx', 'pmRA', 'pmDE'],
    'JAEHNIG21': ["Name", "ra", "dec", 'plx', "pmra", "pmde"],
    'SANTOS21': ["Name", "ra", "dec", 'plx', "pmra", "pmde"],
    'HE21': ["OC", "RA_ICRS", "DE_ICRS", "plx", "pmRA", "pmDE"],
    'CASTRO22': ["Cluster", "RA_ICRS", "DE_ICRS", "plx", "pmRA", "pmDE"],
    'TARRICQ22': ["Cluster", "RA_ICRS", "DE_ICRS", "plx", "pmRA", "pmDE"],
    'LI22': ["id", "ra", "dec", 'plx', "pmRA", "pmDE"],
    'HE22': ["Cluster", "_RA.icrs", "_DE.icrs", "Plx", "pmRA", "pmDE"],
    'HE22_1': ["Cluster", "_RA.icrs", "_DE.icrs", "plx", "pmRA", "pmDE"],
    'HE22_2': ["CWNU_id", "RA_ICRS", "DE_ICRS", "plx", "pmx", "pmy"],
    'HAO22': ["Cluster", "RA_ICRS", "DE_ICRS", 'plx', "pmRA", "pmDE"],
    'HUNT23': ['name', 'ra', 'dec', 'parallax', 'pmra', 'pmdec'],
    'QIN23': ["Name", "RAdeg", "DEdeg", 'plx', "pmRA", "pmDE"],
    'LI23': ["id", "ra", "dec", 'plx', "pmra", "pmdec"],
    'CHI23_2': ["Name", "ra", "dec", 'plx', 'pmra', "pmde"],
    'CHI23': ["name", "RA_ICRS", "DE_ICRS", 'plx', 'pmra', "pmdec"],
    'CHI23_3': ["Id", "ra", "dec", 'plx', 'pmra', "pmdec"],
}

DBs_IDS = {
    'KHARCHENKO12': 1, 'CASTRO18': 2, 'BICA19': 3, 'CASTRO19': 4, 'SIM19': 5,
    'LIUPANG19': 6, 'FERREIRA19': 7, 'CANTAT20': 8, 'CASTRO20': 9,
    'FERREIRA20': 10, 'HAO20': 11, 'DIAS21': 12, 'CASADO21': 13,
    'FERREIRA21': 14, 'HUNT21': 15, 'JAEHNIG21': 16, 'SANTOS21': 17,
    'HE21': 18, 'CASTRO22': 19, 'TARRICQ22': 20, 'LI22': 21, 'HE22': 22,
    'HE22_1': 23, 'HE22_2': 24, 'HAO22': 25, 'HUNT23': 26, 'QIN23': 27,
    'LI23': 28, 'CHI23_2': 29, 'CHI23': 30, 'CHI23_3': 31,
}


def main():
    """
    """
    # data = pd.read_csv(dbs_folder + 'CHI23' + ".csv", index_col=False)
    # ra, dec = lonlat2radec(data['l'].values, data['b'].values)
    # data['RA_ICRS'], data['DE_ICRS'] = ra, dec
    # # data['plx'] = 1000. / data['dist_pc'].values
    # data.to_csv(dbs_folder+'CHI23_3.csv', index=False)
    # breakpoint()

    # # # np.random.seed(12345)
    # keys = list(dbs_used.keys())
    # np.random.shuffle(keys)
    # dbs_used2 = {}
    # for k in keys:
    #     dbs_used2[k] = dbs_used[k]

    DB_data, DB_names_match, DB_names_orig = get_data_and_names(
        dbs_used, DBs_IDS)

    unique_names, unique_names_orig = get_unique_names(
        DB_names_match, DB_names_orig)

    cl_dict = get_matches(
        dbs_used, DB_data, DB_names_match, unique_names, unique_names_orig)

    final_DB = combine_DBs(cl_dict, DBs_IDS)

    print("Writing to file...")
    pd.DataFrame(final_DB).to_csv(
        out_path + '1_combined_DBs.csv', na_rep='nan', index=False)
    
    print("Combined database written to file")


def get_data_and_names(dbs_used, DBs_IDS):
    """
    1. For each DB extract and store all its data --> DB_data
    2. For each cluster in each DB extract and standardize its
       name(s) --> DB_names_match
    3. For each cluster in each DB extract its original unedited
       name(s) --> DB_names_orig
    """
    DB_data, DB_names_match, DB_names_orig, N_all = {}, {}, {}, 0
    for DB, _ in dbs_used.items():
        # Load DB data
        data = pd.read_csv(dbs_folder + DB + ".csv", index_col=False)
        print(DBs_IDS[DB], DB, len(data))
        N_all += len(data)
        # Store in dictionary
        DB_data[DB] = data

        # Extract and standardize all names
        names_all = data[dbs_used[DB][0]]
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

    return DB_data, DB_names_match, DB_names_orig


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
    """
    print("Matching databases...")
    cl_dict = {}
    # For each list of unique names
    for q, cl in enumerate(unique_names):

        # For each name in list
        cl_str = ','.join(unique_names_orig[q])
        cl_dict[cl_str] = {
            'DB': [], 'DB_i': [], 'RA': [], 'DE': [], 'plx': [], 'pmra': [],
            'pmde': []}

        # For each DB
        for DB_ID, names_db in DB_names_match.items():
            df, cols = DB_data[DB_ID], dbs_used[DB_ID]
            ra, de, plx, pmra, pmde = cols[1:]

            # For each name in this list of unique names
            for name in cl:
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


def combine_DBs(cl_dict, DBs_IDS):
    """
    Store unique values for each cluster
    """
    print("Generating final data...")

    db_l, db_i_l, names_l, ra_l, dec_l, glon_l, glat_l, plx_l, pmRA_l,\
        pmDE_l = [[] for _ in range(10)]
    for names, v in cl_dict.items():

        DBs, DBs_i, ra, dec, plx, pmRA, pmDE = v['DB'], v['DB_i'], v['RA'],\
            v['DE'], v['plx'], v['pmra'], v['pmde']

        in_dbs = [DBs_IDS[_] for _ in DBs]
        db_l.append("_".join(str(_) for _ in in_dbs))
        db_i_l.append("_".join(str(_) for _ in DBs_i))

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


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

dbs_folder = '../0_data/databases/'

# Data stored per cluster (per database):
# Name(s), RA, DEC, plx, pmRA, pmDE

dbs_used = {
    'BICA19': ['Name', 'RA_ICRS', 'DE_ICRS', None, None, None],
    'KHARCHENKO12': ['name', 'ra', 'dec', None, 'pm_ra', 'pm_dec'],
    'CANTAT20': ["Cluster", "RA_ICRS", "DE_ICRS", "plx", "pmRA*", "pmDE"],
    'DIAS21': ['Cluster', 'RA_ICRS', 'DE_ICRS', 'Plx', 'pmRA', 'pmDE'],
    'HUNT21': ['Name', 'RA_ICRS', 'DE_ICRS', 'plx', 'pmRA', 'pmDE'],
    'HUNT23': ['name', 'ra', 'dec', 'parallax', 'pmra', 'pmdec'],
    'TARRICQ22': ["Cluster", "RA_ICRS", "DE_ICRS", "plx", "pmRA", "pmDE"],
    'CASADO21': ["Name", "RA_ICRS", "DE_ICRS", "plx", "pmra", "pmde"],
    'CASTRO18': ["Cluster", "_RA.icrs", "_DE.icrs", "Plx", "pmRA*", "pmDE"],
    'CASTRO19': ["Cluster", "RA_ICRS", "DE_ICRS", "Plx", "pmRA*", "pmDE"],
    'CASTRO20': ["Cluster", "RA_ICRS", "DE_ICRS", "plx", "pmRA", "pmDE"],
    'CASTRO22': ["Cluster", "RA_ICRS", "DE_ICRS", "plx", "pmRA", "pmDE"],
    'SIM19': ["ID", "RA_ICRS", "DE_ICRS", "plx", "pmRA*", "pmDE"],
    'LIUPANG19': ["ID", "_RA.icrs", "_DE.icrs", "plx", "pmRA", "pmDE"],
    'FERREIRA19': ["Name", "RA_ICRS", "DE_ICRS", None, "pmRA*", "pmDE"],
    'FERREIRA20': ["Name", "RA_ICRS", "DE_ICRS", "plxcl", "pmRAcl", "pmDEcl"],
    'FERREIRA21': ["Name", "RAJ2000", "DEJ2000", "plx", "pmRA", "pmDE"],
    'HE21': ["OC", "RA_ICRS", "DE_ICRS", "plx", "pmRA", "pmDE"],
    'HE22': ["Cluster", "_RA.icrs", "_DE.icrs", "Plx", "pmRA", "pmDE"],
    'HE22_1': ["Cluster", "_RA.icrs", "_DE.icrs", "plx", "pmRA", "pmDE"],
    'HE22_2': ["CWNU_id", "RA_ICRS", "DE_ICRS", "plx", "pmx", "pmy"],
    'HAO20': ["Name", "ra", "dec", 'plx', "pmra", "pmde"],
    'HAO22': ["Cluster", "RA_ICRS", "DE_ICRS", 'plx', "pmRA", "pmDE"],
    'LI22': ["id", "ra", "dec", 'plx', "pmRA", "pmDE"],
    'QIN23': ["Name", "RAdeg", "DEdeg", 'plx', "pmRA", "pmDE"],
    'JAEHNIG21': ["Name", "ra", "dec", 'plx', "pmra", "pmde"],
    'SANTOS21': ["Name", "ra", "dec", 'plx', "pmra", "pmde"],
    'LI23': ["id", "ra", "dec", 'plx', "pmra", "pmdec"],
    'CHI23_2': ["Name", "ra", "dec", 'plx', 'pmra', "pmde"],
    'CHI23': ["name", "RA_ICRS", "DE_ICRS", 'plx', 'pmra', "pmdec"],
    'CHI23_3': ["Id", "ra", "dec", 'plx', 'pmra', "pmdec"],
}

DBs_IDS = {}
i = 0
for db in dbs_used.keys():
    DBs_IDS[db] = i + 1
    i += 1


def main():
    """
    """

    # data = pd.read_csv(dbs_folder + 'CHI23' + ".csv", index_col=False)
    # ra, dec = lonlat2radec(data['l'].values, data['b'].values)
    # data['RA_ICRS'], data['DE_ICRS'] = ra, dec
    # # data['plx'] = 1000. / data['dist_pc'].values
    # data.to_csv(dbs_folder+'CHI23_3.csv', index=False)
    # breakpoint()

    DB_data, DB_names_match, DB_names = get_data_and_names(dbs_used)

    unique_names, unique_names_orig = get_unique_names(
        DB_names_match, DB_names)

    cl_dict = get_matches(
        DB_data, DB_names_match, unique_names, unique_names_orig)

    final_DB = combine_DBs(cl_dict)

    print("Writing to file...")
    pd.DataFrame(final_DB).to_csv('final_DB.csv', na_rep='nan', index=False)
    
    print("Finished")


def get_data_and_names(dbs_used):
    """
    """
    DB_data, DB_names_match, DB_names = {}, {}, {}
    for DB, cols in dbs_used.items():
        # Load DB data
        data = pd.read_csv(dbs_folder + DB + ".csv", index_col=False)
        print(DBs_IDS[DB], DB, len(data))
        # Store in dictionary
        DB_data[DB] = data

        # Extract and standardize all names
        names_all = data[dbs_used[DB][0]]
        names_final, names_orig = [], []
        for i, names in enumerate(names_all):
            names_l, names_l_orig = [], []
            names_s = names.split(',')
            for name in names_s:
                name = name.strip()

                if name.startswith("FSR"):
                    _, n2 = name.split('FSR')
                    # Remove the leading '_' for CANTAT20 & DIAS21
                    if '_' in n2:
                        name = "FSR_" + str(int(n2[1:]))
                    else:
                        name = "FSR_" + str(int(n2))

                if name.startswith("ESO"):
                    if '-' in name and '_' in name:
                        # This is a HUNT23 cluster
                        _, n2 = name.split('ESO')
                        n3, n4 = n2[1:].split('-')
                    else:
                        _, n2 = name.split('ESO')
                        try:
                            n3, n4 = n2.split('-')
                        except ValueError:
                            # Remove the leading '_' for CANTAT20 & DIAS21
                            n3, n4 = n2[1:].split('_')
                    name = "ESO_" + str(int(n3)) + '_' + str(int(n4))

                names_l.append(name.lower().replace('_', '').replace(
                               ' ', '').replace('-', ''))
                names_l_orig.append(name)

            names_final.append(names_l)
            names_orig.append(names_l_orig)

        DB_names_match[DB] = names_final
        DB_names[DB] = names_orig

    return DB_data, DB_names_match, DB_names


def get_unique_names(DB_names_match, DB_names):
    """
    """
    all_names, all_names_orig = [], []
    for k, v in DB_names_match.items():
        all_names += v
        v_orig = DB_names[k]
        all_names_orig += v_orig

    print("Extracting unique names...")
    unique_names, unique_names_orig, idx_match = [], [], []
    N_all, pp = len(all_names), 10
    for i, cl in enumerate(all_names):
        if i in idx_match:
            # If cluster was already matched, skip
            continue
        if 100 * (i / N_all) > pp:
            print(pp, "%")
            pp += 10
        name_matches = list(cl)
        # For each name for this cluster
        for name in cl:
            # For each remaining list of cluster names
            for j, name_l in enumerate(all_names[i + 1:]):
                # If this name is in this list
                if name in name_l:
                    # Store the index of the match
                    j_match = j + i + 1
                    idx_match.append(j_match)
                    # Add the cluster name(s) to the list
                    name_matches += name_l
        unique_names.append(list(set(name_matches)))
        unique_names_orig.append(all_names_orig[i])

    print(f"N={len(unique_names)} unique names identified")
    return unique_names, unique_names_orig


def get_matches(DB_data, DB_names_match, unique_names, unique_names_orig):
    """
    """
    print("Matching databases...")
    cl_dict = {}
    # For each list of names
    N_all, pp = len(unique_names), 10
    for q, cl in enumerate(unique_names):
        if 100 * (q / N_all) > pp:
            print(pp, "%")
            pp += 10

        # For each name in list
        cl_str = ','.join(unique_names_orig[q])
        cl_dict[cl_str] = {
            'DB': [], 'RA': [], 'DE': [], 'plx': [], 'pmra': [], 'pmde': []}
        
        # For each DB
        for DB_ID, names_db in DB_names_match.items():
            df, cols = DB_data[DB_ID], dbs_used[DB_ID]
            ra, de, plx, pmra, pmde = cols[1:]
            for name in cl:
                # For each name list in DB
                db_match = False
                for i, name_l in enumerate(names_db):
                    # Match found
                    if name in name_l:
                        db_match = True
                        cl_dict[cl_str]['DB'].append(DB_ID)
                        # Extract data
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

    db_l, names_l, ra_l, dec_l, glon_l, glat_l, plx_l, pmRA_l, pmDE_l =\
        [[] for _ in range(9)]
    for names, v in cl_dict.items():

        DBs, ra, dec, plx, pmRA, pmDE = v['DB'], v['RA'], v['DE'], v['plx'],\
            v['pmra'], v['pmde']

        in_dbs = [DBs_IDS[_] for _ in DBs]
        db_l.append("_".join(str(_) for _ in in_dbs))

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
        'DB': db_l, 'ID': names_l, 'RA_ICRS': ra_l, 'DE_ICRS': dec_l,
        'GLON': glon_l, 'GLAT': glat_l, 'plx': plx_l, 'pmRA': pmRA_l,
        'pmDE': pmDE_l
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

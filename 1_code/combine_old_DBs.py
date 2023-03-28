
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

"""
Script used to combine the four 'old' databases (CG20, DIAS21, BICA19, MWSC)
using the clusters' names (ie: not their coordinates).

DIAS21 has bad RA values for some clusters, discarded.
"""

dbs_folder = '../0_data/databases/old/'

# Data stored per cluster (per database):
# ID, RA, DEC, plx, pmRA, pmDE,
# Dist, AV, Age, [Fe/H] <-- Removed, data will be used elsewhere
N_vals = 5

dbs_names = {
    'DIAS21': ['Cluster', 'RA_ICRS', 'DE_ICRS', 'Plx', 'pmRA', 'pmDE'],
    # 'Dist', 'Av', 'logage', '[Fe/H]'],
    'CG20': ["Name", "RA", "DEC", "plx", "pmRA", "pmDE"],
    # 'DistPc', 'AVNN', 'AgeNN', None],
    'MWSC': ['name', 'ra', 'dec', None, 'pm_ra', 'pm_dec'],
    # 'distance_modulus', 'e_bv', 'log_age', 'metallicity'],
    'TARRICQ22': ["Cluster", "RA_ICRS", "DE_ICRS", "plx", "pmRA", "pmDE"],
    # "Dist", None, "logAgeNN", None],
    'HAO21': ["Cluster", "RA_ICRS", "DE_ICRS", "plx", "pmRA*", "pmDE"],
    # None, None, "logt", None],
    'BICA19_all_names': ['Name', 'GLON', 'GLAT', None, None, None],
    # None, None, None, None]
}

# DIAS21, CG20, MWSC, TARRICQ22, BICA19, HAO21
DBs_IDS = np.array([1, 2, 3, 4, 5, 6])


def main():
    """
    """

    DB_data, DB_names_match = {}, {}
    for DB, cols in dbs_names.items():
        data = pd.read_csv(
            dbs_folder + DB + ".csv", sep=',', comment='#', index_col=False)

        print(DB, len(data))

        DB_data[DB] = data
        names = data[dbs_names[DB][0]]
        names = [_.lower().replace('_', '').replace(' ', '').replace('-', '')
                 for _ in names]
        DB_names_match[DB] = np.array(names)

    # Create initial dictionary
    DB_id1, DB_id2 = 'None', 'DIAS21'
    cl_dict = match_DBs(
        {}, dbs_names, DB_names_match, DB_data, DB_id1, DB_id2)

    DB_id1, DB_id2 = 'DIAS21', 'CG20'
    cl_dict = match_DBs(
        cl_dict, dbs_names, DB_names_match, DB_data, DB_id1, DB_id2)
    print(f"Clusters in cross-matched catalog so far: {len(cl_dict)}")

    DB_id1 = DB_id1 + "+" + DB_id2
    DB_id2 = 'MWSC'
    cl_dict = match_DBs(
        cl_dict, dbs_names, DB_names_match, DB_data, DB_id1, DB_id2)
    print(f"Clusters in cross-matched catalog so far: {len(cl_dict)}")

    DB_id1 = DB_id1 + "+" + DB_id2
    DB_id2 = 'TARRICQ22'
    cl_dict = match_DBs(
        cl_dict, dbs_names, DB_names_match, DB_data, DB_id1, DB_id2)
    print(f"Clusters in cross-matched catalog so far: {len(cl_dict)}")

    DB_id1 = DB_id1 + "+" + DB_id2
    DB_id2 = 'HAO21'
    cl_dict = match_DBs(
        cl_dict, dbs_names, DB_names_match, DB_data, DB_id1, DB_id2)
    print(f"Clusters in cross-matched catalog so far: {len(cl_dict)}")

    DB_id1 = DB_id1 + "+" + DB_id2
    DB_id2 = 'BICA19_all_names'
    cl_dict = match_BICA19(cl_dict, dbs_names, DB_data, DB_id1, DB_id2)
    print(f"Clusters in final cross-matched catalog: {len(cl_dict)}")

    final_DB = combine_DBs(cl_dict)

    pd.DataFrame(final_DB).to_csv('OLD_DBs.csv', na_rep='nan', index=False)
    
    print("Finished")


def DIAS21_CG20(cl_dict, dbs_names, DB_names_match, DB_data, DB_id1, DB_id2):
    """
    Find DIAS21 clusters in CG20
    """
    cols1, cols2 = dbs_names[DB_id1], dbs_names[DB_id2]
    db2_idxs_match = []
    for i, cl in enumerate(DB_names_match[DB_id1]):
        # DIAS21
        full_name = DB_data[DB_id1][cols1[0]][i]
        # Save cluster data
        cl_dict[full_name] = [[] for _ in range(N_vals)]

        for q, col in enumerate(cols1[1:]):
            if col is not None:
                cl_dict[full_name][q].append(DB_data[DB_id1][col][i])
            else:
                cl_dict[full_name][q].append(np.nan)

        # CG20
        j = np.where(DB_names_match[DB_id2] == cl)[0]
        if len(j) > 0:
            # Match found
            j = j[0]
            db2_idxs_match.append(j)
            for q, col in enumerate(cols2[1:]):
                if col is not None:
                    cl_dict[full_name][q].append(DB_data[DB_id2][col][j])
                else:
                    cl_dict[full_name][q].append(np.nan)
        else:
            # No match for this cluster in DIAS21
            # print("CL in DIAS21 no match in CG20", i, full_name)
            for q, col in enumerate(cols2[1:]):
                cl_dict[full_name][q].append(np.nan)

    print(f"\nClusters in {DB_id1} matched to clusters in {DB_id2}: "
          + f"{len(db2_idxs_match)}")
    print(f"Clusters in {DB_id1} not matched to clusters in {DB_id2}: "
          + f"{len(DB_data[DB_id1]) - len(db2_idxs_match)}")

    cl_dict = store_no_match(
        dbs_names, DB_data, cl_dict, DB_id1, DB_id2, db2_idxs_match)

    return cl_dict


def match_DBs(
        cl_dict, dbs_names, DB_names_match, DB_data, DB_id1, DB_id2
):
    """
    Match DB1 with DB2
    """
    cl_dict, db2_idxs_match = store_match(
        DB_data, DB_names_match, cl_dict, DB_id1, DB_id2)

    cl_dict = store_no_match(
        dbs_names, DB_data, cl_dict, DB_id1, DB_id2, db2_idxs_match)

    return cl_dict


def store_match(DB_data, DB_names_match, cl_dict, DB_id1, DB_id2):
    """
    """
    N_old = len(cl_dict)
    cols3 = dbs_names[DB_id2]
    db2_idxs_match = []
    for i, (cl, vals) in enumerate(cl_dict.items()):

        cl_match = cl.lower().replace('_', '').replace(' ', '').replace('-', '')
        j = np.where(DB_names_match[DB_id2] == cl_match)[0]
        if len(j) > 0:
            # Match found
            j = j[0]
            db2_idxs_match.append(j)
            for q, col in enumerate(cols3[1:]):
                if col is not None:
                    cl_dict[cl][q].append(DB_data[DB_id2][col][j])
                else:
                    cl_dict[cl][q].append(np.nan)
        else:
            # No match for this cluster in DB2
            # print(f"CL in {DB_id1} no match in {DB_id2}", i, cl_match)
            for q, col in enumerate(cols3[1:]):
                cl_dict[cl][q].append(np.nan)

    print(f"\nClusters in {DB_id1} matched to clusters in {DB_id2}: "
          + f"{len(db2_idxs_match)}")
    print(f"Clusters in {DB_id1} not matched to clusters in {DB_id2}: "
          + f"{N_old - len(db2_idxs_match)}")

    return cl_dict, db2_idxs_match


def store_no_match(
        dbs_names, DB_data, cl_dict, DB_id1, DB_id2, db2_idxs_match
):
    """
    Store clusters from DB1 with no match in DB2
    """
    N_nans = len(DB_id1.split('+'))

    cols2 = dbs_names[DB_id2]
    DB2_all_idxs = np.arange(0, len(DB_data[DB_id2]))
    DB2_no_match = list(set(DB2_all_idxs) - set(db2_idxs_match))
    print(f"Clusters in {DB_id2} not matched to clusters in {DB_id1}: "
          + f"{len(DB2_no_match)}")
    for i in DB2_no_match:
        if DB_id2 == 'BICA19_all_names':
            full_name = DB_data[DB_id2][cols2[0]][i].split('|')[0]
        else:
            full_name = DB_data[DB_id2][cols2[0]][i]
        # print("CL in DB2 no match in DB1", i, full_name)
        # Save cluster data
        cl_dict[full_name] = [
            [np.nan for _ in range(N_nans)] for _ in range(N_vals)]
        for q, col in enumerate(cols2[1:]):
            if col is not None:
                if DB_id1 == 'None':
                    cl_dict[full_name][q] = [DB_data[DB_id2][col][i]]
                else:
                    cl_dict[full_name][q].append(DB_data[DB_id2][col][i])
            else:
                if DB_id1 == 'None':
                    cl_dict[full_name][q] = [np.nan]
                else:
                    cl_dict[full_name][q].append(np.nan)

    return cl_dict


def match_BICA19(cl_dict, dbs_names, DB_data, DB_id1, DB_id2):
    """
    Match DIAS21+CG20+MWSC with BICA19
    """
    N_old = len(cl_dict)
    cols4 = dbs_names[DB_id2]

    DB_names_match = []
    for names in DB_data[DB_id2][cols4[0]]:
        all_names = names.split('|')
        nm = [_.lower().replace('_', '').replace(' ', '').replace('-', '')
              for _ in all_names]
        DB_names_match.append(nm)

    db2_idxs_match = []
    for i, (cl, vals) in enumerate(cl_dict.items()):
        cl_match = cl.lower().replace('_', '').replace(' ', '').replace(
            '-', '')
        j = -1
        for jj, nl in enumerate(DB_names_match):
            if cl_match in nl:
                j = jj
                break
        if j > -1:
            # Match found
            db2_idxs_match.append(j)
            for q, col in enumerate(cols4[1:]):
                if col is not None:
                    val = DB_data[DB_id2][col][j]
                    cl_dict[cl][q].append(val)
                else:
                    cl_dict[cl][q].append(np.nan)
        else:
            # No match for this cluster in BICA19
            for q, col in enumerate(cols4[1:]):
                cl_dict[cl][q].append(np.nan)

    print(f"\nClusters in {DB_id1} matched to clusters in {DB_id2}: "
          + f"{len(db2_idxs_match)}")
    print(f"Clusters in {DB_id1} not matched to clusters in {DB_id2}: "
          + f"{N_old - len(db2_idxs_match)}")

    cl_dict = store_no_match(
        dbs_names, DB_data, cl_dict, DB_id1, DB_id2, db2_idxs_match)

    return cl_dict


def combine_DBs(cl_dict):
    """
    Store unique values for each cluster
    """
    print("\nWriting to file...")

    db_l, names_l, ra_l, dec_l, glon_l, glat_l, plx_l, pmRA_l, pmDE_l, dmod_l,\
        extin_l, logAge_l, FeH_l = [[] for _ in range(N_vals + 1 + 3)]
    for cl, vals in cl_dict.items():
        names_l.append(cl)
        ra, dec, plx, pmRA, pmDE, Dist, extin, logAge, FeH = vals

        msk = ~np.isnan(ra)
        db_l.append("".join(str(_) for _ in DBs_IDS[msk]))

        # Convert BICA19 values
        lon, lat = ra[-1], dec[-1]
        ra_b, dec_b = lonlat2radec(lon, lat)

        # These clusters are only present in DIAS21
        if msk.sum() == 1 and bool(msk[0]) is True:
            ra_m = ra[0]
        else:
            # Don't use DIAS21's RA
            ra_m = np.nanmedian([ra[1], ra[2], ra_b])
        dec_m = np.nanmedian([dec[0], dec[1], dec[2], dec_b])
        ra_l.append(round(ra_m, 5))
        dec_l.append(round(dec_m, 5))

        lon, lat = radec2lonlat(ra_m, dec_m)
        glon_l.append(lon)
        glat_l.append(lat)

        if np.isnan(plx).all():
            plx_m = np.nan
        elif not np.isnan(plx[1]):
            # Prefer CG20 values whenever possible
            plx_m = plx[1]
        else:
            # # Avoid HAO21 Plx data if possible
            # if np.isnan(plx[1:]).all():
            #     plx_m = plx[0]
            # else:
            #     plx_m = np.nanmedian(plx[1:])
            plx_m = np.nanmedian(plx)
        plx_l.append(round(plx_m, 3))

        if np.isnan(pmRA).all():
            pmRA_m = np.nan
        elif not np.isnan(pmRA[1]):
            # Prefer CG20 values whenever possible
            pmRA_m = pmRA[1]
        else:
            # # Avoid HAO21 pmRA data if possible
            # if np.isnan(pmRA[1:]).all():
            #     pmRA_m = pmRA[0]
            # else:
            #     pmRA_m = np.nanmedian(pmRA[1:])
            pmRA_m = np.nanmedian(pmRA)
        pmRA_l.append(round(pmRA_m, 3))

        if np.isnan(pmDE).all():
            pmDE_m = np.nan
        elif not np.isnan(pmDE[1]):
            # Prefer CG20 values whenever possible
            pmDE_m = pmDE[1]
        else:
            # # Avoid HAO21 pmDE data if possible
            # if np.isnan(pmDE[1:]).all():
            #     pmDE_m = pmDE[0]
            # else:
            #     pmDE_m = np.nanmedian(pmDE[1:])
            pmDE_m = np.nanmedian(pmDE)
        pmDE_l.append(round(pmDE_m, 3))

        if np.isnan(Dist).all():
            dmod_l.append(np.nan)
        else:
            # dmod_l.append(round(-5 + 5 * np.log10(np.nanmedian(Dist)), 2))
            dm = []
            for d in Dist:
                if np.isnan(d):
                    dm.append(np.nan)
                else:
                    dm.append(round(-5 + 5 * np.log10(d), 2))
            dmod_l.append("_".join([str(_) for _ in dm]))

        if np.isnan(extin).all():
            extin_l.append(np.nan)
        else:
            extin_l.append("_".join([str(_) for _ in extin]))

        if np.isnan(logAge).all():
            logAge_l.append(np.nan)
        else:
            # logAge_l.append(round(np.nanmedian(logAge), 2))
            la = []
            for a in logAge:
                if np.isnan(a):
                    la.append(np.nan)
                else:
                    la.append(round(a, 2))
            logAge_l.append("_".join([str(_) for _ in la]))

        if np.isnan(FeH).all():
            FeH_l.append(np.nan)
        else:
            # FeH_l.append(round(np.nanmedian(FeH), 2))
            fe = []
            for f in FeH:
                if np.isnan(f):
                    fe.append(np.nan)
                else:
                    fe.append(round(f, 2))
            FeH_l.append("_".join([str(_) for _ in fe]))

    # Store combined databases
    final_DB = {
        'DB': db_l, 'ID': names_l, 'RA_ICRS': ra_l, 'DE_ICRS': dec_l,
        'GLON': glon_l, 'GLAT': glat_l, 'plx': plx_l, 'pmRA': pmRA_l,
        'pmDE': pmDE_l, 'distmod': dmod_l, 'extin': extin_l,
        'logage': logAge_l, 'FeH': FeH_l
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

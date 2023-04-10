
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt

"""
Script used to combine the four 'old' databases (CG20, DIAS21, HAO21, BICA19)
using the clusters' names (ie: not their coordinates).

DIAS21 has bad RA values for some clusters, discarded.
HAO21 has bad pmRA, pmDE, plx values compared to CG20 for some clusters,
discarded whenever possible
"""

dbs_folder = '../0_data/databases/old/'

# Data stored per cluster (per database):
# ID, RA, DEC, plx, pmRA, pmDE,
# Dist, AV, Age, [Fe/H]
N_vals = 9

dbs_names = {
    'HAO21': ['Cluster', 'RA_ICRS', 'DE_ICRS', 'plx', 'pmRA*', 'pmDE',
              None, None, 'logt', None],
    'DIAS21': ['Cluster', 'RA_ICRS', 'DE_ICRS', 'Plx', 'pmRA', 'pmDE',
               'Dist', 'Av', 'logage', '[Fe/H]'],
    'CG20': ["Name", "RA", "DEC", "plx", "pmRA", "pmDE", 'DistPc', 'AVNN',
             'AgeNN', None],
    'BICA19_all_names': [
        'Name', 'GLON', 'GLAT', None, None, None, None, None, None, None]}

# HAO21, DIAS21, CG20, BICA19
DBs_IDS = np.array([1, 2, 3, 4])


def main():
    """
    """
    # simbad = pd.read_csv(
    #     "../0_data/SIMBAD_opc.dat", sep='|', index_col=False, comment='#')
    # simbad.columns = simbad.columns.str.replace(' ', '')
    # for cl in simbad['identifier']:
    #     if not any(x in cl for x in discard):
    #         if ']' in cl:
    #             cl0, cl1 = cl.split(']')
    #             print(cl0, cl1)
    #         else:
    #             print(cl)

    def rem_UBC(data, name_col):
        msk = []
        for cl in data[name_col]:
            if 'UBC' in cl:
                msk.append(False)
            else:
                msk.append(True)
        # Drop rows, update index
        data = data[np.array(msk)]
        data = data.reset_index(drop=True)
        return data

    DB_data, DB_names_match = {}, {}
    for DB, cols in dbs_names.items():
        data = pd.read_csv(
            dbs_folder + DB + ".dat", sep=',', comment='#', index_col=False)

        print(DB, len(data))
        # Remove CASTRO20 duplicates from DIAS21, HAO21, and CG20
        if DB in ('HAO21', 'DIAS21'):
            data = rem_UBC(data, 'Cluster')
        if DB == 'CG20':
            data = rem_UBC(data, 'Name')
        print(DB, len(data))

        DB_data[DB] = data
        names = data[dbs_names[DB][0]]
        names = [_.lower().replace('_', '').replace(' ', '').replace('-', '')
                 for _ in names]
        DB_names_match[DB] = np.array(names)

    cl_dict = HAO21_DIAS21(
        dbs_names, DB_names_match, DB_data, 'HAO21', 'DIAS21')
    print(f"Clusters in cross-matched catalog so far: {len(cl_dict)}")

    cl_dict = HAO21DIAS21_CG20(cl_dict, dbs_names, DB_names_match, DB_data, 'CG20')
    print(f"Clusters in cross-matched catalog so far: {len(cl_dict)}")

    cl_dict = HAO21DIAS21CG20_BICA19(cl_dict, dbs_names, DB_data, 'BICA19_all_names')
    print(f"Clusters in final cross-matched catalog: {len(cl_dict)}")

    final_DB = combine_DBs(cl_dict)

    pd.DataFrame(final_DB).to_csv('OLD_DBs.csv', na_rep='nan', index=False)
    breakpoint()


def combine_DBs(cl_dict):
    """
    Store unique values for each cluster
    """

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
        if msk.sum() == 1 and bool(msk[1]) is True:
            ra_m = ra[1]
        else:
            # Don't use DIAS21's RA
            ra_m = np.nanmedian([ra[0], ra[2], ra_b])
        dec_m = np.nanmedian([dec[0], dec[1], dec[2], dec_b])
        ra_l.append(round(ra_m, 5))
        dec_l.append(round(dec_m, 5))

        lon, lat = radec2lonlat(ra_m, dec_m)
        glon_l.append(lon)
        glat_l.append(lat)

        # Avoid HAO21 Plx data if possible
        if np.isnan(plx).all():
            plx_m = np.nan
        else:
            if np.isnan(plx[1:]).all():
                plx_m = plx[0]
            else:
                plx_m = np.nanmedian(plx[1:])
        plx_l.append(round(plx_m, 3))

        # Avoid HAO21 pmRA data if possible
        if np.isnan(pmRA).all():
            pmRA_m = np.nan
        else:
            if np.isnan(pmRA[1:]).all():
                pmRA_m = pmRA[0]
            else:
                pmRA_m = np.nanmedian(pmRA[1:])
        pmRA_l.append(round(pmRA_m, 3))

        # Avoid HAO21 pmDE data if possible
        if np.isnan(pmDE).all():
            pmDE_m = np.nan
        else:
            if np.isnan(pmDE[1:]).all():
                pmDE_m = pmDE[0]
            else:
                pmDE_m = np.nanmedian(pmDE[1:])
        pmDE_l.append(round(pmDE_m, 3))

        if np.isnan(Dist).all():
            dmod_l.append(np.nan)
        else:
            dmod_l.append(round(-5 + 5 * np.log10(np.nanmedian(Dist)), 2))

        if np.isnan(extin).all():
            extin_l.append(np.nan)
        else:
            extin_l.append(round(np.nanmedian(extin), 2))

        if np.isnan(logAge).all():
            logAge_l.append(np.nan)
        else:
            logAge_l.append(round(np.nanmedian(logAge), 2))

        if np.isnan(FeH).all():
            FeH_l.append(np.nan)
        else:
            FeH_l.append(round(np.nanmedian(FeH), 2))

    # Store combined databases
    final_DB = {
        'DB': db_l, 'ID': names_l, 'RA_ICRS': ra_l, 'DE_ICRS': dec_l,
        'GLON': glon_l, 'GLAT': glat_l, 'plx': plx_l, 'pmRA': pmRA_l,
        'pmDE': pmDE_l, 'distmod': dmod_l, 'extin': extin_l,
        'logage': logAge_l, 'FeH': FeH_l
    }

    return final_DB

    # Check missing data
    values = []
    i_r, i_p = 0, 0
    for cl, vals in cl_dict.items():

        # # Check rad data (no HAO21)
        # dias, cg, bica = vals[5][1:4]
        # if all([np.isnan(bica), np.isnan(dias), np.isnan(cg)]):
        #     i_r += 1

        # # Check nans in Plx that only exist in HAO21
        # hao_p, dias_p, cg_p = vals[2][:3]
        # if all(np.isnan((dias_p, cg_p))) and not np.isnan(hao_p):
        #     print(i_p, cl)
        #     i_p += 1

        # Check nans in pmRA & pmDE in CG20 & DIAS21
        hao_ra, dias_ra, cg_ra = vals[3][:3]
        hao_de, dias_de, cg_de = vals[4][:3]
        # A
        if all(np.isnan((dias_ra, hao_ra, dias_de, hao_de))) and\
                all(~np.isnan((cg_ra, cg_de))):
        # # B
        # if all(np.isnan((cg_ra, hao_ra, cg_de, hao_de))) and\
        #         all(~np.isnan((dias_ra, dias_de))):
        # # C
        # if all(np.isnan((cg_ra, dias_ra, cg_de, dias_de))) and\
        #         all(~np.isnan((hao_ra, hao_de))):
        # # E
        # if all(~np.isnan((hao_ra, dias_ra, hao_de, dias_de))) and\
        #         all(np.isnan((cg_ra, cg_de))):
        # # F
        # if all(~np.isnan((hao_ra, cg_ra, hao_de, cg_de))) and\
        #         all(np.isnan((dias_ra, dias_de))):
        # # D
        # if all(~np.isnan((dias_ra, cg_ra, dias_de, cg_de))) and\
        #         all(np.isnan((hao_ra, hao_de))):
        # G
        # if all(~np.isnan((hao_ra, dias_ra, cg_ra, hao_de, dias_de, cg_de))):
            print(i_p, cl)
            i_p += 1

        # if all(np.isnan((cg_ra, cg_de))):
        # if all(np.isnan((dias_ra, dias_de))):
        # if all(np.isnan((hao_ra, hao_de))):
        #     print(i_p, cl)
        #     i_p += 1

        # # Check nans in Plx & pmRA & pmDE (no Bica)
        # hao_p, dias_p, cg_p = vals[2][:3]
        # hao_ra, dias_ra, cg_ra = vals[3][:3]
        # hao_de, dias_de, cg_de = vals[4][:3]
        # if all(np.isnan((hao_p, dias_p, cg_p, hao_ra, dias_ra, cg_ra, hao_de,
        #                  dias_de, cg_de))):
        #     print(i_p, cl)
        #     i_p += 1

    #     # Store distances to plot
    #     hao_p, dias_p, cg_p = vals[2][:3]
    #     values.append([hao_p, dias_p, cg_p])

    # hao, dias, cg = np.array(values).T
    # # xymin, xymax = 5.5, 10.5
    # xymin, xymax = 0, 3

    # plt.subplot(131)
    # plt.scatter(hao, dias)
    # plt.plot((xymin, xymax), (xymin, xymax), c='r')
    # plt.xlim(xymin, xymax)
    # plt.ylim(xymin, xymax)
    # plt.xlabel("HAO21")
    # plt.ylabel("DIAS21")

    # plt.subplot(132)
    # plt.scatter(hao, cg)
    # plt.plot((xymin, xymax), (xymin, xymax), c='r')
    # plt.xlim(xymin, xymax)
    # plt.ylim(xymin, xymax)
    # plt.xlabel("HAO21")
    # plt.ylabel("CG20")

    # plt.subplot(133)
    # plt.scatter(dias, cg)
    # plt.plot((xymin, xymax), (xymin, xymax), c='r')
    # plt.xlim(xymin, xymax)
    # plt.ylim(xymin, xymax)
    # plt.xlabel("DIAS21")
    # plt.ylabel("CG20")
    # plt.show()


def HAO21_DIAS21(dbs_names, DB_names_match, DB_data, DB_id1, DB_id2):
    """
    Find HAO21 clusters in DIAS21
    """
    cl_dict = {}

    cols1, cols2 = dbs_names[DB_id1], dbs_names[DB_id2]
    dias21_idxs_match = []
    for i, cl in enumerate(DB_names_match[DB_id1]):
        # HAO21
        full_name = DB_data[DB_id1][cols1[0]][i]
        # Save cluster data
        cl_dict[full_name] = [[] for _ in range(N_vals)]

        for q, col in enumerate(cols1[1:]):
            if col is not None:
                cl_dict[full_name][q].append(DB_data[DB_id1][col][i])
            else:
                cl_dict[full_name][q].append(np.nan)

        # DIAS21
        j = np.where(DB_names_match[DB_id2] == cl)[0]
        if len(j) > 0:
            # Match found
            j = j[0]
            dias21_idxs_match.append(j)
            for q, col in enumerate(cols2[1:]):
                if col is not None:
                    cl_dict[full_name][q].append(DB_data[DB_id2][col][j])
                else:
                    cl_dict[full_name][q].append(np.nan)
        else:
            # No match for this cluster in DIAS21
            for q, col in enumerate(cols2[1:]):
                cl_dict[full_name][q].append(np.nan)
    print("\nClusters in HAO21 matched to clusters in DIAS21: "
          + f"{len(dias21_idxs_match)}")
    print("Clusters in HAO21 not matched to clusters in DIAS21: "
          + f"{len(DB_data[DB_id1]) - len(dias21_idxs_match)}")

    # Store clusters from DIAS21 with no match in HAO21
    DIAS21_all_idxs = np.arange(0, len(DB_data[DB_id2]))
    DIAS21_no_match = list(set(DIAS21_all_idxs) - set(dias21_idxs_match))
    print("Clusters in DIAS21 not matched to clusters in HAO21: "
          + f"{len(DIAS21_no_match)}")
    for i in DIAS21_no_match:
        full_name = DB_data[DB_id2][cols2[0]][i]
        # Save cluster data
        cl_dict[full_name] = [[np.nan] for _ in range(N_vals)]
        for q, col in enumerate(cols2[1:]):
            if col is not None:
                cl_dict[full_name][q].append(DB_data[DB_id2][col][i])
            else:
                cl_dict[full_name][q].append(np.nan)

    return cl_dict


def HAO21DIAS21_CG20(cl_dict, dbs_names, DB_names_match, DB_data, DB_id3):
    """
    Match HAO21+DIAS21 with CG20
    """
    N_old = len(cl_dict)
    cols3 = dbs_names[DB_id3]
    cg20_idxs_match = []
    for i, (cl, vals) in enumerate(cl_dict.items()):
        cl_match = cl.lower().replace('_', '').replace(' ', '').replace('-', '')
        # CG20
        j = np.where(DB_names_match[DB_id3] == cl_match)[0]

        if len(j) > 0:
            # Match found
            j = j[0]
            cg20_idxs_match.append(j)
            for q, col in enumerate(cols3[1:]):
                if col is not None:
                    cl_dict[cl][q].append(DB_data[DB_id3][col][j])
                else:
                    cl_dict[cl][q].append(np.nan)
        else:
            # No match for this cluster in CG20
            for q, col in enumerate(cols3[1:]):
                cl_dict[cl][q].append(np.nan)
    print("\nClusters in HAO21+DIAS21 matched to clusters in CG20: "
          + f"{len(cg20_idxs_match)}")
    print("Clusters in HAO21+DIAS21 not matched to clusters in CG20: "
          + f"{N_old - len(cg20_idxs_match)}")

    # Store clusters from CG20 with no match in HAO21+DIAS21
    CG20_all_idxs = np.arange(0, len(DB_data[DB_id3]))
    CG20_no_match = list(set(CG20_all_idxs) - set(cg20_idxs_match))
    print("Clusters in CG20 not matched to clusters in HAO21+DIAS21: "
          + f"{len(CG20_no_match)}")
    for i in CG20_no_match:
        full_name = DB_data[DB_id3][cols3[0]][i]
        # Save cluster data
        cl_dict[full_name] = [[np.nan, np.nan] for _ in range(N_vals)]
        for q, col in enumerate(cols3[1:]):
            if col is not None:
                cl_dict[full_name][q].append(DB_data[DB_id3][col][i])
            else:
                # cl_dict[full_name][q].append([np.nan, np.nan])
                cl_dict[full_name][q].append(np.nan)

    return cl_dict


def HAO21DIAS21CG20_BICA19(cl_dict, dbs_names, DB_data, DB_id4):
    """
    Match HAO21+DIAS21+CG20 with BICA19
    """
    N_old = len(cl_dict)
    cols4 = dbs_names[DB_id4]

    DB_names_match = []
    for names in DB_data[DB_id4][cols4[0]]:
        all_names = names.split('|')
        nm = [_.lower().replace('_', '').replace(' ', '').replace('-', '')
              for _ in all_names]
        DB_names_match.append(nm)

    bica19_idxs_match = []
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
            bica19_idxs_match.append(j)
            for q, col in enumerate(cols4[1:]):
                if col is not None:
                    val = DB_data[DB_id4][col][j]
                    # # Transform diameter to r50 in degrees
                    # if q == 5:
                    #     val = round(.5 * (val / 60.), 3)
                    cl_dict[cl][q].append(val)
                else:
                    cl_dict[cl][q].append(np.nan)
        else:
            # No match for this cluster in BICA19
            for q, col in enumerate(cols4[1:]):
                cl_dict[cl][q].append(np.nan)
    print("\nClusters in HAO21+DIAS21+CG20 matched to clusters in BICA19: "
          + f"{len(bica19_idxs_match)}")
    print("Clusters in HAO21+DIAS21+CG20 not matched to clusters in BICA19: "
          + f"{N_old - len(bica19_idxs_match)}")

    # Store clusters from BICA19 with no match in HAO21+DIAS21+CG20
    BICA19_all_idxs = np.arange(0, len(DB_data[DB_id4]))
    BICA19_no_match = list(set(BICA19_all_idxs) - set(bica19_idxs_match))
    print("Clusters in BICA19 not matched to clusters in HAO21+DIAS21+CG20: "
          + f"{len(BICA19_no_match)}")
    for i in BICA19_no_match:
        full_name = DB_data[DB_id4][cols4[0]][i].split('|')[0]
        # Save cluster data
        cl_dict[full_name] = [[np.nan, np.nan, np.nan] for _ in range(N_vals)]
        for q, col in enumerate(cols4[1:]):
            if col is not None:
                val = DB_data[DB_id4][col][i]
                # # Transform diameter to r50 in degrees
                # if q == 5:
                #     val = round(.5 * (val / 60.), 3)
                cl_dict[full_name][q].append(val)
            else:
                # cl_dict[full_name][q].append([np.nan, np.nan, np.nan])
                cl_dict[full_name][q].append(np.nan)

    return cl_dict


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

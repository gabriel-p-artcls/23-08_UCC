
import numpy as np
import pandas as pd
import csv
from astropy.coordinates import angular_separation
from scipy.spatial.distance import cdist

path_io = "../2_pipeline/"
db_in = 'almost_final_DB.csv'
db_dup = 'duplicates.csv'
db_out = 'final_DB.csv'


def main(N_dups=10):
    """
    """
    # Read file data will all the clusters
    df = pd.read_csv(path_io + db_in)

    idxs_dup = dups_remove(df, N_dups)

    # Store removed duplicates
    df_dup = df.iloc[idxs_dup]
    df_dup.to_csv(path_io + db_dup, na_rep='nan', index=False,
                  quoting=csv.QUOTE_NONNUMERIC)

    # Store final catalogue
    df.drop(df.index[idxs_dup], inplace=True)
    df = df.reset_index(drop=True)
    df.to_csv(path_io + db_out, na_rep='nan', index=False,
              quoting=csv.QUOTE_NONNUMERIC)
    print("\nFinal database written to file")


def dups_remove(df, N_dups):
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
    clusts = []
    for i, cl in enumerate(dist):
        idx = np.argsort(cl)[:N_dups]

        clusts_matched = []
        for j in idx:
            # Angular distance in arcmin (rounded)
            d = round(angular_separation(x[i], y[i], x[j], y[j]) * 60, 2)
            # PMs distance
            pm_d = np.sqrt((pmRA[i]-pmRA[j])**2 + (pmDE[i]-pmDE[j])**2)
            # Parallax distance
            plx_d = abs(plx[i] - plx[j])

            dup_flag = duplicate_find(d, pm_d, plx_d, plx[i])

            if dup_flag:
                clusts_matched.append(
                    [df['DB'][j], df['ID'][j], d, round(pm_d, 2),
                     round(plx_d, 3)])

        if clusts_matched:
            clusts.append([
                i, df['DB'][i], df['ID'][i], clusts_matched])

    idxs_dup, count = [], 1
    for i in range(1, 32):
        print(f"\nID: {i}")
        idxs_dup, count = extract_DB_matches(clusts, str(i), idxs_dup, count)

    return idxs_dup


def extract_DB_matches(clusts, DB_ID, idxs_dup, count):
    j = 0
    for cl in clusts:
        i, DB, ID, clusts_matched = cl
        # Process the catalogue identified by DB_ID
        if DB != DB_ID:
            continue

        # Only search for duplicates in catalogues published *before* this one
        cls_matched = []
        for cl_d in clusts_matched:
            DB_j = cl_d[0].split('_')
            for DB_jj in DB_j:
                if int(DB_jj) < int(DB_ID):
                    cls_matched.append(cl_d)
                    break

        if cls_matched:
            print("{} ({}) {} {}".format(j, count, i + 2, ID))
            idxs_dup.append(i)
            count += 1
            j += 1
            for cl_d in cls_matched:
                print("  ", cl_d)

    return idxs_dup, count


def duplicate_find(d, pm_d, plx_d, plx):
    """
    """
    # dup_prob = 0.
    # if np.isnan(pm_d) or np.isnan(plx_d):
    #     return dup_prob

    if plx >= 4:
        rad, plx_r, pm_r = 15, 0.25, 0.5
    elif 3 <= plx and plx < 4:
        rad, plx_r, pm_r = 10, 0.15, 0.25
    elif 2 <= plx and plx < 3:
        rad, plx_r, pm_r = 5, 0.1, 0.15
    elif 1 <= plx and plx < 2:
        rad, plx_r, pm_r = 2.5, 0.05, 0.1
    else:
        rad, plx_r, pm_r = 1, 0.025, 0.05

    if pm_d < pm_r:
        if plx_d < plx_r:
            if d < rad:
                return True
                # Probability of being a duplicate. Linear trend from 0.5 to 1
                # dup_prob = round(-.5*d/rad + 1, 2)

    return False


if __name__ == '__main__':
    main()

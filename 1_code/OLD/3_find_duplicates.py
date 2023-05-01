
import datetime
import numpy as np
import pandas as pd
import csv
from astropy.coordinates import angular_separation
from scipy.spatial.distance import cdist

path_io = "../2_pipeline/"
db_in = '2_standard_names_DB.csv'
db_out = 'UCC_cat_'


def main(N_dups=10):
    """
    """
    # Read file data will all the clusters
    df = pd.read_csv(path_io + db_in)

    dups_fnames = dups_identify(df, N_dups)

    # Store duplicates in final catalogue
    df['dups_fnames'] = dups_fnames

    # Save to file
    d = datetime.datetime.now()
    date = d.strftime('%Y%m%d')
    df.to_csv(path_io + db_out + date + '.csv', na_rep='nan', index=False,
              quoting=csv.QUOTE_NONNUMERIC)
    print("\nFinal database written to file")


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
    main()


import numpy as np
import pandas as pd
import csv
import astropy.units as u
from astropy.coordinates import SkyCoord


final_dbs_path = "/home/gabriel/Github/UCC/add_New_DB/UCC_cat_20230619.csv"
hunt23_DB_path = "/home/gabriel/Github/UCC/add_New_DB/databases/HUNT23.csv"
hunt23_membs_path = "../0_data/hunt23_members.parquet"


def main(rad_max=1):
    """
    Check the difference between the centers stored in HUNT23 final database
    and those obtained as the median of its members.
    """
    # Read HUNT23 members data
    hunt23_membs = pd.read_parquet(hunt23_membs_path)
    hunt23_names = hunt23_membs['Name']

    final_db = pd.read_csv(final_dbs_path)
    hunt23_DB = pd.read_csv(hunt23_DB_path)

    ours_hist, hunt23_hist, match_hist = [np.zeros(14) for _ in range(3)]
    for i, row in final_db.iterrows():

        if 'HUNT23' not in row['DB']:
            continue

        # Identify cluster in HUNT23 DB
        dbs_idx = row['DB_i'].split(';')
        h23_j = row['DB'].split(';').index('HUNT23')
        hunt23_idx = int(dbs_idx[h23_j])
        hunt23_name = hunt23_DB['Name'][hunt23_idx].strip()

        if 'HSC' not in hunt23_name:
            continue

        # Identify the members for this cluster in HUNT23
        if hunt23_name.startswith('VDBH_'):
            hunt23_name = 'BH_' + hunt23_name.split(',')[0].split('_')[1]
        if hunt23_name.startswith('VDB_'):
            hunt23_name = 'vdBergh_' + hunt23_name.split(',')[0].split('_')[1]
        msk = hunt23_names == hunt23_name
        if msk.sum() == 0:
            print(f"Could not find cluster {hunt23_name} in HUNT23")
            continue

        l1, b1 = hunt23_DB['GLON'][hunt23_idx], hunt23_DB['GLAT'][hunt23_idx]
        l2, b2 = np.median(hunt23_membs['GLON'][msk]),\
            np.median(hunt23_membs['GLAT'][msk])
        dist = np.sqrt((l1-l2)**2+(b1-b2)**2)
        if dist > rad_max:
            print("{}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}".format(
                hunt23_name, l1, b1, l2, b2, dist))

            # Re-write values
            hunt23_DB.at[hunt23_idx, 'GLON'] = l2
            hunt23_DB.at[hunt23_idx, 'GLAT'] = b2
            # Now equatorial
            gc = SkyCoord(l=l2*u.degree, b=b2*u.degree, frame='galactic')
            hunt23_DB.at[hunt23_idx, 'RA_ICRS'] = gc.fk5.ra.value
            hunt23_DB.at[hunt23_idx, 'DE_ICRS'] = gc.fk5.dec.value
    # Update HUNT23 database
    hunt23_DB.to_csv(hunt23_DB_path, na_rep='nan', index=False,
                     quoting=csv.QUOTE_NONNUMERIC)


if __name__ == '__main__':
    main()

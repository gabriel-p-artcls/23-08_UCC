
import numpy as np
from astropy import units as u
from astropy.coordinates import angular_separation
import pandas as pd
UCC_cat = "/home/gabriel/Github/UCC/add_New_DB/UCC_cat_20230619.csv"
GCs_cat = "/home/gabriel/Github/UCC/add_New_DB/databases/globulars.csv"


def main(rad=10):
    """
    """
    df_UCC = pd.read_csv(UCC_cat)
    df_gcs = pd.read_csv(GCs_cat)
    l_gc, b_gc = df_gcs['GLON'].values, df_gcs['GLAT'].values

    for idx, row in df_UCC.iterrows():
        l_ucc, b_ucc = row['GLON'], row['GLAT']

        d_arcmin1 = angular_separation(
            l_ucc*u.deg, b_ucc*u.deg, l_gc*u.deg, b_gc*u.deg).to('deg').value * 60
        j1 = np.argmin(d_arcmin1)
        if d_arcmin1[j1] < rad:
            print(f"{round(d_arcmin1[j1], 2)}, {row['ID']}, {df_gcs['Name'][j1].strip()}, {df_gcs['Nstar'][j1]}")


if __name__ == '__main__':
    # plt.style.use('science')
    main()

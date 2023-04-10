
import numpy as np
import pandas as pd
from astropy.coordinates import angular_separation
from scipy.spatial.distance import cdist


def main():
    """
    Find the closest clusters to all clusters
    """
    # Read file data will all the clusters
    data_all_cls = pd.read_csv("final_DB.csv")

    x, y = data_all_cls['GLON'], data_all_cls['GLAT']
    pmRA, pmDE, plx = data_all_cls['pmRA'], data_all_cls['pmDE'],\
        data_all_cls['plx']
    coords = np.array([x, y]).T

    # Find the distances to all clusters, for all clusters
    dist = cdist(coords, coords)
    # Change distance to itself from 0 to inf
    msk = dist == 0.
    dist[msk] = np.inf

    clusts = []
    for i, cl in enumerate(dist):
        idx_5 = np.argsort(cl)[:5]

        clusts_5 = []
        for j in idx_5:
            # Angular distance in arcmin (rounded)
            d = round(angular_separation(x[i], y[i], x[j], y[j]) * 60, 2)
            # PMs distance
            pm_d = np.sqrt((pmRA[i]-pmRA[j])**2 + (pmDE[i]-pmDE[j])**2)
            # Parallax distance
            plx_d = abs(plx[i] - plx[j])

            flag_append = False
            if d < 5:
                if pm_d < .25 or np.isnan(pm_d):
                    if plx_d < .05 or np.isnan(plx_d):
                        flag_append = True

            if flag_append:
                clusts_5.append(
                    [data_all_cls['DB'][j], data_all_cls['ID'][j], d,
                     round(pm_d, 2), round(plx_d, 3)])

        if clusts_5:
            clusts.append([
                data_all_cls['DB'][i], data_all_cls['ID'][i], clusts_5])

    for i, cl in enumerate(clusts):
        DB, ID, clusts_5 = cl
        print(i, DB, ID)
        for cl_5 in clusts_5:
            print("  ", cl_5)


if __name__ == '__main__':
    main()

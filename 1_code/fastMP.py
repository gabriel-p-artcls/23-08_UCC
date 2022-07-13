
import warnings
from os import listdir
import numpy as np
from astropy.stats import RipleysKEstimator
import pandas as pd
from scipy.spatial.distance import cdist
import time as tm

N_membs = 25
n_clusters = 1000
N_loop = 10
zoom_f = 3 #1.5
N_groups_large_dist = 100
large_dist_perc = 0.75
N_resample = 250

read_PM_center = False
CG_data_file = "../0_data/cantat_gaudin_et_al_2020/cg2020.csv"

out_path = "../2_pipeline/new_method/78_large_fr_flags/"


def main():
    """
    """
    # path = "../0_data/383_largest/303_no_flags/47_larger_5_Mb/"
    # sep = ","

    # path = "../2_pipeline/pyUPMASK/303_no_flags/256_GMM/"
    # path = "../2_pipeline/pyUPMASK/47_larger_5_Mb/33_smaller_10/"
    # path = "../2_pipeline/pyUPMASK/303_no_flags/47_larger_5_Mb/13_larger_10/"
    # path = "../2_pipeline/pyUPMASK/303_no_flags/47_larger_5_Mb/1_larger_50/"
    path = "../0_data/383_largest/78_large_fr_flags/"
    sep = ","

    files = listdir(path)

    for file in files:
        # if file != "Alessi_Teutsch_5.csv":
        #     continue

        print(file)
        start = tm.time()

        data = pd.read_csv(path + file, sep=sep)
        lon, lat = data['GLON'], data['GLAT']

        xmin, xmax = lon.min(), lon.max()
        ymin, ymax = lat.min(), lat.max()
        area = (xmax - xmin) * (ymax - ymin)
        Kest = RipleysKEstimator(
            area=area, x_max=xmax, y_max=ymax, x_min=xmin, y_min=ymin)
        # https://rdrr.io/cran/spatstat/man/Kest.html
        # "For a rectangular window it is prudent to restrict the r values to a
        # maximum of 1/4 of the smaller side length of the rectangle
        # (Ripley, 1977, 1988; Diggle, 1983)"
        lmin = min((xmax - xmin), (ymax - ymin))
        rads = np.linspace(lmin * 0.01, lmin * .25, 10)

        N_stars = data.shape[0]

        data_3 = np.array([data['pmRA'], data['pmDE'], data['Plx']])
        data_err = np.array([data['e_pmRA'], data['e_pmDE'], data['e_Plx']])

        if read_PM_center:
            df = pd.read_csv(CG_data_file)
            idx = df.index[df['Name'] == file.split('.')[0]][0]
            vpd_c = (df['pmRA'][idx], df['pmDE'][idx])
            print("VPD center read")

        idx_selected = []
        for _ in range(N_resample):
            # Gaussian random sample
            grs = np.random.normal(0., 1., N_stars)
            sampled_data = data_3 + grs * data_err
            s_pmRA, s_pmDE, s_Plx = sampled_data

            # Estimate VPD center
            if read_PM_center is False:
                vpd = np.array([s_pmRA, s_pmDE]).T
                vpd_c = centVPD(vpd)

            st_idx = getStars(
                N_stars, Kest, rads, lon, lat, s_pmRA, s_pmDE, s_Plx, vpd_c)

            idx_selected += st_idx

        # Assign probabilities
        values, counts = np.unique(idx_selected, return_counts=True)
        probs = counts / N_resample
        probs_all = np.zeros(N_stars)
        probs_all[values] = probs

        print(len(data), tm.time() - start)

        # Save file
        data['probs_final'] = probs_all
        data.to_csv(out_path + file, index=False)


def centVPD(vpd, N_dim=2):
    """
    Estimate the center of the cluster in the proper motions space.

    This is the most important function in the algorithm If this function
    fails, everything else will fail too.
    """
    for _ in range(N_loop):
        N_stars = vpd.shape[0]

        if N_stars < N_membs * 4:
            break

        N_bins = max(3, int(N_stars / N_membs))
        if N_bins**N_dim > n_clusters:
            N_bins = int(n_clusters**(1 / N_dim))

        H, edges = np.histogramdd(vpd, bins=(N_bins, N_bins))
        edgx, edgy = edges

        # Find center coordinates
        flat_idx = H.argmax()
        cbx, cby = np.unravel_index(flat_idx, H.shape)
        cx = (edgx[cbx + 1] + edgx[cbx]) / 2.
        cy = (edgy[cby + 1] + edgy[cby]) / 2.

        # This makes the process very non-robust
        # Hg = gaussian_filter(H, sigma=5)
        # flat_idx = Hg.argmax()
        # cbx, cby = np.unravel_index(flat_idx, Hg.shape)
        # cx = (edgx[cbx + 1] + edgx[cbx]) / 2.
        # cy = (edgy[cby + 1] + edgy[cby]) / 2.

        # Zoom in
        x, y = vpd.T
        rx, ry = edgx[1] - edgx[0], edgy[1] - edgy[0]
        msk = (x < (cx + zoom_f * rx)) & (x > (cx - zoom_f * rx))\
            & (y < (cy + zoom_f * ry)) & (y > (cy - zoom_f * ry))
        vpd = vpd[msk]

    return (cx, cy)


def getStars(N_stars, Kest, rads, lon, lat, pmRA, pmDE, Plx, vpd_c):
    """
    """
    # Distance to VPD center
    dist = (pmRA - vpd_c[0])**2 + (pmDE - vpd_c[1])**2

    # Closest stars to center
    idx = dist.argsort()[: 2 * N_membs]  # HARDCODED

    # Center in Plx
    plx_c = np.median(Plx[idx])
    cent_all = np.array([vpd_c[0], vpd_c[1], plx_c])
    cent_all = cent_all.reshape(1, 3)
    data_all = np.array([pmRA, pmDE, Plx]).T

    # 3D distance to the estimated center
    dist = cdist(data_all, cent_all).T[0]
    # Sort by smallest
    d_idxs = dist.argsort()

    # Estimate average C_s of large distance stars
    C_S_field = []
    N_low1 = N_stars - N_groups_large_dist * N_membs
    N_low2 = int(N_stars * large_dist_perc)
    N_low = max(N_low1, N_low2)
    step_old = -1
    for step in np.arange(N_stars - N_membs, N_low, -N_membs):
        msk = d_idxs[step:step_old]
        xy = np.array([lon[msk], lat[msk]]).T
        C_s = rkfunc(xy, rads, Kest)
        if not np.isnan(C_s):
            C_S_field.append(C_s)
        step_old = step
    C_thresh = np.median(C_S_field) + np.std(C_S_field)  # HARDCODED

    # Find Ripley's K value where the stars differ from the estimated
    # field value
    step_old = 0
    for step in np.arange(N_membs, N_stars, N_membs):
        msk = d_idxs[step_old:step]
        xy = np.array([lon[msk], lat[msk]]).T
        C_s = rkfunc(xy, rads, Kest)
        if not np.isnan(C_s):
            if C_s < C_thresh:
                break
        step_old = step

    return list(d_idxs[:step])


def rkfunc(xy, rads, Kest):
    """
    Test how similar this cluster's (x, y) distribution is compared
    to a uniform random distribution using Ripley's K.
    https://stats.stackexchange.com/a/122816/10416
    """
    # Avoid large memory consumption if the data array is too big
    if xy.shape[0] > 5000:
        mode = "none"
    else:
        mode = 'translation'

    # Hide RunTimeWarning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        L_t = Kest.Lfunction(xy, rads, mode=mode)

    # Catch all-nans
    if np.isnan(L_t).all():
        C_s = np.nan
    else:
        C_s = np.nanmax(abs(L_t - rads))

    return C_s


if __name__ == '__main__':
    main()

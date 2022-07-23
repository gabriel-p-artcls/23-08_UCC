
import warnings
from os import listdir
import numpy as np
from astropy.stats import RipleysKEstimator
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import multivariate_normal
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import MinMaxScaler
# import matplotlib.pyplot as plt
# import time as tm


"""Summary

Attributes
----------
CG_data_file : str
    Description
large_dist_perc : float
    Description
n_clusters : int
    Description
N_groups_large_dist : int
    Description
N_loop : int
    Description
N_membs : int
    Description
N_resample : int
    Description
out_path : str
    Description
read_PM_center : bool
    Description
zoom_f : int
    Description
"""

N_membs = 25
n_clusters = 1000
N_loop = 10
zoom_f = 3
N_groups_large_dist = 100
large_dist_perc = 0.75
N_resample = 100

read_PM_center = False
CG_data_file = "../0_data/cantat_gaudin_et_al_2020/cg2020.csv"

in_path = "../0_data/383_largest/78_large_fr_flags/"
# in_path = "../2_pipeline/new_method/303_no_flags/256_GMM/"
sep = ","


# import sklearn.cluster as skclust
# model = skclust.MiniBatchKMeans()
# model = skclust.OPTICS()
# model.n_clusters = 1
# from sklearn.covariance import EllipticEnvelope
# from sklearn.neighbors import LocalOutlierFactor
# from sklearn.ensemble import IsolationForest


def main():
    """
    """

    from pathlib import Path
    out_path = "../2_pipeline/new_method/78_large_fr_flags/"
    # out_path = "../2_pipeline/new_method/303_no_flags/256_GMM/M51/"

    i = 58
    # N_loop, N_std = 2, 3
    for N_loop in (2, 3, 4, 5):
        for N_std in (2, 3):
            print(f"M{i}")

            out_path_f = out_path + "M" + str(i) + "/"
            Path(out_path_f).mkdir(parents=True, exist_ok=True)

            run(out_path_f, N_loop, N_std)
            i += 1


def run(out_path_f, N_loop, N_std=3, Rk_std=0):
    """
    """
    # files = listdir(in_path)
    # np.random.shuffle(files)
    # files = files[:15]
    files = (
        'IC_2602.csv', 'UPK_545.csv', 'Collinder_132.csv', 'RSG_8.csv',
        'UPK_495.csv', 'UPK_585.csv', 'LP_5.csv', 'UBC_8.csv', 'UPK_442.csv',
        'Gulliver_9.csv', 'UPK_552.csv', 'Alessi_2.csv', 'UPK_230.csv',
        'Stock_2.csv', 'COIN-Gaia_13.csv', 'ASCC_41.csv', 'UBC_31.csv',
        'NGC_1039.csv', 'ASCC_127.csv', 'NGC_2451A.csv')

    for file in files:
        # if file != "LP_5.csv":
        #     continue
        print(file)

        data = pd.read_csv(in_path + file, sep=sep)
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

        idx_selected = []
        for _ in range(N_resample):
            # Gaussian random sample
            grs = np.random.normal(0., 1., N_stars)
            sampled_data = data_3 + grs * data_err
            s_pmRA, s_pmDE, s_Plx = sampled_data
            data_all = np.array([s_pmRA, s_pmDE, s_Plx]).T

            # Estimate VPD center
            if read_PM_center is False:
                vpd_c = centVPD(np.array([s_pmRA, s_pmDE]).T)
            plx_c = centPlx(s_pmRA, s_pmDE, s_Plx, vpd_c)
            cent_all = np.array([vpd_c[0], vpd_c[1], plx_c])

            # # Normalize
            # if norm_f:
            #     cent_all, data_all = normalize(cent_all, data_all)

            st_idx = getStars(
                Rk_std, N_stars, Kest, rads, lon, lat, cent_all,
                data_all, data_err)

            st_idx = sigmaClip(lon, lat, s_pmRA, s_pmDE, s_Plx, st_idx, N_loop, N_std)

            # # # GUMM cleaning
            # if gumm_f:
            #     clust_xy = np.array([lon[st_idx].values, lat[st_idx].values]).T
            #     clust_xy = MinMaxScaler().fit(clust_xy).transform(clust_xy)
            #     gumm_p = GUMMProbs(clust_xy)
            #     prob_cut = GUMMProbCut(gumm_p)
            #     msk = gumm_p > prob_cut
            #     st_idx = list(np.array(st_idx)[msk])

            # probs = np.zeros(N_stars)
            # probs[st_idx] = gumm_p
            # idx_selected.append(probs)

            idx_selected += st_idx

        # Assign probabilities
        values, counts = np.unique(idx_selected, return_counts=True)
        probs = counts / N_resample
        probs_all = np.zeros(N_stars)
        probs_all[values] = probs

        # probs_all = np.mean(idx_selected, 0)

        # Save file
        data['probs_final'] = probs_all
        data.to_csv(out_path_f + file, index=False)


def centVPD(vpd, N_dim=2):
    """
    Estimate the center of the cluster in the proper motions space.
    
    This is the most important function in the algorithm If this function
    fails, everything else will fail too.
    
    Parameters
    ----------
    vpd : TYPE
        Description
    N_dim : int, optional
        Description
    
    Returns
    -------
    TYPE
        Description
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


def centPlx(pmRA, pmDE, Plx, vpd_c, M=2):
    """
    Estimate the center of the cluster in parallax.
    
    M: multiplier factor that determines how many N_membs are used to estimate
    the parallax center.
    
    Parameters
    ----------
    pmRA : TYPE
        Description
    pmDE : TYPE
        Description
    Plx : TYPE
        Description
    vpd_c : TYPE
        Description
    M : int, optional
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    # Distance to VPD center
    dist = (pmRA - vpd_c[0])**2 + (pmDE - vpd_c[1])**2
    # Closest stars to center
    idx = dist.argsort()[: M * N_membs]  # HARDCODED
    # Center in Plx
    plx_c = np.median(Plx[idx])

    return plx_c


def normalize(cent_all, data_all):
    """
    Parameters
    ----------
    X : TYPE
        Description
    
    Returns
    -------
    TYPE
        Description
    
    """
    # Add center to data
    X = np.append(data_all, [cent_all], 0)
    # Normalize
    data_all = RobustScaler().fit(X).transform(X)
    # Remove center from data an re-write
    cent_all = data_all[-1]
    data_all = data_all[:-1]

    return cent_all, data_all


def getStars(Rk_std, N_stars, Kest, rads, lon, lat, cent, data, data_err):
    """
    Parameters
    ----------
    N_stars : TYPE
        Description
    Kest : TYPE
        Description
    rads : TYPE
        Description
    lon : TYPE
        Description
    lat : TYPE
        Description
    cent : TYPE
        Description
    data : TYPE
        Description
    
    Returns
    -------
    TYPE
        Description
    """
    # 3D distance to the estimated (VPD + Plx) center
    cent = cent.reshape(1, 3)

    # distance of all the stars to the center
    # if e_weight_f:
    #     dataT = data.T
    #     dist = np.sqrt((
    #         dataT[0] - cent[0][0])**2 * data_err[0]
    #         + (dataT[1] - cent[0][1])**2 * data_err[1]
    #         + (dataT[2] - cent[0][2])**2 * data_err[2])
    # else:
    dist = cdist(data, cent).T[0]

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

    # This value is associated to the field stars distribution. When it
    # is achieved, the block below breaks out. The larger this value, the more
    # restrictive the process becomes.
    C_thresh = np.median(C_S_field) + Rk_std * np.std(C_S_field)

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

    return d_idxs[:step]


def rkfunc(xy, rads, Kest):
    """
    Test how similar this cluster's (x, y) distribution is compared
    to a uniform random distribution using Ripley's K.
    https://stats.stackexchange.com/a/122816/10416
    
    Parameters
    ----------
    xy : TYPE
        Description
    rads : TYPE
        Description
    Kest : TYPE
        Description
    
    Returns
    -------
    TYPE
        Description
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


def sigmaClip(lon, lat, s_pmRA, s_pmDE, s_Plx, st_idx, N_loop, N_std):
    """
    Remove outliers in the coordinates space.
    """
    # xy = np.array([lon[st_idx], lat[st_idx]]).T

    xyz = np.array([s_pmRA[st_idx], s_pmDE[st_idx], s_Plx[st_idx]]).T

    # msk = IsolationForest().fit_predict(xy) > 0
    # return list(st_idx[msk])
    # msk = LocalOutlierFactor().fit_predict(xy) > 0
    # return list(st_idx[msk])

    # idxs = np.array(st_idx)
    # for _ in range(N_loop):
    #     xy_d = cdist(xy, xy.mean(0).reshape(-1, 2)).T[0]
    #     xy_std = xy.std(0).mean()
    #     msk_s3 = xy_d < N_std * xy_std
    #     xy, idxs = xy[msk_s3], idxs[msk_s3]

    # return list(idxs)

    idxs = np.array(st_idx)
    for _ in range(N_loop):
        xy_d = cdist(xyz, xyz.mean(0).reshape(-1, 3)).T[0]
        xy_std = xyz.std(0).mean()
        msk_s3 = xy_d < N_std * xy_std
        xyz, idxs = xyz[msk_s3], idxs[msk_s3]

    return list(idxs)

    # xy_d = cdist(xy, xy.mean(0).reshape(-1, 2)).T[0]
    # xy_std = xy.std(0).mean()
    # msk_s3 = xy_d < N_std * xy_std
    # xy, idxs = xy[msk_s3], idxs[msk_s3]
    # print(len(idxs))
    # xy_d = cdist(xy, xy.mean(0).reshape(-1, 2)).T[0]
    # xy_std = xy.std(0).mean()
    # msk_s3 = xy_d < N_std * xy_std
    # idxs = idxs[msk_s3]
    # print(len(idxs))
    # breakpoint()

    # return list(idxs)


# def GUMMProbs(xy, n_epochs=1000, stable_per=.1):
#     """
#     Fit a model composed of a 2D Gaussian and a 2D uniform distribution in a
#     square region with [0., 1.] range.

#     Based on the GMM model implementation shown in:
#     https://towardsdatascience.com/gaussian-mixture-models-explained-6986aaf5a95
#     """
#     # Initialize the 2D Gaussian parameters, and the weights for both
#     # distributions.
#     mu = np.random.uniform(.1, .9, (2,))
#     cov = np.eye(2) * np.random.uniform(.1, .9, (2, 2))
#     pi_u, pi_g = .5, .5
#     N0, N1 = xy.shape

#     lkl_old, nstable = -np.inf, 0
#     for i in range(n_epochs):

#         # expectation_step
#         # Evaluate Gaussian distribution
#         try:
#             gamma_g = pi_g * multivariate_normal(mean=mu, cov=cov).pdf(xy)
#         except np.linalg.LinAlgError:
#             continue
#         # Evaluate uniform distribution (just a constant)
#         gamma_u = pi_u * np.ones(N0)
#         # Normalizing constant
#         gammas_sum = gamma_g + gamma_u
#         # Probabilities for each element
#         gamma_g = gamma_g / gammas_sum
#         gamma_u = gamma_u / gammas_sum
#         # Save for breaking out
#         likelihood = np.sum(np.log(gammas_sum))

#         # maximization_step
#         N_k = gamma_g.sum(0)
#         # Mean
#         mu = (gamma_g[:, np.newaxis] * xy).sum(0) / N_k
#         # Covariance
#         cov = np.zeros((N1, N1))
#         for j in range(N0):
#             diff = (xy[j] - mu).reshape(-1, 1)
#             cov += gamma_g[j] * np.dot(diff, diff.T)
#         cov /= N_k
#         # Weight for the Gaussian distribution
#         pi_g = N_k / N0
#         # Weight for the uniform distribution
#         pi_u = np.sum(gamma_u, axis=0) / N0

#         # Likelihood convergence check
#         if abs(likelihood - lkl_old) / likelihood < .1:
#             nstable += 1
#         if likelihood > lkl_old:
#             lkl_old = likelihood
#         if nstable == int(stable_per * n_epochs):
#             # Converged. Breaking
#             break

#     # Extract probabilities associated to the 2D Gaussian
#     gumm_p = np.array(list(gamma_g.flatten()))

#     return gumm_p


# def GUMMProbCut(gumm_p):
#     """
#     """
#     # Select the probability cut.
#     # Create the percentiles (/100.) vs probabilities array.
#     percentiles = np.arange(.01, .99, .01)
#     perc_probs = np.array([
#         percentiles, np.percentile(gumm_p, percentiles * 100.)]).T

#     # Find 'elbow' where the probabilities start climbing from ~0.
#     prob_cut = rotate(perc_probs)

#     return prob_cut


# def rotate(data):
#     """
#     Rotate a 2d vector.

#     (Very) Stripped down version of the great 'kneebow' package, by Georg
#     Unterholzner. Source:

#     https://github.com/georg-un/kneebow

#     data   : 2d numpy array. The data that should be rotated.
#     return : probability corresponding to the elbow.
#     """
#     # The angle of rotation in radians.
#     theta = np.arctan2(
#         data[-1, 1] - data[0, 1], data[-1, 0] - data[0, 0])

#     # Rotation matrix
#     co, si = np.cos(theta), np.sin(theta)
#     rot_matrix = np.array(((co, -si), (si, co)))

#     # Rotate data vector
#     rot_data = data.dot(rot_matrix)

#     # Find elbow index
#     elbow_idx = np.where(rot_data == rot_data.min())[0][0]

#     # Adding a small percentage to the selected probability improves the
#     # results by making the cut slightly stricter.
#     prob_cut = data[elbow_idx][1] # + 0.05

#     return prob_cut


if __name__ == '__main__':
    main()

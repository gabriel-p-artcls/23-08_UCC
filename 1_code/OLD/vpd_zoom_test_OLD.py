
import numpy as np
import pandas as pd
from os import listdir
import matplotlib.pyplot as plt


in_path1 = "../2_pipeline/new_method/303_no_flags/256_GMM/"
in_path2 = "../0_data/383_largest/78_large_fr_flags/"
in_path3 = "/home/gabriel/Documentos/Sol/test/"
df_CG20 = pd.read_csv("../0_data/cantat_gaudin_et_al_2020/cg2020.csv")
CG20_names = np.array([_.lower() for _ in df_CG20['Name']])


def main():
    """
    """
    files1 = listdir(in_path1)
    files2 = listdir(in_path2)
    files3 = listdir(in_path3)
    files = [in_path1+_ for _ in files1] + [in_path2+_ for _ in files2] + [in_path3+_ for _ in files3]
    # files = [in_path3+_ for _ in files3]

    Nbins_ext, zoom_f_ext = (0, 50, 100), (4,)

    for _ in range(10):
        np.random.shuffle(files)

        dists_all = []
        for file in files[:50]:
            # print(file)
            fname = file.split('/')[-1][:-4].lower()
            idx = df_CG20.index[CG20_names == fname][0]
            vpd_c = (df_CG20['pmRA'][idx], df_CG20['pmDE'][idx])
            cx0, cy0 = vpd_c

            data = pd.read_csv(file)
            msk = np.isnan(data['pmRA']) | np.isnan(data['pmDE'])
            data = data[~msk]
            vpd = np.array([data['pmRA'], data['pmDE']]).T

            dists_ij = np.empty((len(Nbins_ext), len(zoom_f_ext)))
            for i, Nbins_max in enumerate(Nbins_ext):
                for j, zoom_f in enumerate(zoom_f_ext):
                    cx, cy = centVPD(vpd, Nbins_max, zoom_f)
                    dist = np.sqrt((cx - cx0)**2 + (cy - cy0)**2)
                    # dist = np.nan if dist > .5 else dist
                    dists_ij[i][j] = dist
            dists_all.append(dists_ij)

        dists_all = np.array(dists_all)
        median_dists = np.median(dists_all, 0)
        std_dists = np.std(dists_all, 0)
        median_std = median_dists * std_dists
        print(median_std)
        i1, i2 = np.unravel_index(median_std.argmin(), median_std.shape)
        print(Nbins_ext[i1], zoom_f_ext[i2])

        # plt.imshow(
        #     median_std.T, extent=(
        #         min(Nbins_ext), max(Nbins_ext), min(zoom_f_ext), max(zoom_f_ext)),
        #     origin="lower")
        # plt.xlabel("Nbins")
        # plt.ylabel("zoom_f")
        # plt.show()
        # breakpoint()


def centVPD(vpd, Nbins_max, zoom_f, N_zoom=10, N_membs=1, N_dim=2):
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
    # break_flag = False
    # if vpd.shape[0] > N_membs * 4 * 2:
    #     break_flag = True

    for _ in range(N_zoom):
        x, y = vpd.T
        N_stars = len(x)

        # if break_flag and N_stars < N_membs * 4:
        #     break
        if N_stars < N_membs:
            break

        if Nbins_max == 0:
            N_bins = min(
                100, max(25, int(np.sqrt(N_stars / N_membs))))
        else:
            N_bins = Nbins_max

        H, edgx, edgy = np.histogram2d(x, y, bins=N_bins)

        # Find center coordinates
        flat_idx = H.argmax()
        cbx, cby = np.unravel_index(flat_idx, H.shape)
        cx = (edgx[cbx + 1] + edgx[cbx]) / 2.
        cy = (edgy[cby + 1] + edgy[cby]) / 2.

        # Zoom in
        rx, ry = edgx[1] - edgx[0], edgy[1] - edgy[0]
        msk = (x < (cx + zoom_f * rx))\
            & (x > (cx - zoom_f * rx))\
            & (y < (cy + zoom_f * ry))\
            & (y > (cy - zoom_f * ry))
        vpd = vpd[msk]

    return (cx, cy)


if __name__ == '__main__':
    main()

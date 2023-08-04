
import os
import pandas as pd
import numpy as np

path = "/media/gabriel/backup/gabriel/UCC/out_0712/"

# out_path = "/home/gabriel/Descargas/"
# out_path = "/media/gabriel/backup/gabriel/UCC/out_0712/"

"""
Reduce size of large files so that they don't conflict with Github's
large file size limit


The catalogue is composed of ~14000 OCs (so far). I have 8 repos to distribute
these files. The limit per repo is 1 Gb (allowing for images and notebooks etc)

This means each file must be 600 Kb in size MAXIMUM
"""

min_size, max_size = 0, 100  # In Mbs
# The larger this value, the less field stars will be removed
reduce_factor = .15


def main(N_membs_min=25):
    """
    """
    cl_paths = get_size(path)

    for i, large_file in enumerate(cl_paths):
        print(len(cl_paths) - i, large_file)

        cl_d = pd.read_parquet(large_file)

        cl_d = remove_field_stars(cl_d, N_membs_min)

        # out_file = out_path + large_file.split('/')[-1]
        cl_d.to_parquet(large_file, index=False)

        # out_file = out_file.replace('parquet', 'csv.gzip')
        # cl_d.to_csv(out_file, index=False, compression='gzip')
        # breakpoint()


def remove_field_stars(cl_d, N_membs_min):
    """
    """
    probs = cl_d['probs'].values
    msk = probs > 0.5
    Nmembs = N_membs_min
    if msk.sum() < Nmembs:
        # Select the 'N_membs' stars with the largest probabilities
        idx = np.argsort(probs)[::-1][:Nmembs]
        msk = np.full(len(probs), False)
        msk[idx] = True
    membs = cl_d[msk]

    df = membs

    # lmax_m, bmax_m = np.max(membs[['GLON', 'GLAT']].values, 0)
    # lmin_m, bmin_m = np.min(membs[['GLON', 'GLAT']].values, 0)
    # lmax_f, bmax_f = np.max(cl_d[['GLON', 'GLAT']].values, 0)
    # lmin_f, bmin_f = np.min(cl_d[['GLON', 'GLAT']].values, 0)
    # plx_max_m = np.max(membs['Plx'].values)
    # plx_min_m = np.min(membs['Plx'].values)

    # lr = lmax_f - lmax_m
    # lmax = lmax_m + lr * reduce_factor
    # lr = lmin_m - lmin_f
    # lmin = lmin_m - lr * reduce_factor

    # br = bmax_f - bmax_m
    # bmax = bmax_m + br * reduce_factor
    # br = bmin_m - bmin_f
    # bmin = bmin_m - br * reduce_factor

    # plx_max = plx_max_m * (1 + reduce_factor)
    # plx_min = plx_min_m * (1 - reduce_factor)

    # msk = (cl_d['GLON'] < lmax) & (cl_d['GLON'] > lmin) &\
    #     (cl_d['GLAT'] < bmax) & (cl_d['GLAT'] > bmin) &\
    #     (cl_d['Plx'] < plx_max) & (cl_d['Plx'] > plx_min)

    # df = cl_d[msk]

    return df


def get_size(path):
    all_sizes = []
    cl_paths = []
    for dirpath, dirnames, filenames in os.walk(path):
        for f in filenames:
            if f.endswith('parquet'):
                fp = os.path.join(dirpath, f)
                # skip if it is symbolic link
                if not os.path.islink(fp):
                    f_size_bytes = os.path.getsize(fp)
                    f_size = f_size_bytes / 1e6
                    if f_size > min_size and f_size < max_size:
                        # print(dirpath, f, f_size)
                        cl_paths.append(dirpath + '/' + f)
                        all_sizes.append(f_size)
    print(sum(all_sizes))
    import matplotlib.pyplot as plt
    plt.hist(all_sizes, 30)
    plt.show()
    breakpoint()

    return cl_paths


if __name__ == '__main__':
    main()

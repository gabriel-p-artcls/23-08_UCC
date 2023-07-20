
import numpy as np
from scipy.spatial.distance import cdist
import csv
import pandas as pd
import sys
sys.path.insert(1, '/home/gabriel/Github/UCC/add_New_DB/modules/')
import duplicate_probs

"""
Identify possible duplicates (and assign a probability) using
the positions estimated with the most likely members.
"""

# This is the catalogue *after* running the 'fastMP_process' script (that
# calls the 'call_fastMP' module), because it needs to have the positions
# estimated with the most likely members
date = "0712"
clpath = "/media/gabriel/backup/gabriel/UCC/out_" + date
UCC_cat = clpath + "/UCC_cat_20230702_out.csv"


def main(prob_cut=0.25, Nmax=3):
    """
    Assign a 'duplicate probability' for each cluster in 'df' compared to the
    rest of the listed clusters
    """
    df = pd.read_csv(UCC_cat)

    x, y = df['GLON_m'], df['GLAT_m']
    pmRA, pmDE, plx = df['pmRA_m'], df['pmDE_m'], df['plx_m']

    coords = np.array([x, y]).T
    # Find the distances to all clusters, for all clusters
    dist = cdist(coords, coords)
    pmRA, pmDE, plx = df['pmRA'], df['pmDE'], df['plx']

    print("Finding final duplicates and their probabilities...")
    dups_fnames_m, dups_probs_m = [], []
    for i, dists_i in enumerate(dist):

        # Only process relatively close clusters
        rad_max = Nmax * duplicate_probs.max_coords_rad(plx[i])
        msk_rad = dists_i <= rad_max
        idx_j = np.arange(0, len(dists_i))
        dists_i_msk = dists_i[msk_rad]
        j_msk = idx_j[msk_rad]

        dups_fname_i, dups_prob_i = [], []
        for k, d_ij in enumerate(dists_i_msk):
            j = j_msk[k]
            # Skip itself
            if j == i:
                continue

            dup_prob = duplicate_probs.run(x, y, pmRA, pmDE, plx, i, j)
            if dup_prob >= prob_cut:
                # Store just the first fname
                dups_fname_i.append(df['fnames'][j].split(';')[0])
                dups_prob_i.append(str(dup_prob))

        if dups_fname_i:
            dups_fname_i = ";".join(dups_fname_i)
            dups_prob_i = ";".join(dups_prob_i)
        else:
            dups_fname_i, dups_prob_i = 'nan', 'nan'

        dups_fnames_m.append(dups_fname_i)
        dups_probs_m.append(dups_prob_i)

    df['dups_fnames_m'] = dups_fnames_m
    df['dups_probs_m'] = dups_probs_m
    df.to_csv(
        UCC_cat, na_rep='nan', index=False, quoting=csv.QUOTE_NONNUMERIC)
    print(f"File {UCC_cat} updated")


if __name__ == '__main__':
    main()

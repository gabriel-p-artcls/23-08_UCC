
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# date = "0702"
# date = "0710"
date = "0712"
clpath = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/"
final_dbs_path = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/UCC_cat_20230702_out.csv"

hunt23_DB_path = "/home/gabriel/Github/UCC/add_New_DB/databases/HUNT23.csv"
hunt23_membs_path = "/media/gabriel/rest/Dropbox_nosync/Papers/2023/GDR3_members/0_data/hunt23_members.parquet"


def main():
    """
    """
    # Read HUNT23 members data
    print("Reading HUNT23 members...")
    hunt23_membs = pd.read_parquet(hunt23_membs_path)
    hunt23_ids = hunt23_membs['GaiaDR3']
    hunt23_names = hunt23_membs['Name']

    final_db = pd.read_csv(final_dbs_path)
    hunt23_DB = pd.read_csv(hunt23_DB_path)

    print("Processing clusters...")
    ours_hist, hunt23_hist, match_hist = [np.zeros(14) for _ in range(3)]
    for i, row in final_db.iterrows():

        if 'HUNT23' not in row['DB']:
            continue
        fname0 = row['fnames'].split(';')[0]

        # Read fastMP data for this cluster
        Qfold = row['quad']
        try:
            cl_d = pd.read_parquet(
                clpath + Qfold + '/datafiles/' + fname0 + '.parquet')
            probs = cl_d['probs'].values
            msk = probs > 0.5
            Nmembs = 25 # int(row['N_membs'])
            if msk.sum() < Nmembs:
                # Select the 'N_membs' stars with the largest probabilities
                idx = np.argsort(probs)[::-1][:Nmembs]
                msk = np.full(len(probs), False)
                msk[idx] = True
            cl_d = cl_d[msk]
        except:
            print(f"Could not find {Qfold}/{fname0}.parquet")
            continue

        # Identify cluster in HUNT23 DB
        dbs_idx = row['DB_i'].split(';')
        h23_j = row['DB'].split(';').index('HUNT23')
        hunt23_idx = int(dbs_idx[h23_j])
        hunt23_name = hunt23_DB['Name'][hunt23_idx].strip()

        # Identify the members for this cluster in HUNT23
        if hunt23_name.startswith('VDBH_'):
            hunt23_name = 'BH_' + hunt23_name.split(',')[0].split('_')[1]
        if hunt23_name.startswith('VDB_'):
            hunt23_name = 'vdBergh_' + hunt23_name.split(',')[0].split('_')[1]
        msk = hunt23_names == hunt23_name
        if msk.sum() == 0:
            print(f"Could not find cluster {hunt23_name} in HUNT23")
            continue

        # # box_s = np.percentile(hunt23_membs['ra'][msk], 95)-np.percentile(hunt23_membs['ra'][msk], 5)
        # # dec = np.median(hunt23_membs['dec'][msk].values)
        # # cos_dec = np.cos(np.deg2rad(dec))
        # # box_si = cos_dec * box_s / np.sqrt(2)
        # # print(fname0, np.median(hunt23_membs['b'][msk]), np.median(hunt23_membs['parallax'][msk]), box_si)
        # continue

        # Store HUNT23 magnitudes binned distribution
        hunt23_gmag = hunt23_membs['Gmag'][msk].values
        hunt23_hist += np.histogram(hunt23_gmag, bins=range(6, 21))[0]

        # Store our magnitudes binned distribution
        ours_gmag = cl_d['Gmag'].values
        ours_hist += np.histogram(ours_gmag, bins=range(6, 21))[0]

        # Store the matching magnitudes binned distribution
        idx1 = np.where(np.isin(cl_d['Source'], hunt23_ids[msk]))[0]
        match_gmag = ours_gmag[idx1]
        # This is equivalent:
        # idx2 = np.where(np.isin(hunt23_ids[msk], cl_d['Source']))[0]
        # hunt23_gmag[idx2]
        match_hist += np.histogram(match_gmag, bins=range(6, 21))[0]

        N_ours, N_hunt, N_match = ours_gmag.size, hunt23_gmag.size,\
            match_gmag.size
        if N_ours+N_hunt > 0:
            rel_diff = round((N_ours-N_hunt)/(N_ours+N_hunt), 3)
        else:
            rel_diff = 'nan'
        print("{},{},{},{},{},{}".format(
            i, fname0, N_ours, N_hunt, N_match, rel_diff))

    print(match_hist)
    # print(ours_hist)
    print(hunt23_hist)

    x = np.arange(6.5, 20.5, 1)
    plt.subplot(121)
    plt.plot(x, match_hist/ours_hist, c='b')
    plt.plot(x, match_hist/hunt23_hist, c='r')
    plt.subplot(122)
    plt.plot(x, (ours_hist - hunt23_hist) / ours_hist, marker='o', c='b')
    plt.plot(x, (ours_hist - hunt23_hist) / hunt23_hist, marker='o', c='r')
    plt.show()
    breakpoint()


if __name__ == '__main__':
    # plt.style.use('science')
    main()

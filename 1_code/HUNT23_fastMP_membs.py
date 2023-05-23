
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


clpath = "/home/gabriel/Descargas/out/"
final_dbs_path = "/home/gabriel/Github/UCC/add_New_DB/UCC_cat_20230520.csv"
hunt23_DB_path = "/home/gabriel/Github/UCC/add_New_DB/databases/HUNT23.csv"
hunt23_membs_path = "../0_data/hunt23_members.csv.gz"


def main():
    """
    """
    # Read HUNT23 members data
    print("Reading HUNT23 members...")
    hunt23_membs = pd.read_csv(hunt23_membs_path)
    hunt23_ids = hunt23_membs['source_id']
    hunt23_names = hunt23_membs['name']

    final_db = pd.read_csv(final_dbs_path)
    hunt23_DB = pd.read_csv(hunt23_DB_path)

    print("Processing clusters...")
    ours_hist, hunt23_hist, match_hist = [np.zeros(14) for _ in range(3)]
    for i, row in final_db.iterrows():

        if 'HUNT23' not in row['DB']:
            continue
        fname0 = row['fnames'].split(';')[0]
        # if 'vdb' not in row['fnames']:
        #     continue

        # Read fastMP data for this cluster
        Qfold = row['quad']
        try:
            cl_d = pd.read_csv(
                clpath + Qfold + '/datafiles/' + fname0 + '.csv.gz')
            probs = cl_d['probs'].values
            msk = probs > 0.5
            Nmembs = int(row['Nmembs'])
            if msk.sum() < Nmembs:
                # Select the 'N_membs' stars with the largest probabilities
                idx = np.argsort(probs)[::-1][:Nmembs]
                msk = np.full(len(probs), False)
                msk[idx] = True
            cl_d = cl_d[msk]
        except:
            print(f"Could not find {Qfold}/{fname0}.csv.gz file")
            continue

        # Identify cluster in HUNT23 DB
        dbs_idx = row['DB_i'].split(';')
        h23_j = row['DB'].split(';').index('HUNT23')
        hunt23_idx = int(dbs_idx[h23_j])
        hunt23_name = hunt23_DB['name'][hunt23_idx]
        # Identify the members for this cluster in HUNT23

        if hunt23_name.startswith('VDBH_'):
            hunt23_name = 'BH_' + hunt23_name.split(',')[0].split('_')[1]
        if hunt23_name.startswith('VDB_'):
            hunt23_name = 'vdBergh_' + hunt23_name.split(',')[0].split('_')[1]

        msk = hunt23_names == hunt23_name
        if msk.sum() == 0:
            print(f"Could not find cluster {hunt23_name} in HUNT23")
            continue
        # Store HUNT23 magnitudes binned distribution
        hunt23_gmag = hunt23_membs['phot_g_mean_mag'][msk].values
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
            print("{},{},{},{},{},{}".format(
                i, fname0, N_ours, N_hunt, N_match,
                round((N_ours-N_hunt)/(N_ours+N_hunt), 3)))
        else:
            print("{},{},{},{},{}".format(i, fname0, N_ours, N_hunt, 'nan'))

    print(match_hist)
    print(ours_hist)
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

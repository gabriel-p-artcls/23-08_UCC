
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# date = "0702"
# date = "0710"
date = "0712"
clpath = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/"
final_dbs_path = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/UCC_cat_20230702_out.csv"

cantat20_DB_path = "/home/gabriel/Github/UCC/add_New_DB/databases/CANTAT20.csv"
cantat20_membs_path = "/media/gabriel/rest/Dropbox_nosync/Papers/2023/GDR3_members/0_data/CG_2020_members.csv.gz"


def main():
    """
    """
    # Read cantat20 members data
    print("Reading cantat20 members...")
    cantat20_membs = pd.read_csv(cantat20_membs_path)
    cantat20_ids = cantat20_membs['GaiaDR2']
    cantat20_names = cantat20_membs['Cluster']

    final_db = pd.read_csv(final_dbs_path)
    cantat20_DB = pd.read_csv(cantat20_DB_path)

    print("Processing clusters...")
    ours_hist, cantat20_hist, match_hist = [np.zeros(12) for _ in range(3)]
    for i, row in final_db.iterrows():

        if 'CANTAT20' not in row['DB']:
            continue
        fname0 = row['fnames'].split(';')[0]

        # if fname0 != 'ruprecht143':
        #     continue

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
            msk = cl_d['Gmag'].values < 18 # IMPORTANT
            cl_d = cl_d[msk]            
        except:
            print(f"Could not find {Qfold}/{fname0}.parquet")
            continue

        # Identify cluster in cantat20 DB
        dbs_idx = row['DB_i'].split(';')
        h23_j = row['DB'].split(';').index('CANTAT20')
        cantat20_idx = int(dbs_idx[h23_j])
        cantat20_name = cantat20_DB['Cluster'][cantat20_idx].strip()

        # Identify the members for this cluster in cantat20
        if cantat20_name.startswith('FoF_'):
            cantat20_name = cantat20_name.replace('FoF', 'LP')
        if cantat20_name.startswith('VDBH_'):
            cantat20_name = 'BH_' + cantat20_name.split(',')[0].split('_')[1]
        if cantat20_name.startswith('VDB_'):
            cantat20_name = 'vdBergh_' + cantat20_name.split(',')[0].split('_')[1]
        msk = cantat20_names == cantat20_name
        if msk.sum() == 0:
            print(f"Could not find cluster {cantat20_name} in cantat20")
            continue

        # Store CANTAT20 magnitudes binned distribution
        cantat20_gmag = cantat20_membs['Gmag'][msk].values
        cantat20_hist += np.histogram(cantat20_gmag, bins=range(6, 19))[0]

        # Store our magnitudes binned distribution
        ours_gmag = cl_d['Gmag'].values
        ours_hist += np.histogram(ours_gmag, bins=range(6, 19))[0]

        # Store the matching magnitudes binned distribution
        idx1 = np.where(np.isin(cl_d['Source'], cantat20_ids[msk]))[0]
        match_gmag = ours_gmag[idx1]
        # This is equivalent:
        # idx2 = np.where(np.isin(cantat20_ids[msk], cl_d['Source']))[0]
        # cantat20_gmag[idx2]
        match_hist += np.histogram(match_gmag, bins=range(6, 19))[0]

        N_ours, N_hunt, N_match = ours_gmag.size, cantat20_gmag.size,\
            match_gmag.size
        if N_ours+N_hunt > 0:
            rel_diff = round((N_ours-N_hunt)/(N_ours+N_hunt), 3)
        else:
            rel_diff = 'nan'
        print("{},{},{},{},{},{}".format(
            i, fname0, N_ours, N_hunt, N_match, rel_diff))

    print(match_hist)
    # print(ours_hist)
    print(cantat20_hist)

    x = np.arange(6.5, 18.5, 1)
    plt.subplot(121)
    plt.plot(x, match_hist/ours_hist, c='b')
    plt.plot(x, match_hist/cantat20_hist, c='r')
    plt.subplot(122)
    plt.plot(x, (ours_hist - cantat20_hist) / ours_hist, marker='o', c='b')
    plt.plot(x, (ours_hist - cantat20_hist) / cantat20_hist, marker='o', c='r')
    plt.show()
    breakpoint()


if __name__ == '__main__':
    # plt.style.use('science')
    main()

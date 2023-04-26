
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


clpath = "/home/gabriel/Descargas/out_1/"
final_dbs_path = "../2_pipeline/3_final_DB.csv"
hunt23_DB_path = "/home/gabriel/Github/web_sites/UCC/datafiles/databases/HUNT23.csv"
hunt23_membs_path = "../0_data/hunt23_members.csv.gz"

# import os
# files = os.listdir(clpath)
# for index, file in enumerate(files):
#     new_name = file.replace('_', '').replace(' ', '').replace('-', '').lower()
#     os.rename(os.path.join(clpath, file), os.path.join(clpath, new_name))
# breakpoint()


def main():
    """
    """
    # Read HUNT23 members data
    print("Reading HUNT23 members...")
    hunt23_membs = pd.read_csv(hunt23_membs_path)
    hunt23_ids = hunt23_membs['source_id']
    hunt23_names = hunt23_membs['name']
    hunt23_mag = hunt23_membs['phot_g_mean_mag']

    cluster_plot(hunt23_membs, "Blanco_1", "")
    breakpoint()

    final_db = pd.read_csv(final_dbs_path)
    hunt23_DB = pd.read_csv(hunt23_DB_path)
    hunt23_DB_names = [_.replace('_', '').replace(' ', '').replace(
        '-', '').lower() for _ in hunt23_DB['name']]

    print("Processing clusters...")
    ours_hist, hunt23_hist, match_hist = [np.zeros(14) for _ in range(3)]
    for i, cl in enumerate(final_db['fnames']):

        # if cl != 'ngc2516':
        #     continue
        # if i > 500:
        #     break

        if '26' not in final_db['DB'][i]:
            continue

        j = -1
        for name in final_db['ID'][i].split(','):
            name = name.replace(
                '_', '').replace(' ', '').replace('-', '').lower()
            try:
                j = hunt23_DB_names.index(name)
            except ValueError:
                pass
        if j >= 0:
            hunt23_name = hunt23_DB['name'][j]
        else:
            print(f"Could not find {cl} in HUNT23")

        # Read fastMP data for this cluster
        try:
            # cl_d = pd.read_csv(clpath + cl + '.csv.gz')
            cl_d = pd.read_csv(clpath + cl + '.csv')
        except:
            print(f"Could not find {cl} gz file")
            continue

        print(i, cl)

        # Identify the members for this cluster in HUBT23
        msk = hunt23_names == hunt23_name

        # Store our magnitudes binned distribution
        ours_gmag = cl_d['Gmag'].values
        ours_hist += np.histogram(ours_gmag, bins=range(6, 21))[0]

        # Store HUNT23 magnitudes binned distribution
        hunt23_gmag = hunt23_mag[msk].values
        hunt23_hist += np.histogram(hunt23_gmag, bins=range(6, 21))[0]

        # Store the matching magnitudes binned distribution
        idx1 = np.where(np.isin(cl_d['Source'], hunt23_ids[msk]))[0]
        match_gmag = ours_gmag[idx1]
        # This is equivalent:
        # idx2 = np.where(np.isin(hunt23_ids[msk], cl_d['Source']))[0]
        # hunt23_gmag[idx2]
        match_hist += np.histogram(match_gmag, bins=range(6, 21))[0]
        
    x = np.arange(6.5, 20.5, 1)
    plt.subplot(121)
    plt.plot(x, match_hist/ours_hist, c='b')
    plt.plot(x, match_hist/hunt23_hist, c='r')
    ax = plt.subplot(122)
    plt.plot(x, (ours_hist - hunt23_hist) / ours_hist.sum(), marker='o')
    # ax.set_yscale('log')
    plt.show()

    # plot(df, 'NGC_2516', 'Theia_613')

    breakpoint()


def cluster_plot(df, clname1, clname2):
    """
    """
    msk1 = df['name'] == clname1
    clust1 = df[msk1]
    msk2 = df['name'] == clname2
    clust2 = df[msk2]

    plt.suptitle(f"N1={msk1.sum()}, N2={msk2.sum()}")
    plt.subplot(221)
    plt.scatter(clust1['l'], clust1['b'], alpha=.5)
    plt.scatter(clust2['l'], clust2['b'], alpha=.5)
    plt.subplot(222)
    plt.scatter(clust1['pmra'], clust1['pmdec'], alpha=.5)
    plt.scatter(clust2['pmra'], clust2['pmdec'], alpha=.5)
    plt.subplot(223)
    plt.hist(clust1['parallax'], alpha=.5, density=True)
    plt.hist(clust2['parallax'], alpha=.5, density=True)
    plt.subplot(224)
    plt.scatter(clust1['bp_rp'], clust1['phot_g_mean_mag'], alpha=.5)
    plt.scatter(clust2['bp_rp'], clust2['phot_g_mean_mag'], alpha=.5)
    plt.gca().invert_yaxis()
    plt.show()


if __name__ == '__main__':
    # plt.style.use('science')
    main()

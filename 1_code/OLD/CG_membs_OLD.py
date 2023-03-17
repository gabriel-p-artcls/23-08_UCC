
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


r50_mult_factor = 5


def main():
    """
    Helper script. Process the file from "Painting a portrait..",
    Cantat-Gaudin et al (2020; https://doi.org/10.1051/0004-6361/202038192)
    with all the members down to P=0.7.
    """

    # Read file data will all the members
    data_Pmemb = pd.read_csv(
        "../0_data/cantat_gaudin_et_al_2020/CG_2020_members.csv")
    data_all_cls = pd.read_csv("../0_data/cantat_gaudin_et_al_2020/cg2020.csv")

    cl = "NGC_6791"

    msk = data_all_cls['Name'] == cl
    pd.set_option('display.max_columns', None)
    print(data_all_cls[msk])
    msk = data_Pmemb['Cluster'] == cl
    data = data_Pmemb[msk]
    plt.subplot(221)
    plt.title("N={}".format(msk.sum()))
    plt.scatter(data['GLON'], data['GLAT'])
    plt.subplot(222)
    plt.scatter(data['pmRA*'], data['pmDE'])
    plt.subplot(223)
    plt.hist(data['Plx'])
    plt.subplot(224)
    plt.scatter(data['BP-RP'], data['Gmag'])
    plt.gca().invert_yaxis()
    plt.show()
    breakpoint()

    # Some stats
    CG2020Stats(data_all_cls)

    unq_clusts = np.array(list(set(data_Pmemb['Cluster'])))
    print(len(unq_clusts))

    # clust_dont_proc, data_large = clustIsolate(
    #     data_Pmemb, data_all_cls, data_large, unq_clusts)

    # iniFileGen(data_large, clust_dont_proc)


def clustIsolate(data_Pmemb, data_all_cls, data_large, unq_clusts):
    """
    """
    large_cl_names = list(data_large['Name'])

    # N_unq_clusts = []
    clust_dont_proc, box_s = [], {}
    for cl in unq_clusts:

        if cl not in large_cl_names:
            continue

        msk = data_Pmemb['Cluster'] == cl

        # r_lon = data_Pmemb[msk]['GLON'].max() - data_Pmemb[msk]['GLON'].min()
        # r_lat = data_Pmemb[msk]['GLAT'].max() - data_Pmemb[msk]['GLAT'].min()

        r_ra = data_Pmemb[msk]['RA_ICRS'].max()\
            - data_Pmemb[msk]['RA_ICRS'].min()
        r_de = data_Pmemb[msk]['DE_ICRS'].max()\
            - data_Pmemb[msk]['DE_ICRS'].min()

        msk = data_all_cls['Name'] == cl
        r50 = data_all_cls[msk]['r50'].values[0]

        flag = ''
        if (r50 > 2):
            flag += '1'
        if max(r_ra, r_de) > 5 and max(r_ra, r_de) < 100:
            flag += '2'
        # if max(r_lon, r_lat) > 5 and max(r_lon, r_lat) < 100:
        #     flag += '3'
        if flag != '':
            print(
                cl, r50,
                # round(r_lon, 1), round(r_lat, 1), round(max(r_lon, r_lat), 1),
                round(r_ra, 2), round(r_de, 2), round(max(r_ra, r_de), 2),
                flag)
            clust_dont_proc.append(cl)
            box_s[cl] = max(r50 * r50_mult_factor, max(r_ra, r_de))

        # N_unq_clusts.append(msk.sum())

    box_s_order = []
    for cl in data_large['Name']:
        try:
            box_s_order.append(box_s[cl])
        except KeyError:
            box_s_order.append(np.nan)
    data_large['box'] = box_s_order

    # for N in (50, 100, 200):
    #     msk = N_unq_clusts < N
    #     print(N, msk.sum() / N_unq_clusts.size)

    return clust_dont_proc, data_large


def iniFileGen(data_large, clust_dont_proc):
    """
    Generate the 'clusters.ini' input file for the 'GaiaEDR3_db_query' script
    that generates the frames.
    """
    data_out = pd.DataFrame()
    data_out['Name'] = data_large['Name']
    data_out['ra'] = data_large['RA']
    data_out['dec'] = data_large['DEC']
    data_out['box'] = data_large['box']
    d_pc = data_large['D_pc'].values
    msk = d_pc == '---'
    d_pc[msk] = np.nan
    d_pc = d_pc.astype(float)
    data_out['plx_d'] = 1000 / d_pc
    data_out = data_out.sort_values('box', ascending=False)

    # Store only the requested files
    # Don't store those in 'clust_dont_proc'
    # data_out = data_out[~data_out['Name'].isin(clust_dont_proc)]
    # Store those in 'clust_dont_proc'
    data_out = data_out[data_out['Name'].isin(clust_dont_proc)]
    data_out.to_csv("clusters.ini", index=False)


def CG2020Stats(data):
    """
    """
    msk1 = data['r50'] > 0.5
    print(f"N r50 > 0.5: {msk1.sum()}")

    d_pc = data['D_pc'].values
    msk = d_pc == '---'
    d_pc[msk] = -1
    d_pc = d_pc.astype(float)
    msk2 = (d_pc > 0) & (d_pc < 1000)
    print(f"d_pc < 1 Kpc: {msk2.sum()}")

    msk3 = msk1 | msk2
    print(f"N (r50 > 0.5) OR (D_pc < 1 Kpc): {msk3.sum()}")
    print("")

    msk4 = data[~msk3]['r50'] >= 0.1
    N_d = (d_pc[~msk3][msk4] > 0).sum()
    print(f"N r50 >= 0.1: {msk4.sum()}")
    print(f"No dist assigned: {msk4.sum()-N_d}")
    msk41 = (d_pc[~msk3][msk4] > 0) & (d_pc[~msk3][msk4] <= 2000)
    print(f"d_pc <= 2 Kpc: {msk41.sum()}")
    msk42 = (d_pc[~msk3][msk4] > 2000)
    print(f"d_pc > 2 Kpc: {msk42.sum()}")
    print("")

    msk5 = (data[~msk3]['r50'] < 0.1) & (data[~msk3]['r50'] >= 0.05)
    N_d = (d_pc[~msk3][msk5] > 0).sum()
    print(f"N 0.05<=r50<0.1: {msk5.sum()}")
    print(f"No dist assigned: {msk5.sum()-N_d}")
    msk51 = (d_pc[~msk3][msk5] > 0) & (d_pc[~msk3][msk5] <= 2000)
    print(f"d_pc <= 2 Kpc: {msk51.sum()}")
    msk52 = (d_pc[~msk3][msk5] > 2000)
    print(f"d_pc > 2 Kpc: {msk52.sum()}")
    print("")

    msk6 = data[~msk3]['r50'] < 0.05
    N_d = (d_pc[~msk3][msk6] > 0).sum()
    print(f"N r50<0.05: {msk6.sum()}")
    print(f"No dist assigned: {msk6.sum()-N_d}")
    msk61 = (d_pc[~msk3][msk6] > 0) & (d_pc[~msk3][msk6] <= 2000)
    print(f"d_pc <= 2 Kpc: {msk61.sum()}")
    msk62 = (d_pc[~msk3][msk6] > 2000)
    print(f"d_pc > 2 Kpc: {msk62.sum()}")


if __name__ == '__main__':
    main()

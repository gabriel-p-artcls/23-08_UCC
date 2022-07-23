
# from os import listdir
import pandas as pd
# import numpy as np

prob_min = 0.5
max_mag = 18

# not processed with pyUPMASK 'LP_5.csv'
files = (
    'IC_2602.csv', 'UPK_545.csv', 'Collinder_132.csv', 'RSG_8.csv',
    'UPK_495.csv', 'UPK_585.csv', 'UBC_8.csv', 'UPK_442.csv',
    'Gulliver_9.csv', 'UPK_552.csv', 'Alessi_2.csv', 'UPK_230.csv',
    'Stock_2.csv', 'COIN-Gaia_13.csv', 'ASCC_41.csv', 'UBC_31.csv',
    'NGC_1039.csv', 'ASCC_127.csv', 'NGC_2451A.csv')


"""
Compare our outputs with CG2020
"""

path0 = "../0_data/cantat_gaudin_et_al_2020/CG_2020_members.csv"

# pyUPMASK results
path1 = "../2_pipeline/pyUPMASK/78_large_fr_flags/"
sep1 = " "

# New method results
path2 = "../2_pipeline/new_method/78_large_fr_flags/"
sep2 = ","

# files = listdir(path2)
# np.random.shuffle(files)

for file in files:

    # CG2020
    data_Pmemb = pd.read_csv(path0)
    msk0 = data_Pmemb['Cluster'] == file.split('.')[0]
    msk0 = msk0 & (data_Pmemb['Gmag'][msk0] < max_mag)
    data0 = data_Pmemb[msk0]
    N_CG2020 = len(data0)

    data01 = data0[data0['Gmag'] <= 16]
    CG2020_IDs1 = list(data01['GaiaDR2'].values.astype(str))
    N_CG2020_1 = len(CG2020_IDs1)
    data02 = data0[data0['Gmag'] > 16]
    CG2020_IDs2 = list(data02['GaiaDR2'].values.astype(str))
    N_CG2020_2 = len(CG2020_IDs2)

    lonlat_std0 = max(data0['GLON'].std(), data0['GLAT'].std())
    pm_std0 = max(data0['pmRA*'].std(), data0['pmDE'].std())
    plx_std0 = data0['Plx'].std()

    # # pyUPMASK
    # data1 = pd.read_csv(path1 + file, sep=sep1)
    # msk1 = (data1['probs_final'] > 0.7) & (data1['Gmag'] < max_mag)

    # msk11 = data1['Gmag'][msk1] <= 16
    # our_IDs11 = list(data1['EDR3Name'][msk1][msk11])
    # our_IDs11 = [_.replace('Gaia EDR3 ', '') for _ in our_IDs11]
    # perc_missed11 = 100 * (
    #     1 - len(set(CG2020_IDs1) & set(our_IDs11)) / N_CG2020_1)

    # msk12 = data1['Gmag'][msk1] > 16
    # our_IDs12 = list(data1['EDR3Name'][msk1][msk12])
    # our_IDs12 = [_.replace('Gaia EDR3 ', '') for _ in our_IDs12]
    # perc_missed12 = 100 * (
    #     1 - len(set(CG2020_IDs2) & set(our_IDs12)) / N_CG2020_2)

    # N_stars_18_1 = msk1.sum()
    # perc_diff = 100 * (N_CG2020 - N_stars_18_1) / N_CG2020
    # lonlat_std1 = max(data1['GLON'][msk1].std(), data1['GLAT'][msk1].std())
    # pm_std1 = max(data1['pmRA'][msk1].std(), data1['pmDE'][msk1].std())
    # plx_std1 = data1['Plx'][msk1].std()
    # print("{},py,{:.1f},{:.1f},{:.1f},{:.2f},{:.2f},{:.2f}".format(
    #     file, perc_diff, perc_missed11, perc_missed12,
    #     lonlat_std0 - lonlat_std1, pm_std0 - pm_std1,
    #     plx_std0 - plx_std1))

    # New method
    NM1, NM2 = 65, 65
    for mi in range(NM1, NM2 + 1):
        data2 = pd.read_csv(path2 + "M" + str(mi) + "/" + file, sep=sep2)
        msk2 = (data2['probs_final'] > prob_min)\
            & (data2['Gmag'] < max_mag)

        msk21 = data2['Gmag'][msk2] <= 16
        our_IDs21 = list(data2['EDR3Name'][msk2][msk21])
        our_IDs21 = [_.replace('Gaia EDR3 ', '') for _ in our_IDs21]
        perc_missed21 = 100 * (
            1 - len(set(CG2020_IDs1) & set(our_IDs21)) / N_CG2020_1)

        msk22 = data2['Gmag'][msk2] > 16
        our_IDs22 = list(data2['EDR3Name'][msk2][msk22])
        our_IDs22 = [_.replace('Gaia EDR3 ', '') for _ in our_IDs22]
        perc_missed22 = 100 * (
            1 - len(set(CG2020_IDs2) & set(our_IDs22)) / N_CG2020_2)

        N_stars_18_2 = msk2.sum()
        perc_diff = 100 * (N_CG2020 - N_stars_18_2) / N_CG2020
        lonlat_std2 = max(
            data2['GLON'][msk2].std(), data2['GLAT'][msk2].std())
        pm_std2 = max(
            data2['pmRA'][msk2].std(), data2['pmDE'][msk2].std())
        plx_std2 = data2['Plx'][msk2].std()
        print("{},{},{:.1f},{:.1f},{:.1f},{:.2f},{:.2f},{:.2f}".format(
            file, "M" + str(mi), perc_diff,
            perc_missed21, perc_missed22, lonlat_std0 - lonlat_std2,
            pm_std0 - pm_std2, plx_std0 - plx_std2))

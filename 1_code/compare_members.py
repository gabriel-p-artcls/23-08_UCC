
from os import listdir
import pandas as pd
import matplotlib.pyplot as plt


def main():
    """
    Compare our outputs with CG2020
    """

    # pyUPMASK results
    # path1 = "../2_pipeline/pyUPMASK/303_no_flags/256_GMM/"
    # path1 = "../2_pipeline/pyUPMASK/303_no_flags/47_larger_5_Mb/33_smaller_10/"
    # path1 = "../2_pipeline/pyUPMASK/303_no_flags/47_larger_5_Mb/13_larger_10/"
    path1 = "../2_pipeline/pyUPMASK/303_no_flags/47_larger_5_Mb/1_larger_50/"
    sep1 = " "

    # New method results
    # path2 = "../2_pipeline/new_method/303_no_flags/256_GMM/"
    # path2 = "../2_pipeline/new_method/303_no_flags/47_larger_5_Mb/33_smaller_10/"
    # path2 = "../2_pipeline/new_method/303_no_flags/47_larger_5_Mb/13_larger_10/"
    path2 = "../2_pipeline/new_method/303_no_flags/47_larger_5_Mb/1_larger_50/"
    sep2 = ","

    files = listdir(path1)

    print(r"File   % recov pyUPMASK     % recov new")
    for file in files:
        # if file != "Alessi_44.csv":
        #     continue

        # CG2020
        data_Pmemb = pd.read_csv(
            "../0_data/cantat_gaudin_et_al_2020/CG_2020_members.csv")
        msk0 = data_Pmemb['Cluster'] == file.split('.')[0]
        data0 = data_Pmemb[msk0]
        CG2020_IDs = list(data0['GaiaDR2'].values.astype(str))
        N_CG2020 = len(CG2020_IDs)

        # pyUPMASK
        data1 = pd.read_csv(path1 + file, sep=sep1)
        msk1 = data1['probs_final'] > 0.5
        our_IDs1 = list(data1['EDR3Name'][msk1])
        our_IDs1 = [_.replace('Gaia EDR3 ', '') for _ in our_IDs1]
        N_stars_18_1 = (data1['Gmag'][msk1] < 18).sum()
        # N_stars_18_01 = 100 * (N_stars_18_1 - N_stars_18_0) / N_stars_18_0

        # New method
        data2 = pd.read_csv(path2 + file, sep=sep2)
        msk2 = data2['probs_final'] > 0.5
        our_IDs2 = list(data2['EDR3Name'][msk2])
        our_IDs2 = [_.replace('Gaia EDR3 ', '') for _ in our_IDs2]
        N_stars_18_2 = (data2['Gmag'][msk2] < 18).sum()
        # N_stars_18_02 = 100 * (N_stars_18_2 - N_stars_18_0) / N_stars_18_0

        perc_recovered1 = 100 * len(set(CG2020_IDs) & set(our_IDs1)) / N_CG2020
        perc_recovered2 = 100 * len(set(CG2020_IDs) & set(our_IDs2)) / N_CG2020
        print("{} ({}):     {:.0f}% ({})    {:.0f}% ({})".format(
            file, N_CG2020, perc_recovered1, N_stars_18_1, perc_recovered2,
            N_stars_18_2))

        # continue
        # data, msk = data1, msk1
        data, msk = data2, msk2

        lon, lat = data1['GLON'], data1['GLAT']
        xmin, xmax = lon.min(), lon.max()
        ymin, ymax = lat.min(), lat.max()
        plt.subplot(221)
        plt.scatter(
            data0['GLON'], data0['GLAT'], c='b', alpha=.5, label=msk0.sum())
        plt.scatter(
            data['GLON'][msk], data['GLAT'][msk], marker='x', c='r',
            alpha=.5, label=msk.sum())
        plt.legend()
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
        plt.subplot(222)
        plt.scatter(data0['pmRA*'], data0['pmDE'], c='b', alpha=.5)
        plt.scatter(
            data['pmRA'][msk], data['pmDE'][msk], marker='x', c='r', alpha=.5)
        plt.subplot(223)
        plt.hist(data0['Plx'], color='b', alpha=.5)
        plt.hist(data['Plx'][msk], color='r', alpha=.5)
        plt.subplot(224)
        plt.scatter(data0['BP-RP'], data0['Gmag'], c='b', alpha=.5)
        plt.scatter(
            data['BP-RP'][msk], data['Gmag'][msk], marker='x', c='r', alpha=.5)
        plt.gca().invert_yaxis()
        plt.show()
        # breakpoint()


if __name__ == '__main__':
    main()

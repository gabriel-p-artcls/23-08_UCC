
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
Print/plot data on the CG20 clusters.
"""

# CG2020
path0 = "../0_data/cantat_gaudin_et_al_2020/CG_2020_members.csv"
data_CG = pd.read_csv(path0)
def main():
    # func1()
    func2('Stock 11')


def func1():
    """
    """
    names_CG = np.array([_ for _ in data_CG['Cluster']])
    names_unq = list(set(names_CG))

    xyr = []
    num = []
    for name in names_unq:
        msk = names_CG == name
        data_cl = data_CG[msk]

        num.append(len(data_cl))

        xy_c = np.median([data_cl['GLON'], data_cl['GLAT']], 1)
        xyc_dist = np.sqrt(
            (xy_c[0] - data_cl['GLON'])**2 + (xy_c[1] - data_cl['GLAT'])**2)
        rad = np.percentile(xyc_dist, 95)
        frame_l = 2 * rad

        d_kpc = 1. / data_cl['Plx']
        plx_c = np.median(d_kpc)
        plxc_dist = abs(d_kpc - plx_c)
        plx_rad = np.percentile(plxc_dist, 95)

        if frame_l > 50:
            print("Large frame", name, frame_l)

        if plx_rad > 20:
            print("Large plx rad", name, plx_rad)

        xyr.append([plx_c, frame_l, plx_rad])

    print(max(num))

    xyr = np.array(xyr).T
    plt.subplot(121)
    plt.scatter(xyr[0], xyr[1])
    plt.xlabel("plx")
    plt.ylabel("Frame length")
    plt.subplot(122)
    plt.scatter(xyr[0], xyr[2])
    plt.xlabel("plx")
    plt.ylabel("Plx radius")
    plt.show()


def func2(cl="NGC_6791"):
    data_all_cls = pd.read_csv("../0_data/cantat_gaudin_et_al_2020/cg2020.csv")
    msk = data_all_cls['Name'] == cl
    # pd.set_option('display.max_columns', None)
    # print(data_all_cls[msk])
    msk = data_CG['Cluster'] == cl
    data = data_CG[msk]
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


if __name__ == '__main__':
    # plt.style.use('science')
    main()


import pandas as pd
from os import listdir
from scipy import spatial
import matplotlib.pyplot as plt
import numpy as np
import time as tm
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/home/gabriel/Github/fastmp/')
from fastmp import fastMP


in_path = "../0_data/CG_2017/"
out_path = "/home/gabriel/Documentos/fastMP_out/test/"

# already_proc = [_[:-4] for _ in listdir(out_path)]

df_CG20 = pd.read_csv("../0_data/cantat_gaudin_et_al_2020/cg2020.csv")
CG20_names = np.array([_.lower() for _ in df_CG20['Name']])
path0 = "../0_data/cantat_gaudin_et_al_2020/CG_2020_members.csv"
data_CG = pd.read_csv(path0)
names_CG = np.array([_.lower() for _ in data_CG['Cluster']])

# file_process = ()
# file_process = [_.lower() for _ in file_process]


# Close by clusters
xys = np.array([df_CG20['GLON'].values, df_CG20['GLAT'].values]).T
tree = spatial.cKDTree(xys)
inx = tree.query(xys, k=5 + 1)


def main():

    files = listdir(in_path)
    # np.random.shuffle(files)
    # files = np.sort(files)

    for i, file in enumerate(files):
        # if file[:-4] not in file_process:
        #     continue
        # if file[:-4] not in ('ngc_7024',):
        #     continue
        # if file[:-4] in already_proc:
        #     continue

        fname = file.split('/')[-1][:-4]
        # if fname != 'sigma_ori':
        #     continue
        fname_l = fname.lower()

        breakpoint()
        idx = df_CG20.index[CG20_names == fname_l][0]
        N_membs_CG = df_CG20['n07'][idx]
        xy_c = (df_CG20['GLON'][idx], df_CG20['GLAT'][idx])
        vpd_c = (df_CG20['pmRA'][idx], df_CG20['pmDE'][idx])
        plx_c = df_CG20['plx'][idx]

        # Indexes to the 5 closest clusters in XY
        ex_cls_idx = inx[1][idx][1:]
        centers_ex = []
        for ex_i in ex_cls_idx:
            ex_cl_dict = {
                'xy': [df_CG20['GLON'][ex_i], df_CG20['GLAT'][ex_i]],
                'pms': [df_CG20['pmRA'][ex_i], df_CG20['pmDE'][ex_i]],
                'plx': [df_CG20['plx'][ex_i]]
            }
            centers_ex.append(ex_cl_dict)

        # print(file)

        data = pd.read_csv(in_path + file)
        X = np.array([
            data['GLON'].values, data['GLAT'].values, data['pmRA'].values,
            data['pmDE'].values, data['Plx'].values,
            data['e_pmRA'].values, data['e_pmDE'].values,
            data['e_Plx'].values])

        # print("{:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}".format(
        #       xy_c[0], xy_c[1], vpd_c[0], vpd_c[1], plx_c))

        s = tm.time()
        probs_all, N_membs = fastMP(
        # Ns, xy_c1, vpd_c1, plx_c1, xy_c2, vpd_c2, plx_c2, t1, t2 = fastMP(
            xy_c=xy_c,
            # plx_c=plx_c,
            # vpd_c=vpd_c,
            # fixed_centers=True,
            # fix_N_clust=25,
            # centers_ex=centers_ex
        ).fit(X)

        # print("{}, {}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}".format(
        #     file, Ns, xy_c[0], xy_c[1], vpd_c[0], vpd_c[1], plx_c,
        #     xy_c1[0], xy_c1[1], vpd_c1[0], vpd_c1[1], plx_c1,
        #     xy_c2[0], xy_c2[1], vpd_c2[0], vpd_c2[1], plx_c2, t1, t2))
        # continue

        print("{} '{}': Nm_CG: {}, Nm: {}, P>0.5: {}, t={:.2f}".format(
            i, file, N_membs_CG, N_membs, (probs_all > 0.5).sum(),
            tm.time() - s))

        # Save file with IDs and probabilities
        data['probs_final'] = probs_all
        data.to_csv(out_path + '/' + fname_l + '.csv', index=False)

        plot(fname_l, out_path)


def plot(file, path):

    # df1 = pd.read_csv(path + file + '_no.csv')
    # df2 = pd.read_csv(path + file + '.csv')

    # msk1 = df1['probs_final'] > 0.5
    # data1 = df1[msk1]
    # msk2 = df2['probs_final'] > 0.5
    # data2 = df2[msk2]

    # # CG2020
    # msk0 = names_CG == file.lower()
    # data_CGm = data_CG[msk0]

    # plt.figure(figsize=(15, 9))
    # plt.subplot(231)
    # plt.scatter(data2['GLON'], data2['GLAT'], c='r',
    #             alpha=.25, label=f"dim_norm: {msk2.sum()}")
    # plt.scatter(data1['GLON'], data1['GLAT'], c='b',
    #             alpha=.25, label=f"no_dim_norm: {msk1.sum()}")
    # plt.scatter(data_CGm['GLON'], data_CGm['GLAT'], c='k', marker='x', lw=.7,
    #             label=f"CG20: {len(data_CGm)}")
    # plt.legend()

    # plt.subplot(232)
    # plt.scatter(data2['pmRA'], data2['pmDE'], c='r', alpha=.25)
    # plt.scatter(data1['pmRA'], data1['pmDE'], c='b', alpha=.25)
    # plt.scatter(
    #     data_CGm['pmRA*'], data_CGm['pmDE'], c='k', marker='x', lw=.7)

    # plt.subplot(233)
    # plt.title("dim norm")
    # plt.scatter(data2['BP-RP'], data2['Gmag'], c=data2['probs_final'],
    #             alpha=.5, vmin=.5)
    # plt.scatter(
    #     data_CGm['BP-RP'], data_CGm['Gmag'], c='k', marker='x', lw=.7,
    #     alpha=.5, zorder=3)
    # plt.colorbar()
    # plt.gca().invert_yaxis()

    # plt.subplot(234)
    # plt.hist(data2['Plx'], color='r', alpha=.5, density=True)
    # plt.hist(data1['Plx'], color='b', alpha=.5, density=True)
    # plt.hist(data_CGm['Plx'], color='k', alpha=.5, density=True)

    # plt.subplot(235)
    # # plt.hist(data2['probs_final'], color='r', alpha=.5)
    # # plt.hist(data1['probs_final'], color='b', alpha=.5)
    # pms = np.array([data2['pmRA'], data2['pmDE']]).T
    # vpd_c = np.median(pms, 0)
    # pm_rad_2 = spatial.distance.cdist(pms, np.array([vpd_c])).T[0]
    # plt.scatter(data2['Plx'], pm_rad_2, c='r', alpha=.25)
    # pms = np.array([data1['pmRA'], data1['pmDE']]).T
    # vpd_c = np.median(pms, 0)
    # pm_rad_1 = spatial.distance.cdist(pms, np.array([vpd_c])).T[0]
    # plt.scatter(data1['Plx'], pm_rad_1, c='b', alpha=.25)

    # plt.subplot(236)
    # plt.title("no dim norm")
    # plt.scatter(data1['BP-RP'], data1['Gmag'], c=data1['probs_final'],
    #             alpha=.5, vmin=.5)
    # plt.scatter(
    #     data_CGm['BP-RP'], data_CGm['Gmag'], c='k', marker='x', lw=.7,
    #     alpha=.5, zorder=3)
    # plt.colorbar()
    # plt.gca().invert_yaxis()

    # plt.suptitle(file)
    # plt.show()

    df = pd.read_csv(path + file + '.csv')
    msk = df['probs_final'] > 0.5
    data = df[msk]

    msk_gr_18 = data['Gmag'] > 18.
    msk_le_18 = data['Gmag'] <= 18.
    data_gr_18 = data[msk_gr_18]
    data_le_18 = data[msk_le_18]

    N_G_gr_18 = msk_gr_18.sum()
    N_G_lt_18 = len(data) - N_G_gr_18

    # CG2020
    msk0 = names_CG == file.lower()
    data_CGm = data_CG[msk0]

    plt.figure(figsize=(15, 9))
    plt.subplot(231)
    plt.scatter(data_le_18['GLON'], data_le_18['GLAT'], c='b',
                alpha=.25, label=f"G<=18: {N_G_lt_18}")
    plt.scatter(data_gr_18['GLON'], data_gr_18['GLAT'], c='r',
                alpha=.25, label=f"G>18: {N_G_gr_18}")
    plt.scatter(data_CGm['GLON'], data_CGm['GLAT'], c='k', marker='x', lw=.7,
                label=f"CG20: {len(data_CGm)}")
    plt.legend()

    plt.subplot(232)
    plt.scatter(data_le_18['pmRA'], data_le_18['pmDE'], c='b', alpha=.25)
    plt.scatter(data_gr_18['pmRA'], data_gr_18['pmDE'], c='r', alpha=.25)
    plt.scatter(
        data_CGm['pmRA*'], data_CGm['pmDE'], c='k', marker='x', lw=.7)

    plt.subplot(233)
    plt.scatter(data['BP-RP'], data['Gmag'], c=data['probs_final'],
                alpha=.5, vmin=.5)
    plt.scatter(
        data_CGm['BP-RP'], data_CGm['Gmag'], c='k', marker='x', lw=.7,
        alpha=.5, zorder=3)
    plt.colorbar()
    plt.gca().invert_yaxis()

    plt.subplot(234)
    plt.hist(data_le_18['Plx'], color='b', alpha=.5, density=True)
    plt.hist(data_gr_18['Plx'], color='r', alpha=.5, density=True)
    plt.hist(data_CGm['Plx'], color='k', alpha=.5, density=True)

    plt.subplot(235)
    plt.hist(data_le_18['probs_final'], color='b', alpha=.5, density=True)
    plt.hist(data_gr_18['probs_final'], color='r', alpha=.5, density=True)

    plt.suptitle(file)
    plt.show()


if __name__ == '__main__':
    # plt.style.use('science')
    main()


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Read file data will all the members
data_Pmemb = pd.read_csv(
    "../0_data/cantat_gaudin_et_al_2020/CG_2020_members.csv")
data_all_cls = pd.read_csv("../0_data/cantat_gaudin_et_al_2020/cg2020.csv")
fastMPfold = "/home/gabriel/Documentos/fastMP_out/all/"
cldata = "/media/gabriel/rest/Dropbox_nosync/Papers/2023/GDR3_members/0_data/"
folders = ('S1', 'S2', 'M1', 'L1', 'L2', 'L3')

# for i, name in enumerate(data_all_cls['Name']):
#     print(i, name)
# breakpoint()


# name = 'Juchert_10'
# data_probs = pd.read_csv(fastMPfold + name.lower() + '.csv')
# msk_p = data_probs['probs_final'] > 0.7
# for fold in folders:
#     try:
#         data = pd.read_csv(cldata + fold + '/' + name.lower() + '.csv')
#     except:
#         pass
# data = data[msk_p]
# plt.subplot(221)
# # plt.title("N={}".format(msk.sum()))
# plt.scatter(data['GLON'], data['GLAT'])
# plt.subplot(222)
# plt.scatter(data['pmRA'], data['pmDE'])
# plt.subplot(223)
# plt.hist(data['Plx'])
# plt.subplot(224)
# plt.scatter(data['BP-RP'], data['Gmag'])
# plt.gca().invert_yaxis()
# plt.show()


def readAllMembs():

    allMembs, allMags, allPlx = [[], []], [[], []], []
    for i, name in enumerate(data_all_cls['Name']):
        # if name != 'UPK_126':
        #     continue
        print(i, name)

        msk = data_Pmemb['Cluster'] == name
        if msk.sum() <= 1:
            print("*** CG <= 1 stars")
            allMembs[0].append([np.nan, np.nan, np.nan])
            allMembs[1].append([np.nan, np.nan, np.nan])
            allPlx.append([np.nan, np.nan])
            continue

        # glon, glat = data_Pmemb['GLON'][msk], data_Pmemb['GLAT'][msk]
        # pmra, pmde = data_Pmemb['pmRA*'][msk], data_Pmemb['pmDE'][msk]
        plx = data_Pmemb['Plx'][msk]
        Gmag = data_Pmemb['Gmag'][msk]

        data_probs = pd.read_csv(fastMPfold + name.lower() + '.csv')
        msk_p = data_probs['probs_final'] > 0.7
        if msk_p.sum() <= 1:
            print("*** fastMP <= 1 stars")
            allMembs[0].append([np.nan, np.nan, np.nan])
            allMembs[1].append([np.nan, np.nan, np.nan])
            allPlx.append([np.nan, np.nan])
            continue

        for fold in folders:
            try:
                data = pd.read_csv(cldata + fold + '/' + name.lower() + '.csv')
            except:
                pass
        data_m = data[msk_p]

        # Number of members
        Nmembs_CG, Nmembs = [], []
        gold = 0
        for gm in (16, 17, 18):
            Nmembs_CG.append(((Gmag > gold) & (Gmag < gm)).sum())
            Nmembs.append(
                ((data_m['Gmag'] > gold) & (data_m['Gmag'] < gm)).sum())
            gold = gm
        allMembs[0].append(Nmembs_CG)
        allMembs[1].append(Nmembs)

        # # Size of frame
        # CG_fr_s = max(np.ptp(glon), np.ptp(glat))
        # data_fr_s = max(np.ptp(data['GLON']), np.ptp(data['GLAT']))
        # allFrsize.append(data_fr_s / CG_fr_s)

        # Plx difference
        fMP_plx = np.median(data_m['Plx'])
        delta_plx = np.median(plx) - fMP_plx
        allPlx.append([fMP_plx, delta_plx])

        allMags[0] += list(Gmag)
        allMags[1] += list(data_m['Gmag'])

    allPlx_df = pd.DataFrame(
        np.array(allPlx), columns=('fMP_plx', 'delta_plx'))
    allPlx_df.to_csv('allPlx_df.csv')

    allMags_df = pd.DataFrame(allMags[0], columns=('CG',))
    allMags_df.to_csv('allMags_CG.csv')
    allMags_df = pd.DataFrame(allMags[1], columns=('fastMP',))
    allMags_df.to_csv('allMags_fMP.csv')

    allMembs = np.array(allMembs)
    membs_ratio = np.nanmedian((allMembs[1] + 1) / (allMembs[0] + 1), 0)
    print(membs_ratio)


def compareMembs():
    """
    """
    recov_ratio = []
    for i, name in enumerate(data_all_cls['Name']):
        print(i, name)

        msk = data_Pmemb['Cluster'] == name
        CG_ids = data_Pmemb['GaiaDR2'][msk]

        data_probs = pd.read_csv(fastMPfold + name.lower() + '.csv')
        msk_p = data_probs['probs_final'] > 0.7
        fMP_ids = list(data_probs['EDR3Name'][msk_p])
        fMP_ids = [int(_.replace('Gaia EDR3 ', '')) for _ in fMP_ids]
        N_recovered = len(list(set(CG_ids) & set(fMP_ids)))

        recov_ratio.append(N_recovered / len(CG_ids))

    return recov_ratio


# readAllMembs()
# recov_ratio = compareMembs()

mags_CG = pd.read_csv('allMags_CG.csv')
mags_fMP = pd.read_csv('allMags_fMP.csv')
all_plx = pd.read_csv('allPlx_df.csv')

breakpoint()

plt.style.use('science')
plt.title(r"Identified members with $P>0.7$")
plt.hist(mags_fMP['fastMP'], 25, alpha=.35, label="This work")
plt.hist(mags_CG['CG'], 25, alpha=.35, label="CG2020")
plt.legend()
plt.xlim(5.1, 19.6)
plt.xlabel("Gmag")
plt.ylabel("N stars")
plt.savefig("mags_histo.png", dpi=300)

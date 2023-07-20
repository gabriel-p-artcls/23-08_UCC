
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use('science')

"""
"""

# # USe this block to generate the 'fatMP_hist' array with the latest data
# import pandas as pd
# # date = "0702"
# date = "0712"
# clpath = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/"
# final_dbs_path = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/UCC_cat_20230702_out.csv"
# final_db = pd.read_csv(final_dbs_path)
# fastMP_gmag = []
# for i, row in final_db.iterrows():
#     if i in [500, 1000, 2000, 4000, 6000, 8000, 10000, 13000]:
#         print(i)

#     fname0 = row['fnames'].split(';')[0]
#     # Read fastMP data for this cluster
#     Qfold = row['quad']
#     try:
#         cl_d = pd.read_parquet(
#             clpath + Qfold + '/datafiles/' + fname0 + '.parquet')
#         probs = cl_d['probs'].values
#         msk = probs > 0.5
#         Nmembs = 25 #int(row['N_membs'])
#         if msk.sum() < Nmembs:
#             # Select the 'N_membs' stars with the largest probabilities
#             idx = np.argsort(probs)[::-1][:Nmembs]
#             msk = np.full(len(probs), False)
#             msk[idx] = True
#         cl_d = cl_d[msk]
#     except:
#         print(f"Could not find {Qfold}/{fname0}.parquet")
#         continue
#     fastMP_gmag += list(cl_d['Gmag'].values)

# fatMP_hist = np.histogram(fastMP_gmag, bins=range(6, 21))[0]
# print(fatMP_hist)
# breakpoint()

# # 0702, fastMP using KDE removal
# fatMP_hist = np.array([
#     366, 787, 1864, 3933, 7787, 14310, 27428, 50510,
#     86757, 132869, 184912, 207091, 160102, 71437])
# #
# h23_match_hist = np.array([
#     195., 485., 1090., 2399., 4700., 8820., 16205., 29484., 49015., 70120.,
#     88084., 85144., 52312., 13807.])
# hunt23_hist = np.array([
#     308., 740., 1621., 3543., 6878., 12533., 22679., 40560., 66774.,
#     97199., 127087., 138377., 112699., 54748.])
# c20_match_hist = np.array([
#     107., 238., 609., 1332., 2428., 5413., 10222., 18233.,
#     29487., 38388., 42490., 30385.])
# cantat20_hist = np.array([
#     151., 328.,  779., 1763., 3489., 6943., 13009., 22841.,
#     36644., 49080., 56561., 40395.])

# 0712, fastMP NOT using KDE removal
fatMP_hist = np.array([
    406, 890, 2090, 4433, 8722, 16146, 31061, 57787, 100936, 159287,
    230105, 271588, 224368, 107015])
#
h23_match_hist = np.array([
    213., 525.,  1164.,  2549. , 5004. , 9351., 17174. ,31223. ,51856. ,74482.,
    94385., 93656., 62084., 19608.])
hunt23_hist = np.array([
    308., 740., 1621., 3543., 6878., 12533., 22679., 40560., 66774.,
    97199., 127087., 138377., 112699., 54748.])
#
c20_match_hist = np.array([
    114., 261., 635., 1404., 2531., 5694., 10702., 19001., 30692., 40109.,
    44608., 32277.])
cantat20_hist = np.array([
    151., 328., 779., 1763., 3489, 6943., 13009.,
    22841., 36644., 49080., 56561., 40395.])

print(fatMP_hist.sum()/13500, hunt23_hist.sum()/6200, c20_match_hist.sum()/2000)


fig = plt.figure(figsize=(3, 8.8))
ax = plt.subplot()

x1 = np.arange(6.5, 20.5, 1)
x2 = np.arange(6.5, 18.5, 1)

ax = plt.subplot(311)
plt.plot(x1, fatMP_hist, label='UCC', lw=2.5, alpha=.75)
plt.plot(x1, hunt23_hist, label='HUNT23', lw=2.5, alpha=.75)
plt.plot(x2, c20_match_hist, label='CANTAT20', lw=2.5, alpha=.75)
plt.legend()
ax.set_xticklabels([])
plt.yscale('log')
plt.ylabel("N")
# plt.xlabel("G [mag]")

ax = plt.subplot(312)
plt.plot(x1, fatMP_hist/13500, label='UCC', lw=2.5, alpha=.75)
plt.plot(x1, hunt23_hist/6200, label='HUNT23', lw=2.5, alpha=.75)
plt.plot(x2, c20_match_hist/2000, label='CANTAT20', lw=2.5, alpha=.75)
plt.legend()
ax.set_xticklabels([])
ax.set_xticklabels([])
plt.ylabel("N (normalized)")
# plt.xlabel("G [mag]")

ax = plt.subplot(313)
plt.plot(x1, h23_match_hist/hunt23_hist, label="HUNT23", lw=2.5, alpha=.75, c='green')
plt.plot(x2, c20_match_hist/cantat20_hist, label="CANTAT20", lw=2.5, alpha=.75, c='orange')
plt.legend()
plt.ylabel(r"\% of matched members")
# plt.xlabel("G [mag]")
ax.set_xlabel("G [mag]", labelpad=-1)

plt.subplots_adjust(hspace=.025)
plt.savefig("membs_compare.png", dpi=300)

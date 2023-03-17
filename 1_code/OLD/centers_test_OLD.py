
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


"""
name,N,xcg,ycg,pmracg,pmdecg,plxcg,x1,y1,pmra1,pmde1,plx1,x2,y2,pmra2,pmde2,plx2,t1,t2

* xy_c given, NN instead of KDE
1.  N_cent=200, Ndens=5, 5D data <-- worse than OLD
2.  N_cent=200, Ndens=25, 5D data <-- not better than OLD
3.  N_cent=200, Ndens=25, 4D data (no plx) <-- worse than OLD
4.  N_cent=100, Ndens=25, 5D data <-- worse than OLD
5.  N_cent=100, Ndens=25, 5D data, cent=median(10) <-- worse than OLD
6.  N_cent=500, Ndens=5, 5D data, cent=median(10) <-- marginally better t/ OLD
7.  N_cent=500, Ndens=5, 4D data, cent=median(10) <-- equal to 6
8.  N_cent=500, Ndens=25, 4D data, cent=median(10) <-- marginally better than 6
9.  N_cent=500, Ndens=25, 5D data, cent=max(dens) <-- worse than 8
10. N_cent=500, Ndens=25, 4D data, NN_dist=median(), cent=median(10) <-- worse than 8
11. N_cent=500, Ndens=10, 4D data, NN_dist=median(), cent=median(10) <-- worse than 8
17: 8 but N_cent=1000 in both <-- no improv in either over 16
18: 8 but N_cent=100 in both <-- both worse
--
16. 8 but using N_cent=500 in OLD too <-- still slightly better than OLD
20. 16 but old method using 4D in KDE <-- not better than before


* vpd_c given, NN instead of KDE
12. 8 <-- marginally better than OLD
--
15. 12 but using N_cent=500 in OLD too <-- still slightly better than OLD
21. 15 but old method using 4D in KDE <-- not better than before


* plx_c given, NN instead of KDE
--
19: 8 <-- marginally better than OLD
22. 19 but old method using 4D in KDE <-- not better than before


* xy_c+vpd_c given, NN instead of KDE
13. 8 <-- marginally better than OLD
--
14: 13 but using N_cent=500 in OLD too <-- still slightly better than OLD
23. 14 but old method using 4D in KDE <-- not better than before

"""

for f_ID in (16, 20, 15, 21, 19, 22, 14, 23): #(16, 15, 19, 14):
    f_ID = str(f_ID)
    df = pd.read_csv(
        f"/home/gabriel/Github/fastmp/fastmp/centers_data_{f_ID}.txt")
    print('\n', f_ID)

    xy_diff1 = np.sqrt((df['xcg'] - df['x1'])**2 + (df['ycg'] - df['y1'])**2)
    xy_diff2 = np.sqrt((df['xcg'] - df['x2'])**2 + (df['ycg'] - df['y2'])**2)
    pm_diff1 = np.sqrt(
        (df['pmracg'] - df['pmra1'])**2 + (df['pmdecg'] - df['pmde1'])**2)
    pm_diff2 = np.sqrt(
        (df['pmracg'] - df['pmra2'])**2 + (df['pmdecg'] - df['pmde2'])**2)
    plx_diff1, plx_diff2 = abs(df['plxcg'] - df['plx1']),\
        abs(df['plxcg'] - df['plx2'])

    print('Old method  {:.0f}  | New method  {:.0f}'.format(
        df['t1'].sum(), df['t2'].sum()))
    print("{:.4f}, {:.4f} ({:.4f}) | {:.4f}, {:.4f} ({:.4f})".format(
        np.median(xy_diff1), np.percentile(xy_diff1, 95), max(xy_diff1),
        np.median(xy_diff2), np.percentile(xy_diff2, 95), max(xy_diff2))
    )
    print("{:.4f}, {:.4f} ({:.4f}) | {:.4f}, {:.4f} ({:.4f})".format(
        np.median(pm_diff1), np.percentile(pm_diff1, 95), max(pm_diff1),
        np.median(pm_diff2), np.percentile(pm_diff2, 95), max(pm_diff2))
    )
    print("{:.4f}, {:.4f} ({:.4f}) | {:.4f}, {:.4f} ({:.4f})".format(
        np.median(plx_diff1), np.percentile(plx_diff1, 95), max(plx_diff1),
        np.median(plx_diff2), np.percentile(plx_diff2, 95), max(plx_diff2))
    )

    continue

    plt.subplot(221)
    plt.scatter(df['N'], df['t1'], c='b', alpha=.5, label='Old ({:.0f})'.format(df['t1'].sum()))
    plt.scatter(df['N'], df['t2'], c='r', alpha=.5, label='New ({:.0f})'.format(df['t2'].sum()))
    # plt.ylabel('T_old / T_new')
    plt.legend()
    plt.subplot(222)
    plt.scatter(df['N'], xy_diff1, c='b', alpha=.5, label='Old')
    plt.scatter(df['N'], xy_diff2, c='r', alpha=.5, label='New')
    # plt.legend()
    plt.subplot(223)
    plt.scatter(df['N'], pm_diff1, c='b', alpha=.5)
    plt.scatter(df['N'], pm_diff2, c='r', alpha=.5)
    plt.subplot(224)
    plt.scatter(df['plxcg'], plx_diff1, c='b', alpha=.5)
    plt.scatter(df['plxcg'], plx_diff2, c='r', alpha=.5)
    plt.show()

# plt.subplot(321)
# plt.hist(x_diff1, 30, color='b', alpha=.5)
# plt.hist(x_diff2, 30, color='r', alpha=.5)
# plt.subplot(322)
# plt.hist(y_diff1, 30, color='b', alpha=.5)
# plt.hist(y_diff2, 30, color='r', alpha=.5)
# plt.subplot(323)
# plt.hist(pmra_diff1, 30, color='b', alpha=.5)
# plt.hist(pmra_diff2, 30, color='r', alpha=.5)
# plt.subplot(324)
# plt.hist(pmde_diff1, 30, color='b', alpha=.5)
# plt.hist(pmde_diff2, 30, color='r', alpha=.5)
# plt.subplot(325)
# plt.hist(plx_diff1, 30, color='b', alpha=.5)
# plt.hist(plx_diff2, 30, color='r', alpha=.5)
# plt.show()

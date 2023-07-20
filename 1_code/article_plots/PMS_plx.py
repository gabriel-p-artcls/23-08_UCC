
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scienceplots
plt.style.use('science')

date = "0712"


def main():
    """
    """
    data = pd.read_csv("classif_data_" + date + ".csv")

    # plt.figure(figsize=(8, 4))
    fig, (ax1, ax2) = plt.subplots(2, figsize=(12, 8))

    s1 = np.clip(np.log(data['N_ucc'])**2.5, a_min=10, a_max=500)
    msk1 = ~np.isnan(s1)
    s2 = np.clip(np.log(data['N_h23'])**2.5, a_min=10, a_max=500)
    msk2 = ~np.isnan(s2)
    s3 = np.clip(np.log(data['N_cg20'])**2.5, a_min=10, a_max=500)
    msk3 = ~np.isnan(s3)

    d_pc_ucc = np.clip(1000/data['ucc_plx'][msk1], a_min=0.01, a_max=np.inf)
    d_pc_h23 = np.clip(1000/data['h23_plx'][msk2], a_min=0.01, a_max=np.inf)
    d_pc_cg20 = np.clip(1000/data['cg20_plx'][msk3], a_min=0.01, a_max=np.inf)

    ax1.scatter(
        d_pc_ucc, data['ucc_pm_disp'][msk1], s=s1[msk1],
        alpha=.5, ec='w', lw=.5, label='UCC')
    ax1.scatter(
        d_pc_h23, data['h23_pm_disp'][msk2], s=s2[msk2],
        alpha=.5, marker='v', ec='w', lw=.5, label='HUNT23')
    ax1.scatter(
        d_pc_cg20, data['cg20_pm_disp'][msk3], s=s3[msk3],
        alpha=.5, marker='^', ec='w', lw=.5, label='CANTAT20')

    d_lim_kpc = 1 / (4.74/5/np.sqrt(2))  # 1491.786
    xmin, xmax = 40, 20000

    d_pc1 = np.arange(xmin, d_lim_kpc * 1000, 100)
    d_pc2 = np.arange(d_lim_kpc * 1000, xmax, 100)
    pm_lim1 = 5 * 1.414 / (d_pc1/1000 * 4.74)
    pm_lim2 = np.ones(len(d_pc2))
    ax1.plot(d_pc1, pm_lim1, ls=':', lw=2, c='k')
    ax1.plot(d_pc2, pm_lim2, ls=':', lw=2, c='k')

    # Fill
    y2 = np.array(list(pm_lim1) + list(pm_lim2))
    x = np.linspace(xmin, xmax, len(y2))
    y1 = np.ones(len(x)) * 1000
    ax1.fill_between(x, y1, y2, color='grey', alpha=0.075)

    ax1.annotate(
        'NGC 7826', xy=(1150, 5.2), xytext=(1500, 13), fontsize=12,
        horizontalalignment="center",
        # Custom arrow
        arrowprops=dict(arrowstyle='->', lw=.7))

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(.02, 90)
    # ax1.set_xlabel("Distance [pc]")
    ax1.set_xticklabels([])
    ax1.set_ylabel("Total PM dispersion [mas/yr]", fontsize=14)
    lgnd = ax1.legend(loc="lower left", fontsize=14)
    # lgnd = plt.legend(loc="lower left", scatterpoints=1, fontsize=10)
    for handle in lgnd.legend_handles:
        handle.set_sizes([80])

    ax2.scatter(
        d_pc_ucc, data['r_50_pc'][msk1], s=s1[msk1],
        alpha=.5, ec='w', lw=.5)
    ax2.scatter(
        d_pc_h23, data['h23_r50'][msk2], s=s2[msk2],
        alpha=.5, marker='v', ec='w', lw=.5)
    ax2.scatter(
        d_pc_cg20, data['cg20_r50'][msk3], s=s3[msk3],
        alpha=.5, marker='^', ec='w', lw=.5)
    ax2.axhline(20, ls=':', lw=2, c='k')

    # Fill
    x = np.linspace(xmin, xmax, 100)
    y2 = np.ones(len(x)) * 20
    y1 = np.ones(len(x)) * 1000
    ax2.fill_between(x, y1, y2, color='grey', alpha=0.075)

    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(.2, 70)
    ax2.set_xlabel("Distance [pc]", fontsize=14)
    ax2.set_ylabel(r"$r_{50}$ [pc]", fontsize=14)

    plt.subplots_adjust(hspace=.025)
    # plt.show()
    plt.savefig("pms_plx.png", dpi=300)


if __name__ == '__main__':
    main()

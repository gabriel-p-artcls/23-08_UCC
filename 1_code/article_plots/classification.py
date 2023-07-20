
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import scienceplots
plt.style.use('science')

date = "0712"


def main():
    """
    """
    df = pd.read_csv("classif_data_" + date + ".csv")
    names = df['name'].values

    # df['CST'] = np.clip(df['CST'], a_min=0, a_max=40) / 40
    # plt.scatter(df['CST'], df['C_dens'], alpha=.5)
    # plt.xlabel("CST")
    # plt.ylabel("C_dens")
    # plt.show()

    ucc_plx = df['ucc_plx'].values
    msk_plx = ucc_plx <= 0
    ucc_plx[msk_plx] = np.nan
    ucc_plx = np.clip(ucc_plx, a_min=0.01, a_max=np.inf)
    d_cut_pc = 1491.6
    d_pc = 1000/ucc_plx

    # plt.scatter(df['ucc_pm_disp']/df['ucc_plx'], df['r_50_pc'])
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.show()

    pm_y = 5 * 1.414 / (d_pc/1000 * 4.74)
    dist_1 = df['ucc_pm_disp'] - pm_y
    dist_2 = df['ucc_pm_disp'] - 1.

    dist = np.array(dist_1)
    msk_plx = 1000/ucc_plx > d_cut_pc
    dist[msk_plx] = dist_2[msk_plx]

    vmin, vmax = -1, .1

    dist = np.clip(dist, a_min=vmin, a_max=vmax)

    fig, ax = plt.subplots(1, figsize=(6, 6))

    size = np.log10(1 + df['N_50'])**5
    # size = np.log10(df['r_50_pc'] * 10)**5

    msk1 = df['r_50_pc'] < 20

    im = ax.scatter(
        df['C_lkl'][msk1], df['C_dens'][msk1], c=-dist[msk1], alpha=.7, cmap='plasma',
        s=size[msk1], lw=.1, ec='k', vmin=-vmax, vmax=-vmin)

    ax.scatter(
        df['C_lkl'][~msk1], df['C_dens'][~msk1], c=-dist[~msk1], alpha=.7, cmap='plasma',
        s=size[~msk1], lw=.1, ec='k', marker='s', vmin=-vmax, vmax=-vmin)

    plt.xlabel(r"$C_{phot}$")
    plt.ylabel(r"$C_{dens}$")
    plt.xlim(0.01, 1.04)
    plt.ylim(0.01, 1.04)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad='5%')
    cbar = fig.colorbar(im, cax=cax)

    plt.savefig("classif.png", dpi=300)


if __name__ == '__main__':
    main()

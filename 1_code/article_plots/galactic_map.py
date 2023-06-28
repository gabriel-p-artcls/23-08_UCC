
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
# import astropy.coordinates as coord
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np


def main():
    """
    Gridspec idea: http://www.sc.eso.org/~bdias/pycoffee/codes/20160407/
                   gridspec_demo.html
    """

    # gc_frame = coord.Galactocentric()

    df = pd.read_csv("/home/gabriel/Github/UCC/add_New_DB/UCC_cat_20230619_in.csv")

    d_pc = 1000 / df['plx'].values
    d_pc = np.clip(d_pc, 1, 20000)

    sizes = []
    for i, cl in df.iterrows():
        sizes.append(len(cl['DB'].split(';')))
    sizes = np.array(sizes)

    # Galactic coordinates.
    eq = SkyCoord(ra=df['RA_ICRS']*u.degree, dec=df['DE_ICRS']*u.degree, frame='icrs')
    lb = eq.transform_to('galactic')
    lon = lb.l.wrap_at(180 * u.deg).radian * u.radian
    lat = lb.b.radian * u.radian
    # coords = SkyCoord(l=lon, b=lat, distance=d_pc*u.pc, frame='galactic')
    # # Galactocentric coordinates.
    # c_glct = coords.transform_to(gc_frame)
    # x_pc, y_pc, z_pc = c_glct.x, c_glct.y, c_glct.z
    # x_kpc, y_kpc, z_kpc = x_pc.value/1000, y_pc.value/1000, z_pc.value/1000

    mnan = np.isnan(d_pc)
    m1 = abs(d_pc) <= 1500
    m2 = (1500 < abs(d_pc)) & (abs(d_pc) <= 3000)
    m3 = abs(d_pc) > 3000

    # invert lon
    i_lon = -lon

    fig = plt.figure(figsize=(25, 25))
    gs = gridspec.GridSpec(6, 6)
    #     6, 6, height_ratios=[1, 1, 1, 0.04, .7, .7],
    #     width_ratios=[1, 1, 1, 1, 1, 1])
    # gs.update(hspace=0.03, wspace=.3)

    def set_ticks(ax):
        ax.set_xticklabels([])
        ax.text(.495, .04, r'0{\textdegree}', fontsize=10,
                transform=ax.transAxes)
        ax.text(.41, .02, r'90{\textdegree}', fontsize=10,
                transform=ax.transAxes, rotation=20)
        ax.text(.55, .01, r'270{\textdegree}', fontsize=10,
                transform=ax.transAxes, rotation=-20)
        plt.yticks(fontsize=15)
        ax.grid(True)

    ax = plt.subplot(gs[0:1, 0:2], projection="aitoff")
    set_ticks(ax)
    # ax.invert_xaxis() # does not work
    ax.scatter(i_lon[mnan], lat[mnan], s=np.exp(sizes[mnan]), alpha=.7, facecolor='none',
               ec='grey', zorder=5)
    plt.title(f'd=nan, N={mnan.sum()}', fontsize=15, x=.55)

    ax = plt.subplot(gs[0:1, 2:4], projection="aitoff")
    set_ticks(ax)
    ax.scatter(i_lon[m1], lat[m1], s=np.exp(sizes[m1]), alpha=.7, facecolor='none',
               ec='green', zorder=5)
    plt.title(r'$d\leq1.5$ [kpc], ' + f'N={m1.sum()}', fontsize=15, x=.59)

    ax = plt.subplot(gs[1:2, 0:2], projection="aitoff")
    set_ticks(ax)
    ax.scatter(i_lon[m2], lat[m2], s=np.exp(sizes[m2]), alpha=.7, facecolor='none',
               ec='blue', zorder=5)
    plt.title(r'$1.5<d\leq3$ [kpc], ' + f'N={m2.sum()}', fontsize=15, x=.62)

    ax = plt.subplot(gs[1:2, 2:4], projection="aitoff")
    set_ticks(ax)
    ax.scatter(i_lon[m3], lat[m3], s=np.exp(sizes[m3]), alpha=.7, facecolor='none',
               ec='orange', zorder=5)
    plt.title(r'$d>3$ [kpc], ' + f'N={m3.sum()}', fontsize=15, x=.59)

    # # Sun's coords according to the Galactocentric frame.
    # x_sun, z_sun = gc_frame.galcen_distance, gc_frame.z_sun
    # s_xys = SkyCoord(
    #     -x_sun, 0., z_sun, unit='kpc', representation_type='cartesian')

    # ax = plt.subplot(gs[4:6, 0:2])
    # plt.xlabel("X [kpc]", fontsize=40)
    # plt.ylabel("Y [kpc]", fontsize=40)
    # plt.scatter(x_kpc, y_kpc, alpha=.25)
    # plt.scatter(s_xys.x, s_xys.y, c='yellow', s=50, edgecolor='k')
    # plt.scatter(0., 0., c='k', marker='x', s=70)
    # # Plot spiral arms
    # plt.xlim(-15, 5)
    # plt.ylim(-10, 10)
    # plt.legend(fontsize=10)

    # ax = plt.subplot(gs[4:6, 2:4])
    # plt.xlabel("X [kpc]", fontsize=40)
    # plt.ylabel("Z [kpc]", fontsize=40)
    # plt.scatter(x_kpc, z_kpc, alpha=.25)
    # plt.scatter(s_xys.x, s_xys.z, c='yellow', s=50, edgecolor='k')
    # plt.scatter(0., 0., c='k', marker='x', s=70)
    # plt.xlim(-15, 5)
    # plt.ylim(-3, 4)

    # ax = plt.subplot(gs[4:6, 4:6])
    # plt.xlabel("Y [kpc]", fontsize=40)
    # plt.ylabel("Z [kpc]", fontsize=40)
    # plt.scatter(y_kpc, z_kpc, alpha=.25)
    # plt.scatter(s_xys.y, s_xys.z, c='yellow', s=50, edgecolor='k')
    # plt.xlim(-10, 10)
    # plt.ylim(-3, 4)

    # fig.tight_layout()
    fig.savefig('galactic_map.png', dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    import scienceplots
    plt.style.use('science')
    main()

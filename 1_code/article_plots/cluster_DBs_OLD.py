
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scienceplots
plt.style.use('science')


final_dbs_path = "/home/gabriel/Github/UCC/add_New_DB/UCC_cat_20230702_out.csv"
hunt23_membs_path = "/media/gabriel/rest/Dropbox_nosync/Papers/2023/GDR3_members/0_data/hunt23_members.parquet"
print("Reading HUNT23 members...")
hunt23_membs = pd.read_parquet(hunt23_membs_path)
print("Reading CANTAT20 members...")
cg20_path = "/media/gabriel/rest/Dropbox_nosync/Papers/2023/GDR3_members/0_data/CG_2020_members.csv.gz"
data_CG = pd.read_csv(cg20_path)
print("Reading fastMP output file...")
fastMP_db = pd.read_csv(final_dbs_path)
fnames = [_.split(';') for _ in fastMP_db['fnames']]

clpath = "/media/gabriel/rest/out/"


def main():
    cl_name = "NGC_3114"
    make_plot(cl_name, 900)


def make_plot(cl_name, Nmembs_manual):
    """
    """
    fastMP_name = cl_name.lower().replace('_', '').replace(
        ' ', '').replace('-', '').replace('.', '').replace('+', 'p')

    # if cl_name.startswith('VDBH_'):
    #     hunt23_name = 'BH_' + cl_name.split('_')[1]
    # if hunt23_name.startswith('VDB_'):
    #     hunt23_name = 'vdBergh_' + cl_name.split('_')[1]

    # Read HUNT23 members data
    msk1 = hunt23_membs['Name'] == cl_name
    hunt23_cl = hunt23_membs[msk1]

    # Read CANTAT20 data
    msk = data_CG['Cluster'] == cl_name
    cg20_cl = data_CG[msk]
    cg20_cl = cg20_cl.rename(columns={'pmRA*': 'pmRA'})

    # Read fastMP output catalogue
    i = np.nan
    for i, fname in enumerate(fnames):
        if fastMP_name in fname:
            break
    if np.isnan(i):
        print("Cluster not found")
        return

    row = fastMP_db.iloc[i]
    fname0 = row['fnames'].split(';')[0]
    # Search for the cluster file
    file_name = clpath + row['quad'] + '/datafiles/' + fname0 + '.parquet'
    try:
        cl_d = pd.read_parquet(file_name)
    except:
        print(f"Not found: {file_name}")

    probs = cl_d['probs'].values

    # Select the 'Nmembs_manual' stars with the largest probabilities
    idx = np.argsort(probs)[::-1][:Nmembs_manual]
    msk = np.full(len(probs), False)
    msk[idx] = True
    cl_manual = cl_d[msk]

    msk = probs > .5
    Nmembs = int(row['N_membs'])
    if msk.sum() < Nmembs:
        # Select the 'N_membs' stars with the largest probabilities
        idx = np.argsort(probs)[::-1][:Nmembs]
        msk = np.full(len(probs), False)
        msk[idx] = True
    fastMP_cl = cl_d[msk]

    l_min = min(min(hunt23_cl['GLON']), min(cg20_cl['GLON']), min(fastMP_cl['GLON']), min(cl_manual['GLON']))
    l_max = max(max(hunt23_cl['GLON']), max(cg20_cl['GLON']), max(fastMP_cl['GLON']), max(cl_manual['GLON']))
    b_min = min(min(hunt23_cl['GLAT']), min(cg20_cl['GLAT']), min(fastMP_cl['GLAT']), min(cl_manual['GLAT']))
    b_max = max(max(hunt23_cl['GLAT']), max(cg20_cl['GLAT']), max(fastMP_cl['GLAT']), max(cl_manual['GLAT']))
    pmRA_min = min(min(hunt23_cl['pmRA']), min(cg20_cl['pmRA']), min(fastMP_cl['pmRA']), min(cl_manual['pmRA']))
    pmRA_max = max(max(hunt23_cl['pmRA']), max(cg20_cl['pmRA']), max(fastMP_cl['pmRA']), max(cl_manual['pmRA']))
    pmDE_min = min(min(hunt23_cl['pmDE']), min(cg20_cl['pmDE']), min(fastMP_cl['pmDE']), min(cl_manual['pmDE']))
    pmDE_max = max(max(hunt23_cl['pmDE']), max(cg20_cl['pmDE']), max(fastMP_cl['pmDE']), max(cl_manual['pmDE']))
    BPRP_min = min(min(hunt23_cl['BP-RP']), min(cg20_cl['BP-RP']), min(fastMP_cl['BP-RP']), min(cl_manual['BP-RP']))
    BPRP_min -= .2
    BPRP_max = max(max(hunt23_cl['BP-RP']), max(cg20_cl['BP-RP']), max(fastMP_cl['BP-RP']), max(cl_manual['BP-RP']))
    G_min = min(min(hunt23_cl['Gmag']), min(cg20_cl['Gmag']), min(fastMP_cl['Gmag']), min(cl_manual['Gmag']))
    G_min -= 0.5
    G_max = max(max(hunt23_cl['Gmag']), max(cg20_cl['Gmag']), max(fastMP_cl['Gmag']), max(cl_manual['Gmag']))
    G_max += 0.5
    plx_min = min(min(hunt23_cl['Plx']), min(cg20_cl['Plx']), min(fastMP_cl['Plx']), min(cl_manual['Plx']))
    plx_max = max(max(hunt23_cl['Plx']), max(cg20_cl['Plx']), max(fastMP_cl['Plx']), max(cl_manual['Plx']))

    plt.figure(figsize=(20, 11))
    plot(1, 'green', hunt23_cl, 'HUNT23', l_min, l_max, b_min, b_max, pmRA_min, pmRA_max, pmDE_min, pmDE_max, BPRP_min, BPRP_max, G_min, G_max, plx_min, plx_max)
    plot(5, 'orange', cg20_cl, 'CANTAT20', l_min, l_max, b_min, b_max, pmRA_min, pmRA_max, pmDE_min, pmDE_max, BPRP_min, BPRP_max, G_min, G_max, plx_min, plx_max)
    plot(9, 'blue', fastMP_cl, 'fastMP', l_min, l_max, b_min, b_max, pmRA_min, pmRA_max, pmDE_min, pmDE_max, BPRP_min, BPRP_max, G_min, G_max, plx_min, plx_max)
    plot(13, 'purple', cl_manual, 'fastMP (sup)', l_min, l_max, b_min, b_max, pmRA_min, pmRA_max, pmDE_min, pmDE_max, BPRP_min, BPRP_max, G_min, G_max, plx_min, plx_max)
    plt.subplots_adjust(hspace=.025)
    plt.savefig(fastMP_name + ".png", dpi=300)


def plot(q, color, clust, lbl, l_min, l_max, b_min, b_max, pmRA_min, pmRA_max,
    pmDE_min, pmDE_max, BPRP_min, BPRP_max, G_min, G_max, plx_min, plx_max):
    """
    """
    ax = plt.subplot(4, 4, q)
    plt.scatter(clust['GLON'], clust['GLAT'], ec='w', lw=.5, alpha=.5, c=color)
    if q >= 13:
        plt.xlabel("lon [deg]")
    else:
        ax.set_xticklabels([])
    plt.ylabel("lat [deg]")
    plt.xlim(l_min, l_max)
    plt.ylim(b_min, b_max)

    ax = plt.subplot(4, 4, q + 1)
    plt.scatter(clust['pmRA'], clust['pmDE'], ec='w', lw=.5, alpha=.5, c=color)
    if q >= 13:
        plt.xlabel("pmRA [mas/yr]")
    else:
        ax.set_xticklabels([])
    plt.ylabel("pmDE [mas/yr]")
    plt.xlim(pmRA_min, pmRA_max)
    plt.ylim(pmDE_min, pmDE_max)

    ax = plt.subplot(4, 4, q + 2)
    plt.hist(clust['Plx'], alpha=.5, color=color, label=lbl) #, density=True
    plt.legend()
    if q >= 13:
        plt.xlabel("Plx [mas]")
    else:
        ax.set_xticklabels([])
    plt.ylabel("N")
    plt.xlim(plx_min, plx_max)

    ax = plt.subplot(4, 4, q + 3)
    plt.scatter(clust['BP-RP'], clust['Gmag'], ec='w', lw=.5, alpha=.5, c=color)
    # plt.gca().invert_yaxis()
    if q >= 13:
        plt.xlabel("BP-RP [mag]")
    else:
        ax.set_xticklabels([])
    plt.ylabel("G [mag]")
    plt.xlim(BPRP_min, BPRP_max)
    plt.ylim(G_max, G_min)


if __name__ == '__main__':
    main()


import csv
import numpy as np
import pandas as pd

date = "0712"
clpath = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/"
final_dbs_path = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/UCC_cat_20230702_out.csv"

cantat20_DB_path = "/home/gabriel/Github/UCC/add_New_DB/databases/CANTAT20.csv"
cantat20_membs_path = "/media/gabriel/rest/Dropbox_nosync/Papers/2023/GDR3_members/0_data/CG_2020_members.csv.gz"
hunt23_DB_path = "/home/gabriel/Github/UCC/add_New_DB/databases/HUNT23.csv"
hunt23_membs_path = "/media/gabriel/rest/Dropbox_nosync/Papers/2023/GDR3_members/0_data/hunt23_members.parquet"


def main():
    """
    """
    final_db = pd.read_csv(final_dbs_path)

    print("Reading CANTAT20 members...")
    cantat20_membs = pd.read_csv(cantat20_membs_path)
    cantat20_names = cantat20_membs['Cluster']
    cantat20_DB = pd.read_csv(cantat20_DB_path)

    print("Reading HUNT23 members...")
    hunt23_membs = pd.read_parquet(hunt23_membs_path)
    hunt23_names = hunt23_membs['Name']
    hunt23_DB = pd.read_csv(hunt23_DB_path)

    print("Processing clusters...")
    data = {
        'name': [], 'N_ucc': [], 'ucc_plx': [], 'ucc_pm_disp': [],
        'N_h23': [], 'h23_plx': [], 'h23_pm_disp': [], 'h23_r50': [],
        'CMDCl50': [],
        'N_cg20': [], 'cg20_plx': [], 'cg20_pm_disp': [], 'cg20_r50': [],
        'CST': [], 'C_lkl': [], 'C_dens': [], 'N_50': [], 'r_50_pc': []
    }

    for i, row in final_db.iterrows():

        data['name'].append(row['fnames'].split(';')[0])

        N_ucc, ucc_plx, ucc_pm_disp, C1, C2, N_50, r_50, r_50_pc = [np.nan] * 8
        N_ucc, ucc_plx, ucc_pm_disp, C1, C2, N_50, r_50, r_50_pc =\
            fastMP_vals(row)

        data['N_ucc'].append(N_ucc)
        data['ucc_plx'].append(ucc_plx)
        data['ucc_pm_disp'].append(ucc_pm_disp)
        data['C_lkl'].append(C1)
        data['C_dens'].append(C2)
        data['N_50'].append(N_50)
        data['r_50_pc'].append(r_50_pc)

        N_cg20, cg20_plx, cg20_pm_disp, cg20_r50 = [np.nan] * 4
        if 'CANTAT20' in row['DB']:
            N_cg20, cg20_plx, cg20_pm_disp, cg20_r50 = CG20_vals(
                row, cantat20_DB, cantat20_names, cantat20_membs)

        N_h23, h23_plx, h23_pm_disp, CST, h23_r50, CMDCl50 = [np.nan] * 6
        if 'HUNT23' in row['DB']:
            N_h23, h23_plx, h23_pm_disp, CST, h23_r50, CMDCl50 = H23_vals(
                row, hunt23_DB, hunt23_names, hunt23_membs)

        data['N_cg20'].append(N_cg20)
        data['cg20_plx'].append(cg20_plx)
        data['cg20_pm_disp'].append(cg20_pm_disp)
        data['cg20_r50'].append(cg20_r50)
        data['N_h23'].append(N_h23)
        data['h23_plx'].append(h23_plx)
        data['h23_pm_disp'].append(h23_pm_disp)
        data['h23_r50'].append(h23_r50)
        data['CMDCl50'].append(CMDCl50)
        data['CST'].append(CST)

        line = []
        for k, v in data.items():
            line.append(v[i])
        print(line)

    df = pd.DataFrame(data)
    df.to_csv("classif_data_" + date + ".csv", na_rep='nan', index=False,
              quoting=csv.QUOTE_NONNUMERIC)


def fastMP_vals(row, N_membs_min=25, plx_offset=0.029):
    """
    """
    # Read fastMP data for this cluster
    fname0 = row['fnames'].split(';')[0]
    Qfold = row['quad']
    try:
        cl_d = pd.read_parquet(
            clpath + Qfold + '/datafiles/' + fname0 + '.parquet')
        probs = cl_d['probs'].values
        msk = probs > 0.5
        Nmembs = N_membs_min
        if msk.sum() < Nmembs:
            # Select the 'N_membs' stars with the largest probabilities
            idx = np.argsort(probs)[::-1][:Nmembs]
            msk = np.full(len(probs), False)
            msk[idx] = True
        cl_d = cl_d[msk]            
    except:
        print(f"Could not find {Qfold}/{fname0}.parquet")
        return np.nan, np.nan
    N_membs = len(cl_d)

    # Estimate total PM dispersion
    pmRA = cl_d['pmRA'].values
    pmDE = cl_d['pmDE'].values
    pm_tot_disp = np.sqrt(np.std(pmRA)**2 + np.std(pmDE)**2)

    # N_membs, pm_tot_disp = np.nan, 0

    C1, C2, N_50, r_50 = row['C1'], row['C2'], row['N_50'], row['r_50']

    # Parallax as median
    # plx = np.median(cl_d['Plx'].values)
    plx = row['plx_m'] + plx_offset
    if plx <= 0:
        plx = np.nan

    d_pc = 1000/plx
    r_50_rad = np.deg2rad(r_50/60)  # Arcmin to degrees to rad
    r_50_pc = round(d_pc * np.tan(r_50_rad), 2)

    return N_membs, round(plx, 3), round(pm_tot_disp, 3), C1, C2, N_50, r_50,\
        r_50_pc


def CG20_vals(row, cantat20_DB, cantat20_names, cantat20_membs):
    """
    """
    # Identify cluster in CANTAT20 DB
    dbs_idx = row['DB_i'].split(';')
    h23_j = row['DB'].split(';').index('CANTAT20')
    cantat20_idx = int(dbs_idx[h23_j])
    cantat20_name = cantat20_DB['Cluster'][cantat20_idx].strip()

    # Identify the members for this cluster in CANTAT20
    if cantat20_name.startswith('FoF_'):
        cantat20_name = cantat20_name.replace('FoF', 'LP')
    if cantat20_name.startswith('VDBH_'):
        cantat20_name = 'BH_' + cantat20_name.split(',')[0].split('_')[1]
    if cantat20_name.startswith('VDB_'):
        cantat20_name = 'vdBergh_' + cantat20_name.split(',')[0].split('_')[1]
    msk = cantat20_names == cantat20_name
    N_cg20 = msk.sum()
    if N_cg20 == 0:
        print(f"Could not find cluster {cantat20_name} in CANTAT20")
        return np.nan, np.nan, np.nan, np.nan
    # Estimate total PM dispersion
    pmRA = cantat20_membs['pmRA*'][msk].values
    pmDE = cantat20_membs['pmDE'][msk].values
    pm_tot_disp = np.sqrt(np.std(pmRA)**2 + np.std(pmDE)**2)
    # Parallax as median
    plx = np.median(cantat20_membs['Plx'][msk].values)

    d_pc = 1000/cantat20_DB['plx'][cantat20_idx]
    r_50 = cantat20_DB['r50'][cantat20_idx]
    r_50_rad = np.deg2rad(r_50)  # Degrees to rad
    cg20_r50 = round(d_pc * np.tan(r_50_rad), 2)

    # N_cg20, pm_tot_disp, plx = 0, 0, 0

    return N_cg20, round(plx, 3), round(pm_tot_disp, 3), cg20_r50


def H23_vals(row, hunt23_DB, hunt23_names, hunt23_membs):
    """
    """
    # Identify cluster in HUNT23 DB
    dbs_idx = row['DB_i'].split(';')
    h23_j = row['DB'].split(';').index('HUNT23')
    hunt23_idx = int(dbs_idx[h23_j])
    hunt23_name = hunt23_DB['Name'][hunt23_idx].strip()

    # Identify the members for this cluster in HUNT23
    if hunt23_name.startswith('VDBH_'):
        hunt23_name = 'BH_' + hunt23_name.split(',')[0].split('_')[1]
    if hunt23_name.startswith('VDB_'):
        hunt23_name = 'vdBergh_' + hunt23_name.split(',')[0].split('_')[1]
    msk = hunt23_names == hunt23_name
    N_h23 = msk.sum()
    if msk.sum() == 0:
        print(f"Could not find cluster {hunt23_name} in HUNT23")
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    # Estimate total PM dispersion
    pmRA = hunt23_membs['pmRA'][msk].values
    pmDE = hunt23_membs['pmDE'][msk].values
    pm_tot_disp = np.sqrt(np.std(pmRA)**2 + np.std(pmDE)**2)
    # Parallax as median
    plx = np.median(hunt23_membs['Plx'][msk].values)

    CST = hunt23_DB['CST'][hunt23_idx]
    h23_r50 = hunt23_DB['r50pc'][hunt23_idx]
    CMDCl50 = hunt23_DB['CMDCl50'][hunt23_idx]

    # N_h23, plx, pm_tot_disp, CST = 0, 0, 0, 0

    return N_h23, round(plx, 3), round(pm_tot_disp, 3), CST, h23_r50, CMDCl50


if __name__ == '__main__':
    main()


from os import listdir
import pandas as pd
import numpy as np


# CG2020
path_CG = "../0_data/cantat_gaudin_et_al_2020/cg2020.csv"
data_CG = pd.read_csv(path_CG)
names_CG = list(np.array([_.lower() for _ in data_CG['Name']]))

folders = ('S1', 'S2', 'M1', 'L1', 'L2', 'L3')
# folders = ('L2', 'L2', 'L3')

for folder in folders:
    path = f"../0_data/{folder}/"
    files = listdir(path)

    for file in files:
        # if file[:-4] != 'Alessi_2':
        #     continue

        df = pd.read_csv(path + file)
        rlon, rlat = np.ptp(df['GLON']), np.ptp(df['GLAT'])

        fname = file[:-4]
        idx = names_CG.index(fname.lower())
        r50 = data_CG['r50'][idx]

        box_s = max(rlon, rlat)

        plx_d = float(data_CG['plx'][idx])
        if plx_d > 10:
            plx_p = .5
        elif plx_d <= 10 and plx_d > 5:
            plx_p = .3
        else:
            plx_p = .2
        plx_min = plx_d - plx_d * plx_p
        if df['Plx'].min() < plx_min - plx_min * 0.05:
            plx_method = 'full'
        else:
            plx_method = 'plx_min'

        print("{}, {}, {:.1f}, {:.1f}, {}".format(
            fname, folder, box_s, box_s / r50, plx_method))


import numpy as np
import pandas as pd


"""
Count the number of clusters with no data for a given parameter
"""
df = pd.read_csv("final_DB.csv")

no_plx, no_pmRA, no_pmDE, no_plxpm, just_hao, msk = 0, 0, 0, 0, 0, []
for i, cl in df.iterrows():
    if np.isnan(cl['plx']):
        no_plx += 1
    if np.isnan(cl['pmRA']):
        no_pmRA += 1
    if np.isnan(cl['pmDE']):
        no_pmDE += 1
    if np.isnan(cl['plx']) and np.isnan(cl['pmRA']) and np.isnan(cl['pmDE']):
        no_plxpm += 1
        if cl['ID'].startswith('MWSC'):
            print(cl['DB'], cl['ID'])

print(no_plx, no_pmRA, no_pmDE, no_plxpm, just_hao)
# msk = np.array(msk)

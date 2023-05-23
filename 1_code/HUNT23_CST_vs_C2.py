
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


fastMP_out_path = "/home/gabriel/Github/UCC/add_New_DB/UCC_cat_20230520.csv"
HUNT23_path = "/home/gabriel/Github/UCC/add_New_DB/databases/HUNT23.csv"

db_fastMP = pd.read_csv(fastMP_out_path)
db_hunt23 = pd.read_csv(HUNT23_path)

fnames = [_.split(';') for _ in db_fastMP['fnames']]

CST_C2 = []
for i, row in db_hunt23.iterrows():
    fname_H23 = row['name'].lower().replace('_', '')
    # Read gastMP output catalogue
    j = np.nan
    for j, fname in enumerate(fnames):
        if fname_H23 in fname:
            break
    if np.isnan(j):
        print("Cluster not found", fname_H23)
        continue

    CST = row['cst']
    C2 = db_fastMP['C2'][j]
    CST_C2.append([CST, C2])
    if CST > 90:
        print(fname_H23, fnames[j], CST, C2)

CST_C2 = np.array(CST_C2).T

plt.scatter(CST_C2[0], CST_C2[1])
plt.show()

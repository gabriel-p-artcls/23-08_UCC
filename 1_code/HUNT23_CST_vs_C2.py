
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


fastMP_out_path = "/home/gabriel/Github/UCC/add_New_DB/UCC_cat_20230620_out.csv"
HUNT23_path = "/home/gabriel/Github/UCC/add_New_DB/databases/HUNT23.csv"

db_fastMP = pd.read_csv(fastMP_out_path)
db_hunt23 = pd.read_csv(HUNT23_path)

fnames = [_.split(';') for _ in db_fastMP['fnames']]

CST_C12 = []
for i, row in db_hunt23.iterrows():
    hunt23_name = row['Name']
    if ',' in hunt23_name:
        hunt23_name = hunt23_name.split(',')[0]
    fname_H23 = hunt23_name.lower().replace('_', '').replace(
        ' ', '').replace('-', '').replace('.', '').replace('+', 'p')

    # Read fastMP output catalogue
    j_old = np.nan
    for j, fname in enumerate(fnames):
        if fname_H23 in fname:
            j_old = j
            break
    if np.isnan(j_old):
        print("Cluster not found", hunt23_name)
        continue

    CST = row['CST']
    C1, C2 = db_fastMP['C1'][j], db_fastMP['C2'][j]
    CST_C12.append([CST, C1, C2, db_fastMP['N_membs'][j]])

    # if abs(C3-C4)>0.5:
    #     print(i, fname_H23, j, fnames[j], C3, C4, CST)

CST_C12 = np.array(CST_C12).T

# plt.subplot(131)
plt.scatter(CST_C12[1], CST_C12[2], s=CST_C12[-1]/10, alpha=.5)
# plt.xlabel("C1")
# plt.ylabel("C2")
# plt.subplot(132)
# plt.scatter(CST_C12[0], CST_C12[1])
# plt.xlabel("CST")
# plt.ylabel("C1")
# plt.subplot(133)
# plt.scatter(CST_C12[0], CST_C12[2])
# plt.xlabel("CST")
# plt.ylabel("C2")
plt.show()

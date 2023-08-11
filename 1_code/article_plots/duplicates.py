
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from string import digits
import scienceplots
plt.style.use('science')

"""
Count duplicates in selected IDs
"""

path_to_UCC_Db = "/home/gabriel/Github/UCC/updt_UCC/UCC_cat_230702.csv"

db_years = {
    'KHARCHENKO12': 1200, 'LOKTIN17': 1700, 'CASTRO18': 1800,
    'BICA19': 1901, 'CASTRO19': 1907, 'SIM19': 1910, 'LIUPANG19': 1912,
    'FERREIRA19': 1903,
    'CANTAT20': 2008, 'CASTRO20': 2003, 'FERREIRA20': 2008, 'HAO20': 2003,
    'DIAS21': 2106, 'CASADO21': 2106, 'FERREIRA21': 2103, 'HUNT21': 2102,
    'JAEHNIG21': 2112, 'SANTOS21': 2111, 'HE21': 2105,
    'CASTRO22': 2205, 'TARRICQ22': 2203, 'LI22': 2203, 'HE22': 2205,
    'HE22_1': 2209, 'HAO22': 2204,
    'HE23': 2301, 'HUNT23': 2305, 'QIN23': 2303, 'LI23': 2303,
    'CHI23_2': 2303, 'CHI23': 2306, 'CHI23_3': 2306
}

# Total number of entries for these IDs
ids_total = {
    'ocsn': 0, 'cwwdl': 0, 'hsc': 0, 'cwnu': 0, 'oc': 0, 'fof': 0,
    'lisciii': 0, 'hxhwl': 0, 'lisc': 0, 'hxwhb': 0, 'cwwl': 0
}
ids_dup = {
    'ocsn': 0, 'cwwdl': 0, 'hsc': 0, 'cwnu': 0, 'oc': 0, 'fof': 0,
    'lisciii': 0, 'hxhwl': 0, 'lisc': 0, 'hxwhb': 0, 'cwwl': 0  # 'xdocc': 0, 
}

df = pd.read_csv(path_to_UCC_Db)
fnames = [_.split(';')[0] for _ in df['fnames']]

# For each cluster in the DB
for i, row in df.iterrows():

    # Check that this is one of the OCs we want to check for duplicates
    cid_found = 'nan'
    fname = row['fnames'].split(';')[0]
    fname_digits = str.maketrans('', '', digits)
    fname_nodigs = fname.translate(fname_digits)
    for cid in ids_total.keys():
        if fname_nodigs == cid:
            cid_found = cid
    if cid_found == 'nan':
        continue
    ids_total[cid_found] += 1

    # This cluster has no duplicates
    if str(row['dups_fnames_m']) == 'nan':
        continue

    # Check to see if at least one probability is >0.5
    dup_j = []
    probs = row['dups_probs_m'].split(';')
    for j, prob in enumerate(probs):
        if float(prob) > 0.5:
            dup_j.append(j)  # Index of duplicate with P>0.5
    if len(dup_j) == 0:
        continue

    # Database(s) where this cluster is listed
    cl_dbs = row['DB'].split(';')
    cl_years = [db_years[_] for _ in cl_dbs]
    # Select the minimum year
    cl_year = min(cl_years)

    # Database(s) where each duplicate exists
    dup_dbs = []
    names_dups = row['dups_fnames_m'].split(';')
    for j in dup_j:
        dup_fname = names_dups[j]
        k = fnames.index(dup_fname)
        dup_dbs += df['DB'][k].split(';')
    dup_dbs = list(set(dup_dbs))
    dup_years = [db_years[_] for _ in dup_dbs]
    dup_year = min(dup_years)

    if dup_year < cl_year:
        # This cluster is a duplicate of a previously presented object
        ids_dup[cid_found] += 1

cl_percentage = {}
for cl_id, Nd in ids_dup.items():
    cl_percentage[cl_id] = Nd / ids_total[cl_id]

ids_sorted = dict(sorted(
    cl_percentage.items(), key=lambda x: x[1], reverse=True))
print(ids_sorted)

# Select first elements
cl_groups = {k: ids_sorted[k] for k in list(ids_sorted)[:5]}

xticks = np.arange(len(cl_groups))

heights = np.array(list(cl_groups.values())) * 100

# fig = plt.figure(figsize=(5, 5))
ax = plt.subplot()
ax.bar(xticks, heights, alpha=.5)
ax.set_ylabel(r"\% of probable duplicates")
ax.set_xticks(xticks)
xlabels = [_.upper() for _ in cl_groups.keys()]
xlabels = [_.replace('FOF', 'FoF') for _ in xlabels]
ax.set_xticklabels(xlabels, rotation=45)
ax.tick_params(axis='both', which='minor', bottom=False)
plt.savefig("duplicates.png", dpi=300)


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots
plt.style.use('science')


# dd = {'A':1, 'B':.75, 'C':.5, 'D':.25}
# for c1 in ('A', 'B', 'C', 'D'):
#     for c2 in ('A', 'B', 'C', 'D'):
#         print(c1+c2, dd[c1]+dd[c2])
# AA 2
# AB 1.75
# BA 1.75
# --
# AC 1.5
# CA 1.5
# BB 1.5
# --
# AD 1.25
# DA 1.25
# BC 1.25
# CB 1.25
# --
# BD 1.0
# DB 1.0
# CC 1.0
# --
# CD 0.75
# DC 0.75
# --
# DD 0.5

date = "0712"
clpath = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/"
final_dbs_path = "/media/gabriel/backup/gabriel/UCC/out_" + date + "/UCC_cat_20230702_out.csv"

final_db = pd.read_csv(final_dbs_path)

unique_c = list(set(final_db['C3']))

lst = {}
for c in unique_c:
    N = (final_db['C3']==c).sum()
    lst[c] = N
    print(c, N)

k_order = [
    'AA', 'AB', 'BA', 'AC', 'CA', 'BB', 'AD', 'DA', 'BC', 'CB', 'BD', 'DB',
    'CC', 'CD', 'DC', 'DD']

height = []
for c in k_order:
    height.append(lst[c])

x, y = k_order, height


def rescale(x):
    return (x - np.min(x)) / (np.max(x) - np.min(x))


def my_barchart(my_df, my_cmap):
    # rescale = lambda y: (y - np.min(y)) / (np.max(y) - np.min(y))

    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    ax.bar(my_df['days'], my_df['y'], color=my_cmap(rescale(my_df['y'])))
    return fig


my_cmap = plt.get_cmap("RdYlGn_r")
x1 = np.arange(0, 17)

plt.figure(figsize=(6, 3))

plt.bar(x, y, color=my_cmap(rescale(x1)))
plt.xticks(rotation=30)
plt.ylabel("N")

plt.savefig("classif_bar.png", dpi=300)

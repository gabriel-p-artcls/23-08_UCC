
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
Plot the results for the methods vs CG2020 tests
"""

path = "../2_pipeline/"
data = pd.read_csv(path + "methods_res.csv")
N_tot = len(list(set(data['name'])))

msk_py = data['method'] == 'py'

# methods = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
#            25, 26, 27, 28, 29, 30, 31, 32)
# methods = (1, 2, 3, 17, 18, 19)
# methods = range(1, 17)
# # (1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15, 8, 16)
# methods = range(17, 21)
# methods = (1, 17, 2, 18, 5, 19, 6, 20, 9, 21, 10, 22, 13, 23, 14, 24)
methods = (51, 59, 61, 63)

all_tests = (
    'perc_diff', 'perc_missed_-16', 'perc_missed_+16', 'delta_coords_std',
    'delta_pms_std', 'delta_plx_std')
for i, test in enumerate(all_tests):
    dd = data[test]
    sbpl = int("23" + str(i + 1))
    plt.subplot(sbpl)
    plt.title(test)
    box_data = [dd[msk_py]]
    for mi in methods:
        msk = data['method'] == 'M' + str(mi)
        msk_nan = np.isnan(dd[msk].values)
        if msk_nan.sum() > 0:
            print("nans in ", data['name'][msk][msk_nan].values[0], mi)
            ddd = dd[msk][~msk_nan]
        else:
            ddd = dd[msk]
        box_data.append(ddd)
    plt.boxplot(box_data, positions=range(0, len(box_data)),
                labels=['py'] + list(methods))
    plt.axhline(0, ls=':', c='k')

plt.show()

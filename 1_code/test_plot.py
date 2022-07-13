
from os import listdir
import pandas as pd
import matplotlib.pyplot as plt

# See which clusters are badly generates (ie: Alessi_20)

path = '../2_pipeline/pyUPMASK/303_no_flags/256_GMM/'
files = listdir(path)

for cl in files:
    print(cl)
    data = pd.read_csv(path + cl, sep=' ')
    msk = data['probs_final'] > 0.95
    cluster = data[msk]
    field = data[~msk]

    plt.figure()
    plt.subplot(221)
    plt.scatter(field['GLON'], field['GLAT'], alpha=.25)
    plt.scatter(cluster['GLON'], cluster['GLAT'], alpha=.5, zorder=3)
    plt.subplot(222)
    # plt.scatter(field['pmRA'], field['pmDE'], alpha=.25)
    plt.scatter(cluster['pmRA'], cluster['pmDE'], alpha=.5, zorder=3)
    plt.subplot(223)
    plt.hist(cluster['Plx'], 50, zorder=3)
    plt.subplot(224)
    plt.scatter(field['BP-RP'], field['Gmag'], alpha=.25)
    plt.scatter(cluster['BP-RP'], cluster['Gmag'], alpha=.5, zorder=3)
    plt.gca().invert_yaxis()
    plt.savefig(cl.replace('csv', 'png'), dpi=200)
    plt.close()

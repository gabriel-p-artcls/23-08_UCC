
from os import listdir
import pandas as pd
import numpy as np
from scipy.stats import median_abs_deviation as MAD
import matplotlib.pyplot as plt


# Path to clusters datafiles
path = "/home/gabriel/Descargas/out_1/"


def main():
    """
    """

    df = pd.read_csv("/home/gabriel/Descargas/table_1.csv")

    x = df['N_membs'].values
    sz = np.clip(10 / df['plx'].values, a_min=10, a_max=200)

    abcd_c = df['Class'].values
    class_2_numb = {'A': 1., 'B': 0.5, 'C': 0.25, 'D': 0.1}
    y = []
    for ccl in abcd_c:
        c1, c2, c3 = ccl
        y.append((class_2_numb[c1] + class_2_numb[c2] + class_2_numb[c3])/3)

    makePlot(x, y, df['ID'].values, sz)
    breakpoint()

    QCD = []
    for cluster in listdir(path):
        df = pd.read_csv(path + cluster)

        cxy = (np.median(df['GLON']), np.median(df['GLAT']))
        dxy = np.sqrt((df['GLON'] - cxy[0])**2 + (df['GLAT'] - cxy[1])**2)

        cpm = (np.median(df['pmRA']), np.median(df['pmDE']))
        dpm = np.sqrt((df['pmRA'] - cpm[0])**2 + (df['pmDE'] - cpm[1])**2)

        # Qxy = np.mean([quart_coeff_disp(df['GLON']), quart_coeff_disp(df['GLAT'])])
        # Qpm = np.mean([quart_coeff_disp(df['pmRA']), quart_coeff_disp(df['pmDE'])])

        # Qxy = quart_coeff_disp(dxy)
        # Qpm = quart_coeff_disp(dpm)

        Qxy = MAD(dxy) / np.median(dxy)
        Qpm = MAD(dpm) / np.median(dpm)

        QCD.append([Qxy, Qpm, 10/np.median(df['Plx']), len(df)])

        print("{}: {:.2f}, {:.2f} ({})".format(cluster, Qxy, Qpm, len(df)))

    QCD = np.array(QCD).T

    plt.subplot(121)
    plt.scatter(QCD[0], QCD[1], s=QCD[2], alpha=.5)
    plt.subplot(122)
    sz = np.clip(QCD[3], a_min=10, a_max=200)
    plt.scatter(QCD[0], QCD[1], s=sz, alpha=.5)
    plt.show()


def quart_coeff_disp(x):
    """
    """
    q25, q75 = np.percentile(x, 25), np.percentile(x, 75)
    return (q75 - q25) / (q75 + q25)


def makePlot(x, y, names, sz):
    """
    Source: https://stackoverflow.com/a/47166787/1391441
    """
    fig, ax = plt.subplots()

    sc = plt.scatter(x, y, s=sz, alpha=.5)

    annot = ax.annotate(
        "", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
        bbox=dict(boxstyle="round", fc="w"), arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
        annot.set_text(text)
        # annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)
    plt.show()


if __name__ == '__main__':
    main()

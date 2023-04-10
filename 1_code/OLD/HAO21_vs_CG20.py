
"""
Compare the clusters in the HAO21 database with those matched in the CG20
database.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


input_folder = '../0_data/databases/old/'


def main():
    """
    """
    hao21 = pd.read_csv(input_folder + 'HAO21.csv')
    cg20 = pd.read_csv(input_folder + 'CG20.csv')

    names_hao21 = [_.lower().replace('_', '').replace(' ', '').replace('-', '')
                   for _ in hao21['Cluster']]
    names_cg20 = [_.lower().replace('_', '').replace(' ', '').replace('-', '')
                  for _ in cg20['Name']]

    names, cl_match_hao, cl_match_cg = [], [], []
    for i, cl_a in enumerate(names_hao21):
        try:
            j = names_cg20.index(cl_a)
            # RA, DEC values are in good agreement so I removed them
            names.append([hao21['Cluster'][i], cg20['Name'][j]])
            cl_match_hao.append(list(hao21[['plx', 'pmRA*', 'pmDE']].values[i]))
            cl_match_cg.append(list(cg20[['plx', 'pmRA', 'pmDE']].values[j]))
        except:
            pass

    names = np.array(names)
    print("Matched: ", len(names))
    cl_match_hao = np.array(cl_match_hao).T
    cl_match_cg = np.array(cl_match_cg).T

    print("\nRel diff > 25%")
    for i, col in enumerate(('plx', 'pmRA', 'pmDE')):
        print(col, (abs(cl_match_hao[i] - cl_match_cg[i]) / cl_match_hao[i]
              > .25).sum())

    print("\nRel diff > 50%")
    for i, col in enumerate(('plx', 'pmRA', 'pmDE')):
        print(col, (abs(cl_match_hao[i] - cl_match_cg[i]) / cl_match_hao[i]
              > .5).sum())

    print("\nRel diff > 75%")
    for i, col in enumerate(('plx', 'pmRA', 'pmDE')):
        print(col, (abs(cl_match_hao[i] - cl_match_cg[i]) / cl_match_hao[i]
              > .75).sum())

    print("\nRel diff > 90%")
    for i, col in enumerate(('plx', 'pmRA', 'pmDE')):
        print(col, (abs(cl_match_hao[i] - cl_match_cg[i]) / cl_match_hao[i]
              > .90).sum())

    names = names.T[0]
    fig, ax = plt.subplots()
    makePlot(fig, ax, cl_match_hao[0], cl_match_cg[0], names)
    ax.set_xlabel("HAO21 plx")
    ax.set_ylabel("CG20 plx")
    plt.show()

    fig, ax = plt.subplots()
    makePlot(fig, ax, cl_match_hao[1], cl_match_cg[1], names)
    ax.set_xlabel("HAO21 pmRA")
    ax.set_ylabel("CG20 pmRA")
    plt.show()

    fig, ax = plt.subplots()
    makePlot(fig, ax, cl_match_hao[2], cl_match_cg[2], names)
    ax.set_xlabel("HAO21 pmDE")
    ax.set_ylabel("CG20 pmDE")
    plt.show()

    # plt.subplot(132)
    # plt.scatter(cl_match_hao[1], cl_match_cg[1])
    # plt.xlabel("HAO21 pmRA")
    # plt.ylabel("CG20 pmRA")
    # plt.subplot(133)
    # plt.scatter(cl_match_hao[2], cl_match_cg[2])
    # plt.xlabel("HAO21 pmDE")
    # plt.ylabel("CG20 pmDE")
    # plt.show()
    breakpoint()


def makePlot(fig, ax, x, y, names):
    """
    Source: https://stackoverflow.com/a/47166787/1391441
    """
    sc = plt.scatter(x, y, s=5)

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
    # plt.show()


if __name__ == '__main__':
    # plt.style.use('science')
    main()

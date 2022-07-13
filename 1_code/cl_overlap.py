
import numpy as np
import pandas as pd
from astropy.stats import sigma_clipped_stats
from scipy.spatial import Voronoi, ConvexHull, voronoi_plot_2d
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt


def main(vorplot=False):
    """
    Helper script. Check the overlap between clusters catalogued in CG2020
    """

    # Read file data will all the clusters
    data_all_cls = pd.read_csv("../0_data/cantat_gaudin_et_al_2020/cg2020.csv")
    # names = np.array(list(data_all_cls['Name']))

    x, y = data_all_cls['GLON'], data_all_cls['GLAT']
    coords = np.array([x, y]).T

    if vorplot:
        # Print and plot voronoi cells
        vor, vols = voronoi_volumes(coords)
        idxs = np.argsort(vols)
        for i in idxs:
            cl = data_all_cls['Name'][i]
            area_arcmin = np.sqrt(vols[i] * 3600)
            print(cl, x[i], y[i], area_arcmin)
        names = data_all_cls['Name']
        makePlot(x, y, names, vor)

    # Find the distances to all clusters, for all clusters
    dist = cdist(coords, coords)

    for i, xy_c in enumerate(coords):
        breakpoint()
        dd = dist[i]
        msk = (xy_c[0])

    # match_found = []
    # for i, d in enumerate(dist):
    #     msk1 = (d < 1) & (d > 0.)
    #     if msk1.sum() > 1:
    #         names_i = names[msk1]
    #         name_all = ''
    #         for name in names_i:
    #             if name != names[i]:
    #                 name_comb = names[i] + name
    #                 if name_comb not in match_found:
    #                     match_found.append(name_comb)
    #                     name_all += ' ' + name
    #         if name_all != '':
    #             print(names[i] + name_all)
    #     # breakpoint()


def makePlot(x, y, names, vor):
    """
    Source: https://stackoverflow.com/a/47166787/1391441
    """
    fig, ax = plt.subplots()

    voronoi_plot_2d(
        vor, show_vertices=False, line_colors='orange', line_width=2,
        line_alpha=0.6, point_size=0, ax=ax)

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
    plt.show()


def voronoi_volumes(points, Nsigma=3.):
    """
    For an D-dimensional dataset obtain its Voronoi diagram, and calculate
    the volume associated to each region. Unbounded regions are assigned the
    95th percentile volume. Outlier regions (with large volumes) are clipped to
    median volumes.
    """

    Ndim = points.shape[1]
    v = Voronoi(points)

    vol = np.zeros(v.npoints)
    for i, reg_num in enumerate(v.point_region):

        # Indices of the Voronoi vertices forming each Voronoi region.
        # -1 indicates vertex outside the Voronoi diagram.
        indices = v.regions[reg_num]

        # If the region is not bounded, assign a 'nan' volume.
        if -1 in indices:
            vol[i] = np.nan
        else:
            # # Clip vertexes outside of the frame's boundaries
            # vertxs = np.clip(v.vertices[indices], a_min=0., a_max=1.)

            # Coordinates of the Voronoi vertices.
            vertxs = v.vertices[indices]
            # Obtain volume for this region
            vol[i] = area_of_polygon(vertxs, Ndim)

    # For points with unbounded regions, assign the 95th percentile volume.
    vol[np.isnan(vol)] = np.nanpercentile(vol, 95)

    # Clip volumes of N-sigma outliers to the median volume of the dataset.
    mean, median, std = sigma_clipped_stats(vol)
    vol[vol > mean + Nsigma * std] = median

    return v, vol


def area_of_polygon(points, Ndim):
    """
    Calculates the area of an arbitrary polygon given its vertices

    Source: http://stackoverflow.com/a/4682656/1391441
    """
    if Ndim > 2:
        # N-dimensional approach (slower)
        p_area = ConvexHull(points).volume
    else:
        # For Ndim=2 this is faster than using ConvexHull()
        x, y = zip(*points)
        area = 0.0
        for i in range(-1, len(x) - 1):
            area += x[i] * (y[i + 1] - y[i - 1])
        p_area = abs(area) / 2.0

    return p_area


if __name__ == '__main__':
    main()

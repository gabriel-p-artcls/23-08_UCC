
import numpy as np
import pandas as pd
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import angular_separation
from scipy.spatial import Voronoi, ConvexHull, voronoi_plot_2d
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt


def main(vorplot=True):
    """
    Helper script. Plot the Voronoi volumes of the clusters and find the
    smallest radius that delimits the region with no other clusters inside.
    """
    # Read file data will all the clusters
    # data_all_cls = pd.read_csv("../2_pipeline/final_catalog.csv")

    df_new = pd.read_csv("../1_code/NEW_DBs.csv")
    df_old = pd.read_csv("../1_code/OLD_DBs.csv")
    data_all_cls = pd.concat([df_old, df_new])
    data_all_cls = data_all_cls.reset_index()

    x, y = data_all_cls['GLON'], data_all_cls['GLAT']
    pmRA, pmDE, plx = data_all_cls['pmRA'], data_all_cls['pmDE'], data_all_cls['plx']
    coords = np.array([x, y]).T

    # Find the distances to all clusters, for all clusters
    dist = cdist(coords, coords)
    msk = dist == 0.
    dist[msk] = np.inf
    dist_idx = np.argmin(dist, 0)

    cl_dists, clusts = [], []
    for i, j in enumerate(dist_idx):
        # d = round(dist[i][j] * 60, 3)
        d = round(angular_separation(x[i], y[i], x[j], y[j]) * 60, 3)
        pm_d = np.sqrt((pmRA[i]-pmRA[j])**2 + (pmDE[i]-pmDE[j])**2)
        plx_d = abs(plx[i] - plx[j])
        if d < 1 and pm_d < 1 and plx_d < .05:
            cl_dists.append(dist[i][j])
            clusts.append(
                [data_all_cls['DB'][i], data_all_cls['ID'][i],
                 data_all_cls['DB'][j], data_all_cls['ID'][j], d,
                 round(pm_d, 2), round(plx_d, 3), round(plx[i], 3)])
    idx = np.argsort(cl_dists)
    clusts = np.array(clusts)

    for cl in clusts[idx]:
        print(cl)
    breakpoint()

    # Print and plot voronoi cells
    vor, vols = voronoi_volumes(coords)
    idxs = np.argsort(vols)
    for i in idxs:
        cl = data_all_cls['ID'][i]
        area_arcmin = np.sqrt(vols[i] * 3600)
        # print(cl, x[i], y[i], round(area_arcmin, 1))
        print(cl, round(area_arcmin, 1))

    if vorplot:
        names = data_all_cls['ID']
        makePlot(x, y, names, vor)


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

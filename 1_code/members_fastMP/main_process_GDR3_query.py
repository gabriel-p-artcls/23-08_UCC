
import gzip
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np


# ******** IMPORTANT ********
#
# Clusters that wrap around the edges of the (ra, dec) coordinates are not
# still properly process; e.g.: Blanco 1
#
# ******** IMPORTANT ********


# # Gaia EDR3 zero points. Sigmas are already squared here.
Zp_G, sigma_ZG_2 = 25.6873668671, 0.00000759
Zp_BP, sigma_ZBP_2 = 25.3385422158, 0.000007785
Zp_RP, sigma_ZRP_2 = 24.7478955012, 0.00001428


def run(frames_path, fdata, c_ra, c_dec, box_s_eq, plx_min, max_mag):
    """
    box_s_eq: Size of box to query (in degrees)
    """
    # print("  ({:.3f}, {:.3f}); Box size: {:.2f}, Plx min: {:.2f}".format(
    #       c_ra, c_dec, box_s_eq, plx_min))

    data_in_files, xmin_cl, xmax_cl, ymin_cl, ymax_cl = findFrames(
        c_ra, c_dec, box_s_eq, fdata)

    all_frames = query(
        c_ra, c_dec, box_s_eq, frames_path, max_mag, data_in_files, xmin_cl,
        xmax_cl, ymin_cl, ymax_cl, plx_min)

    # print("Adding uncertainties")
    all_frames = uncertMags(all_frames)
    all_frames = all_frames.drop(columns=[
        'FG', 'e_FG', 'FBP', 'e_FBP', 'FRP', 'e_FRP'])

    return all_frames


def findFrames(c_ra, c_dec, box_s_eq, fdata):
    """
    """
    # # These are the points that determine the range of *all* the frames
    # l_min, l_max = fdata['l_min'], fdata['l_max']
    # b_min, b_max = fdata['b_min'], fdata['b_max']

    # xl, yl = box_s_eq * .5, box_s_eq * .5

    # # Limits of the cluster's region in Equatorial
    # xmin_cl, xmax_cl = c_lon - xl, c_lon + xl
    # ymin_cl, ymax_cl = c_lat - yl, c_lat + yl

    # These are the points that determine the range of *all* the frames
    ra_min, ra_max = fdata['ra_min'], fdata['ra_max']
    dec_min, dec_max = fdata['dec_min'], fdata['dec_max']

    # frame == 'galactic':
    box_s_eq = np.sqrt(2) * box_s_eq
    # Correct size in RA
    box_s_x = box_s_eq / np.cos(np.deg2rad(c_dec))

    xl, yl = box_s_x * .5, box_s_eq * .5

    # Limits of the cluster's region in Equatorial
    xmin_cl, xmax_cl = c_ra - xl, c_ra + xl
    ymin_cl, ymax_cl = c_dec - yl, c_dec + yl

    # Identify which frames contain the cluster region
    l2 = (xmin_cl, ymax_cl)  # Top left
    r2 = (xmax_cl, ymin_cl)  # Bottom right
    frame_intersec = []
    for i, xmin_fr_i in enumerate(ra_min):
        l1 = (xmin_fr_i, dec_max[i])  # Top left
        r1 = (ra_max[i], dec_min[i])  # Bottom right
        frame_intersec.append(doOverlap(l1, r1, l2, r2))
    frame_intersec = np.array(frame_intersec)

    data_in_files = list(fdata[frame_intersec]['filename'])
    # print(f"  Cluster is present in {len(data_in_files)} frames")

    return data_in_files, xmin_cl, xmax_cl, ymin_cl, ymax_cl


def doOverlap(l1, r1, l2, r2):
    """
    Given two rectangles, find if the given two rectangles overlap or not.
    Note that a rectangle can be represented by two coordinates, top left and
    bottom right.
    l1: Top Left coordinate of first rectangle.
    r1: Bottom Right coordinate of first rectangle.
    l2: Top Left coordinate of second rectangle.
    r2: Bottom Right coordinate of second rectangle.

    Source: https://www.geeksforgeeks.org/find-two-rectangles-overlap/
    """
    min_x1, max_y1 = l1
    max_x1, min_y1 = r1
    min_x2, max_y2 = l2
    max_x2, min_y2 = r2
    # If one rectangle is on left side of other
    if (min_x1 > max_x2 or min_x2 > max_x1):
        return False
    # If one rectangle is above other
    if(min_y1 > max_y2 or min_y2 > max_y1):
        return False
    return True


def query(
    c_ra, c_dec, box_s_eq, frames_path, max_mag, data_in_files, xmin_cl,
    xmax_cl, ymin_cl, ymax_cl, plx_min
):
    """
    """
    # Mag (flux) filter
    min_G_flux = 10**((max_mag - Zp_G) / (-2.5))

    all_frames = []
    for i, file in enumerate(data_in_files):
        with gzip.open(frames_path + file) as f:
            data = pd.read_csv(f, index_col=False)

            # lon, lat = radec2lonlat(data['ra'].values, data['dec'].values)
            # # print(lon, lat, data['l'], data['b'])
            # # lon, lat = data['l'].values, data['b'].values
            # mx = (lon >= xmin_cl) & (lon <= xmax_cl)
            # my = (lat >= ymin_cl) & (lat <= ymax_cl)

            mx = (data['ra'] >= xmin_cl) & (data['ra'] <= xmax_cl)
            my = (data['dec'] >= ymin_cl) & (data['dec'] <= ymax_cl)
            m_plx = data['parallax'] > plx_min
            m_gmag = data['phot_g_mean_flux'] > min_G_flux
            msk = (mx & my & m_plx & m_gmag)

            # print(f"{i+1}, {file} contains {msk.sum()} cluster stars")
            if msk.sum() == 0:
                continue

            all_frames.append(data[msk])
    all_frames = pd.concat(all_frames)

    # print(f"  {len(all_frames)} stars retrieved")

    c_ra, c_dec = c_ra, c_dec
    box_s_h = box_s_eq * .5
    gal_cent = radec2lonlat(c_ra, c_dec)

    if all_frames['l'].max() - all_frames['l'].min() > 180:
        # print("Frame wraps around 360 in longitude. Fixing..")

        lon = all_frames['l'].values
        if gal_cent[0] > 180:
            msk = lon < 180
            lon[msk] += 360
        else:
            msk = lon > 180
            lon[msk] -= 360
        all_frames['l'] = lon

    xmin_cl, xmax_cl = gal_cent[0] - box_s_h, gal_cent[0] + box_s_h
    ymin_cl, ymax_cl = gal_cent[1] - box_s_h, gal_cent[1] + box_s_h
    mx = (all_frames['l'] >= xmin_cl) & (all_frames['l'] <= xmax_cl)
    my = (all_frames['b'] >= ymin_cl) & (all_frames['b'] <= ymax_cl)
    msk = (mx & my)
    all_frames = all_frames[msk]

    all_frames = all_frames.rename(columns={
        'source_id': 'Source', 'ra': 'RA_ICRS', 'dec': 'DE_ICRS',
        'parallax': 'Plx', 'parallax_error': 'e_Plx',
        'pmra': 'pmRA', 'pmra_error': 'e_pmRA', 'b': 'GLAT',
        'pmdec': 'pmDE', 'pmdec_error': 'e_pmDE', 'l': 'GLON',
        'phot_g_mean_flux': 'FG', 'phot_g_mean_flux_error': 'e_FG',
        'phot_bp_mean_flux': 'FBP', 'phot_bp_mean_flux_error': 'e_FBP',
        'phot_rp_mean_flux': 'FRP', 'phot_rp_mean_flux_error': 'e_FRP',
        'radial_velocity': 'RV', 'radial_velocity_error': 'e_RV'}
    )

    # print(f"  {len(all_frames)} stars final")
    return all_frames


def uncertMags(data):
    """
    # Gaia DR3 zero points:
    https://www.cosmos.esa.int/web/gaia/dr3-passbands

    "The GBP (blue curve), G (green curve) and GRP (red curve) passbands are
    applicable to both Gaia Early Data Release 3 as well as to the full Gaia
    Data Release 3"
    """
    I_G, e_IG = data['FG'].values, data['e_FG'].values
    I_BP, e_IBP = data['FBP'].values, data['e_FBP'].values
    I_RP, e_IRP = data['FRP'].values, data['e_FRP'].values

    data['Gmag'] = Zp_G + -2.5 * np.log10(I_G)
    BPmag = Zp_BP + -2.5 * np.log10(I_BP)
    RPmag = Zp_RP + -2.5 * np.log10(I_RP)
    data['BP-RP'] = BPmag - RPmag

    e_G = np.sqrt(sigma_ZG_2 + 1.179 * (e_IG / I_G)**2)
    data['e_Gmag'] = e_G
    e_BP = np.sqrt(sigma_ZBP_2 + 1.179 * (e_IBP / I_BP)**2)
    e_RP = np.sqrt(sigma_ZRP_2 + 1.179 * (e_IRP / I_RP)**2)
    data['e_BP-RP'] = np.sqrt(e_BP**2 + e_RP**2)

    return data


def radec2lonlat(ra, dec):
    gc = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    lb = gc.transform_to('galactic')
    lon, lat = lb.l.value, lb.b.value
    return lon, lat

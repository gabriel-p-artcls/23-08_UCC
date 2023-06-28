
import sys
import os
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, GalacticLSR

UCC_path = "/home/gabriel/Github/UCC/add_New_DB/"
sys.path.insert(1, UCC_path + 'modules/')
import main_process_GDR3_query as G3Q


def main():
    """
    """
    GAIADR3_path = '/media/gabriel/backup/gabriel/GaiaDR3/'
    frames_path = GAIADR3_path + 'datafiles_parquet_G20/'
    BJ21 = pd.read_parquet("/home/gabriel/Descargas/BJ21.parquet")
    min_plx = 4

    df_all = pd.DataFrame()
    for i, file in enumerate(os.listdir(frames_path)):
        print(file)

        # Read Gaia file
        df = pd.read_parquet(frames_path + file)
        msk = df['parallax'].values > min_plx * 1.5
        df = df[msk]

        # Use BJ21 distances
        msk = np.isin(df['source_id'].values, BJ21['Source'].values)
        df = df[msk]
        idx = np.where(np.isin(BJ21['Source'].values, df['source_id'].values))[0]
        d_pc = BJ21['rpgeo'].values[idx]
        df['r_med_photogeo'] = d_pc

        x, y, z, v_x, v_y, v_z = H23_transf(df)
        df['X'], df['Y'], df['Z'] = x, y, z
        df['VX'], df['VY'], df['VZ'] = v_x, v_y, v_z

        # plx = 1000 / d_pc
        # ra, dec = df['ra'].values, df['dec'].values

        # # Velocities
        # pmra_t, pmdec_t = 4.74*df['pmra'].values/plx, 4.74*df['pmdec'].values/plx
        # df['pmRA_t'] = pmra_t
        # df['pmDE_t'] = pmdec_t
        # # Galactic coordinates.
        # eq = SkyCoord(ra=df['ra'].values*u.degree, dec=df['dec'].values*u.degree, frame='icrs')
        # lb = eq.transform_to('galactic')
        # lon = lb.l.wrap_at(180 * u.deg).radian * u.radian
        # lat = lb.b.radian * u.radian
        # coords = SkyCoord(l=lon, b=lat, distance=d_pc*u.pc, frame='galactic')
        # # Galactocentric coordinates.
        # c_glct = coords.transform_to(gc_frame)
        # x_pc, y_pc, z_pc = c_glct.x, c_glct.y, c_glct.z
        # # x_kpc, y_kpc, z_kpc = x_pc.value/1000, y_pc.value/1000, z_pc.value/1000
        # df['X'] = x_pc.value #_kpc
        # df['Y'] = y_pc.value #_kpc_kpc
        # df['Z'] = z_pc.value #_kpc_kpc

        df_all = pd.concat([df_all, df])
        if i in (500, 1000, 1500, 2000, 2500):
            df_all.to_parquet(GAIADR3_path + f"GaiaDR3_plx_{i}.parquet", index=False)
            df_all = pd.DataFrame()
    # Last file
    df_all.to_parquet(GAIADR3_path + f"GaiaDR3_plx_{i}.parquet", index=False)
    print("Finished")


def H23_transf(df):
    """
    """
    # Initialise a SkyCoord object for raw Gaia data
    coords = SkyCoord(ra=df['ra'].to_numpy() << u.deg, 
                      dec=df['dec'].to_numpy() << u.deg, 
                      pm_ra_cosdec=df['pmra'].to_numpy() << (u.mas / u.yr), 
                      pm_dec=df['pmdec'].to_numpy() << (u.mas / u.yr), 
                      distance=df['r_med_photogeo'].to_numpy() << u.pc,
                      radial_velocity=np.zeros(len(df)) << (u.m / u.s),
                      frame="icrs")

    # Transform to galactic standard of rest 
    # (may not be necessary to use the LSR, which simply adds a mean shift to
    # every velocity)
    coords_lsr = coords.transform_to(GalacticLSR())

    # Create a coordinate frame with these galactic coordinates
    # We assume that RVs are all zero, i.e. we get out tangential velocities
    # in 2D projected into 3 dimensions.
    coords_lsr_no_rv = SkyCoord(coords_lsr.l, 
                                coords_lsr.b, 
                                distance=coords_lsr.distance,
                                pm_l_cosb=coords_lsr.pm_l_cosb, 
                                pm_b=coords_lsr.pm_b, 
                                radial_velocity=np.zeros(len(df)) << (u.km / u.s),
                                frame=GalacticLSR())

    # Extract raw numpy arrays from these objects (without units)
    xyz = coords.cartesian
    x, y, z = xyz.x.value, xyz.y.value, xyz.z.value

    velocities = coords_lsr_no_rv.velocity
    v_x, v_y, v_z = velocities.d_x.value, velocities.d_y.value, velocities.d_z.value

    return x, y, z, v_x, v_y, v_z


def generate_clust():
    """
    """
    # Generate HSC_2971
    ra, dec, plx = 253.19315, -23.51814, 9.334
    d_pc = 1000 / plx

    cc = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, distance=d_pc*u.pc)
    x, y, z = cc.cartesian.x, cc.cartesian.y, cc.cartesian.z
    x, y, z = x.value, y.value, z.value
    print(x, y, z)

    xmin, xmax = x - 50, x + 50
    ymin, ymax = y - 50, y + 50
    zmin, zmax = z - 50, z + 50

    GAIADR3_path = '/media/gabriel/backup/gabriel/GaiaDR3/'
    frames_path = GAIADR3_path + 'datafiles_parquet_G20_plx_cut/'

    df_all = pd.DataFrame()
    for i in (500, 1000, 1500, 2000, 2500, 3385):
        print(i)
        df = pd.read_parquet(frames_path + f"GaiaDR3_plx_{i}.parquet")
        msk = (df['X'] > xmin) & (df['X'] < xmax) & (df['Y'] > ymin)\
            & (df['Y'] < ymax) & (df['Z'] < zmax) & (df['Z'] > zmin)
        df_all = pd.concat([df_all, df[msk]])

    df_all = df_all.rename(columns={
        'source_id': 'Source', 'ra': 'RA_ICRS', 'dec': 'DE_ICRS',
        'parallax': 'Plx', 'parallax_error': 'e_Plx',
        'pmra': 'pmRA', 'pmra_error': 'e_pmRA', 'b': 'GLAT',
        'pmdec': 'pmDE', 'pmdec_error': 'e_pmDE', 'l': 'GLON',
        'phot_g_mean_flux': 'FG', 'phot_g_mean_flux_error': 'e_FG',
        'phot_bp_mean_flux': 'FBP', 'phot_bp_mean_flux_error': 'e_FBP',
        'phot_rp_mean_flux': 'FRP', 'phot_rp_mean_flux_error': 'e_FRP',
        'radial_velocity': 'RV', 'radial_velocity_error': 'e_RV'}
    )
    df_all = G3Q.uncertMags(df_all)
    df_all.to_parquet("/home/gabriel/Descargas/hsc2971.parquet", index=False)


if __name__ == '__main__':
    generate_clust()


import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd

"""
Combine all the new databases with the old databases (already merged)
into a single catalog

Keys:


-QIN23: Name,OC-Flag,GLON,GLAT,RAdeg,DEdeg,pmRA,e_pmRA,pmDE,e_pmDE,plx,e_plx,
RV,e_RV,f_RV,o_RV,N,m-M,logAge,E(B-V),rh,rc,e_rc,rt,e_rt,ro,e_ro,re,e_re,
Ref-Name,Ref

-LI23: id,ra,sig_ra,dec,sig_dec,l,b,plx,pmra,pmdec,rsc,nstar

-LI22: id,ra,dec,plx,sig_plx,pmRA,sig_pmRA,pmDE,sig_pmDE,rsc,m_M,E(V-I),Z,t,
f_bin,f_rot

-HE22_2: 'CWNU_id', 'l', 'b', 'sig_l', 'sig_b', 'n', 'plx', 'sig_plx', 'pmx',
'sig_pmx', 'pmy', 'sig_pmy', 'rv', 'sig_rv', 'nrv', 'logage', 'A0',
'm-M', 'class'

-HE22_1: 'Cluster', 'GLON', 'e_LON', 'GLAT', 'e_LAT', 'Num', 'plx', 'e_plx',
'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'logA', 'A0', '_RA.icrs', '_DE.icrs'

-HE22: 'Cluster', 'GLON', 's_GLON', 'GLAT', 's_GLAT', 'Plx', 's_Plx', 'pmRA',
's_pmRA', 'pmDE', 's_pmDE', 'RV', 's_RV', 'o_RV', 'm-M', 'AG', 'logAge',
'Z', 'N70', 'minProb', 'NminProb', 'Class', '_RA.icrs', '_DE.icrs'

-HAO22: 'Cluster', 'RA_ICRS', 'e_RA_ICRS', 'DE_ICRS', 'e_DE_ICRS', 'GLON',
'e_GLON', 'GLAT', 'e_GLAT', 'adeg', 'plx', 'e_plx', 'pmRA', 'e_pmRA',
'pmDE', 'e_pmDE', 'age', 'e_age', 'AG', 'Z', 'RV', 's_RV', 'N', 'o_RV',
'_RA.icrs', '_DE.icrs'

-CASTRO22: 'Seq', 'Cluster', 'RA_ICRS', 's_RA_ICRS', 'DE_ICRS', 's_DE_ICRS',
'GLON', 's_GLON', 'GLAT', 's_GLAT', 'plx', 's_plx', 'pmRA', 's_pmRA',
'pmDE', 's_pmDE', 'RV', 's_RV', 'Nmemb', 'NmembRV', 'Flag', 'logAge',
'Dist', 'AV', '_RA.icrs', '_DE.icrs'

-HE21: "OC","GLON","GLAT","r","Nstar","Ncore","plx","s_plx","pmRA","s_pmRA",
"pmDE","s_pmDE","m-M","logAge","AG","Z"

-CASTRO20: 'Cluster', 'RA_ICRS', 'e_RA_ICRS', 'DE_ICRS', 'e_DE_ICRS', 'GLON',
'e_GLON', 'GLAT', 'e_GLAT', 'plx', 'e_plx', 'pmRA', 'e_pmRA', 'pmDE',
'e_pmDE', 'Note', '_RA.icrs', '_DE.icrs'

-FERREIRA20: "Name","GLON","GLAT","Dist","e_Dist","m-M","e_m-M","logt",
"e_logt","E(B-V)","e_E(B-V)","rlim","rc","rt","pmRAcl","pmDEcl","plxcl","Nmemb"

-CASTRO19: "Cluster","RA_ICRS","e_RA_ICRS","DE_ICRS","e_DE_ICRS","GLON",
"e_GLON","GLAT","e_GLAT","Plx","e_Plx","pmRA*","e_pmRA*","pmDE","e_pmDE",
"Simbad"

-LIUPANG19: 'ID', 'GLON', 'e_GLON', 'GLAT', 'e_GLAT', 'r', 'plx', 'e_plx',
'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'nTot', 'Age', 'e_Age', 'Z', '_RA.icrs',
'_DE.icrs'

-SIM19: "ID","lon","lat","pmRA*","e_pmRA*","pmDE","e_pmDE","dist_pc","e_dist",
"rGMM","rc","e_rc","N","e_N","log(age)"

-CASTRO18: Cluster,GLON,e_GLON,GLAT,e_GLAT,Plx,e_Plx,pmRA*,e_pmRA*,pmDE,
e_pmDE,RV,e_RV,N,NRV,Simbad,_RA.icrs,_DE.icrs
"""

# Data stored per cluster (per database):
# ID, RA, DEC, GLON, GLAT, plx, pmRA, pmDE,
# Dist, AV, Age, [Fe/H],
# radec_flag, dmod_flag, age_flag, feh_flag
#
# radec_flag: True (RA, DEC) / False (lon, lat)
# dmod_flag: mag / kpc / pc
# age_flag: log / Gyr
# feh_flag: True (FeH) / False (Z)


dbs_folder = "../0_data/databases/new/"

dbs_names = {
    'QIN23': ['Name', 'RAdeg', 'DEdeg', 'plx', 'pmRA', 'pmDE', 'm-M', 'E(B-V)',
              'logAge', None, [True, 'mag', 'log', True]],
    # 'LI23': ['id', 'ra', 'dec', 'plx', 'pmra', 'pmdec', None, None,
    #          None, None, [True, 'mag', 'Gyr', False]],
    'LI22': ['id', 'ra', 'dec', 'plx', 'pmRA', 'pmDE', 'm_M', 'E(V-I)',
             't', 'Z', [True, 'mag', 'Gyr', False]],
    'HE22_2': ['CWNU_id', 'l', 'b', 'plx', 'pmx', 'pmy', 'm-M', 'A0', 'logage',
               None, [False, 'mag', 'log', True]],
    'HE22_1': ['Cluster', 'GLON', 'GLAT', 'plx', 'pmRA', 'pmDE', None, 'A0',
               'logA', None, [False, 'mag', 'log', True]],
    'HE22': ['Cluster', 'GLON', 'GLAT', 'Plx', 'pmRA', 'pmDE', 'm-M', 'AG',
             'logAge', None, [False, 'mag', 'log', True]],
    'HAO22': ['Cluster', 'RA_ICRS', 'DE_ICRS', 'plx', 'pmRA', 'pmDE', None,
              'AG', 'age', None, [True, 'mag', 'log', True]],
    'CASTRO22': ['Cluster', 'RA_ICRS', 'DE_ICRS', 'plx', 'pmRA', 'pmDE',
                 'Dist', 'AV', 'logAge', None, [True, 'pc', 'log', True]],
    'HE21': ['OC', 'GLON', 'GLAT', 'plx', 'pmRA', 'pmDE', "m-M", "logAge",
             'AG', "Z", [False, 'mag', 'log', False]],
    'CASTRO20': ['Cluster', 'RA_ICRS', 'DE_ICRS', 'plx', 'pmRA', 'pmDE', None,
                 None, None, None, [True, 'mag', 'log', True]],
    'FERREIRA20': ['Name', 'GLON', 'GLAT', 'plxcl', 'pmRAcl', 'pmDEcl', "m-M",
                   "logt", "E(B-V)", None, [False, 'mag', 'log', True]],
    'CASTRO19': ['Cluster', 'RA_ICRS', 'DE_ICRS', 'Plx', 'pmRA*', 'pmDE', None,
                 None, None, None, [True, 'mag', 'log', True]],
    'LIUPANG19': ['ID', 'GLON', 'GLAT', 'plx', 'pmRA', 'pmDE', None, None,
                  'Age', 'Z', [False, 'mag', 'Gyr', True]],
    'SIM19': ['ID', 'lon', 'lat', None, 'pmRA*', 'pmDE', 'dist_pc', None,
              'log(age)', None, [False, 'pg', 'log', True]],
    'CASTRO18': ['Cluster', 'GLON', 'GLAT', 'Plx', 'pmRA*', 'pmDE', None,
                 None, None, None, [False, 'mag', 'log', True]],
}

DB_ids = {
    'CASTRO18': 7, 'CASTRO19': 8, 'SIM19': 9, 'LIUPANG19': 8, 'CASTRO20': 'a',
    'FERREIRA20': 'b', 'HE21': 'c', 'CASTRO22': 'd', 'HAO22': 'e', 'HE22': 'f',
    'HE22_1': 'g', 'HE22_2': 'g', 'LI22': 'i',  'QIN23': 'j',
}


def main():
    """
    """
    final_DB = {
        'DB': [], 'ID': [], 'RA_ICRS': [], 'DE_ICRS': [], 'GLON': [],
        'GLAT': [], 'plx': [], 'pmRA': [], 'pmDE': [], 'distmod': [],
        'extin': [], 'logage': [], 'FeH': []
    }

    for DB, cols in dbs_names.items():
        data = pd.read_csv(dbs_folder + DB + ".csv", index_col=False)
        print(DB, len(data))

        radec_flag, dmod_flag, age_flag, feh_flag = cols[-1]

        final_DB['DB'] += [DB_ids[DB] for _ in range(len(data))]

        if DB == "HE22_2":
            names = ["CWNU" + str(_) for _ in data[cols[0]]]
        elif DB in ('LIUPANG19', 'HE21', 'HE22'):
            names = [DB[:4] + str(_) for _ in data[cols[0]]]
        else:
            names = list(data[cols[0]])
        final_DB['ID'] += names

        # Transfor to (ra, dec) if required
        if radec_flag is False:
            lon, lat = data[cols[1]], data[cols[2]]
            ra, dec = lonlat2radec(lon, lat)
        else:
            ra, dec = data[cols[1]], data[cols[2]]
            lon, lat = radec2lonlat(ra, dec)

        final_DB['RA_ICRS'] += list(ra)
        final_DB['DE_ICRS'] += list(dec)
        final_DB['GLON'] += list(lon)
        final_DB['GLAT'] += list(lat)

        if cols[3] is None:
            final_DB['plx'] += [np.nan for _ in range(len(data))]
        else:
            final_DB['plx'] += list(data[cols[3]])
        final_DB['pmRA'] += list(data[cols[4]])
        final_DB['pmDE'] += list(data[cols[5]])

        if cols[6] is None:
            final_DB['distmod'] += [np.nan for _ in range(len(data))]
        else:
            if dmod_flag == 'mag':
                final_DB['distmod'] += list(data[cols[6]])
            else:
                distmod = dist2dm(data[cols[6]], dmod_flag)
                final_DB['distmod'] += list(distmod)

        if cols[7] is None:
            final_DB['extin'] += [np.nan for _ in range(len(data))]
        else:
            final_DB['extin'] += list(data[cols[7]])

        if cols[8] is None:
            final_DB['logage'] += [np.nan for _ in range(len(data))]
        else:
            if age_flag == 'log':
                final_DB['logage'] += list(data[cols[8]])
            else:
                logage = round(np.log10(data[cols[8]] * 10**9), 2)
                final_DB['logage'] += list(logage)

        if cols[9] is None:
            final_DB['FeH'] += [np.nan for _ in range(len(data))]
        else:
            if feh_flag is True:
                final_DB['FeH'] += list(data[cols[9]])
            else:
                FeH = round(np.log10(data[cols[9]] / 0.0152), 2)
                final_DB['FeH'] += list(FeH)

    # Save to file
    pd.DataFrame(final_DB).to_csv('NEW_DBs.csv', na_rep='nan', index=False)
    breakpoint()


def lonlat2radec(lon, lat):
    gc = SkyCoord(l=lon * u.degree, b=lat * u.degree, frame='galactic')
    rd = gc.transform_to('fk5')
    ra, dec = rd.ra.value, rd.dec.value
    return np.round(ra, 5), np.round(dec, 5)


def radec2lonlat(ra, dec):
    gc = SkyCoord(ra=ra * u.degree, dec=dec * u.degree)
    lb = gc.transform_to('galactic')
    lon, lat = lb.l.value, lb.b.value
    return np.round(lon, 5), np.round(lat, 5)


def dist2dm(dist, unit='kpc'):
    """
    Transform distances to distance modulus
    """
    if unit == 'kpc':
        dist *= 1000
    distmod = np.round(-5 + 5 * np.log10(dist), 2)
    return distmod


if __name__ == '__main__':
    # plt.style.use('science')
    main()

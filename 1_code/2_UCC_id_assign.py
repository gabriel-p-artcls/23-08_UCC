
import numpy as np
from string import ascii_lowercase
import pandas as pd
import csv


path_io = "../2_pipeline/"
db_in = '1_combined_DBs.csv'
db_out = '2_standard_names_DB.csv'


def main():
    """
    Standardize all names + add UCC ids
    """
    dbs = pd.read_csv(path_io + db_in)

    no_dup_names = rm_dup_names(dbs['ID'])
    dbs['ID'] = no_dup_names

    fnames = assign_fname(dbs['ID'])
    dbs['fnames'] = fnames
    dup_check(fnames, dbs)

    ucc_ids = assign_UCC_ids(dbs['GLON'], dbs['GLAT'])
    dbs['UCC_ID'] = ucc_ids
    dup_check(ucc_ids, dbs)

    dbs.to_csv(path_io + db_out, na_rep='nan', index=False,
               quoting=csv.QUOTE_NONNUMERIC)
    print("Final database written to file")


def rm_dup_names(all_names):
    """
    This removes duplicates such as "NGC_2516" and "NGC 2516" 
    """
    no_dup_names = []
    for names in all_names:
        names = names.split(',')
        names_temp = []
        for name in names:
            name = name.strip()
            name = name.replace('_', ' ')
            names_temp.append(name)
        # Equivalent to set() but maintains order
        no_dup_names.append(','.join(list(dict.fromkeys(names_temp))))

    return no_dup_names


def assign_fname(all_names):
    """
    Assign names used for files and urls
    """
    fnames = []
    for names in all_names:
        names = names.split(',')
        names_temp = []
        for name in names:
            name = name.strip()
            # We replace '+' with 'p' to avoid duplicating names for clusters
            # like 'Juchert J0644.8-0925' and 'Juchert_J0644.8+0925'
            name = name.lower().replace('_', '').replace(' ', '').replace(
                '-', '').replace('.', '').replace('+', 'p')
            names_temp.append(name)

        fname = preferred_names(names_temp)

        fnames.append(fname)

    return fnames


def preferred_names(names_temp):
    """
    Use naming conventions according to this list of preferred names
    """
    names_lst = (
        'blanco', 'ngc', 'melotte', 'trumpler', 'ruprecht', 'berkeley',
        'pismis', 'vdbh', 'loden', 'kronberger', 'collinder', 'harvard',
        'eso', 'ascc')

    fname = names_temp[0]
    if len(names_temp) > 1:
        for id_prefer in names_lst:
            for name in names_temp:
                if id_prefer in name:
                    fname = name
                    return fname
    return fname


def assign_UCC_ids(glon, glat):
    """
    Format: UCC GXXX.X+YY.Y
    """
    lonlat = np.array([glon, glat]).T
    lonlat = trunc(lonlat)
    
    ucc_ids = []
    for idx, ll in enumerate(lonlat):
        lon, lat = str(ll[0]), str(ll[1])

        if ll[0] < 10:
            lon = '00' + lon
        elif ll[0] < 100:
            lon = '0' + lon

        if ll[1] >= 10:
            lat = '+' + lat
        elif ll[1] < 10 and ll[1] > 0:
            lat = '+0' + lat
        elif ll[1] == 0:
            lat = '+0' + lat.replace('-', '')
        elif ll[1] < 0 and ll[1] >= -10:
            lat = '-0' + lat[1:]
        elif ll[1] < -10:
            pass

        ucc_id = 'UCC G' + lon + lat

        i = 0
        while True:
            if i > 25:
                ucc_id += "ERROR"
                print("ERROR NAMING")
                break
            if ucc_id in ucc_ids:
                if i == 0:
                    # Add a letter to the end
                    ucc_id += ascii_lowercase[i]
                else:
                    # Replace last letter
                    ucc_id = ucc_id[:-1] + ascii_lowercase[i]
                i += 1
            else:
                break
        # if i > 0:
        #     print(ucc_id, dbs['ID'][idx], dbs['GLON'][idx], dbs['GLAT'][idx])

        ucc_ids.append(ucc_id)

    return ucc_ids


def trunc(values, decs=1):
    return np.trunc(values*10**decs)/(10**decs)


def dup_check(names, db):
    """
    Check for duplicates in 'names' list
    """
    for i, cl0 in enumerate(names):
        for j, cl1 in enumerate(names[i + 1:]):
            if cl0 == cl1:
                print(i, i + 1 + j, cl0, dbs['ID'][i], dbs['ID'][i + 1 + j])
                break
    print("End duplicates check")


if __name__ == '__main__':
    main()

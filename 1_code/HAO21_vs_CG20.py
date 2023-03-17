
"""
Compare the clusters in the HAO21 database with those matched in the CG20
database.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii


input_folder = '/home/gabriel/Github/sc_full_cat/input/'


def main():
    """
    """
    hao21 = ascii.read(input_folder + 'HAO21.dat', delimiter=',')
    cg20 = ascii.read(input_folder + 'CG20.dat', delimiter=',')

    names, cl_match_hao, cl_match_cg = [], [], []
    for i, cl_a in enumerate(hao21):
        print(cl_a['Cluster'])
        for j, cl_b in enumerate(cg20):
            cl_a_n, cl_b_n = cl_a['Cluster'], cl_b['Name']
            r = similar(cl_a_n, cl_b_n)
            if r is True:
                names.append([cl_a['Cluster'], cl_b['Name']])
                cl_match_hao.append(list(
                    cl_a['RA_ICRS', 'DE_ICRS', 'plx', 'pmRA*', 'pmDE']))
                cl_match_cg.append(list(
                    cl_b['RA', 'DEC', 'plx', 'pmRA', 'pmDE']))
    names = np.array(names)
    cl_match_hao = np.array(cl_match_hao).T
    cl_match_cg = np.array(cl_match_cg).T

    print((abs(cl_match_hao[2] - cl_match_cg[2]) / cl_match_hao[2] > .5).sum())
    breakpoint()


def similar(a, b):
    a1 = a.lower().replace('_', '').replace(' ', '')
    b1 = b.lower().replace('_', '').replace(' ', '')
    if a1 == b1:
        return True

    a1 = a.replace('vdBergh-Hagen_', 'BH').replace('_', '').replace(' ', '')
    b1 = b.replace('vdBergh-Hagen_', 'BH').replace('_', '').replace(' ', '')
    if a1 == b1:
        return True

    return False


if __name__ == '__main__':
    # plt.style.use('science')
    main()

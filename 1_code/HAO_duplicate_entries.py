
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii


def main():
    """
    """
    db = ascii.read(
        '/home/gabriel/Github/sc_full_cat/input/HAO21.dat', delimiter=',')

    cl_idx = 1
    for i, cl_a in enumerate(db):
        for j, cl_b in enumerate(db[cl_idx:]):
            cl_a_n, cl_b_n = cl_a['Cluster'], cl_b['Cluster']
            r = similar(cl_a_n, cl_b_n)
            if r is True:
                print()
                print(i, list(cl_a['Cluster','RA_ICRS','DE_ICRS','GLON','GLAT','plx','pmRA*','pmDE','N','logt','SimbadName']))
                print(cl_idx+j+1, list(cl_b['Cluster','RA_ICRS','DE_ICRS','GLON','GLAT','plx','pmRA*','pmDE','N','logt','SimbadName']))
        cl_idx += 1


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

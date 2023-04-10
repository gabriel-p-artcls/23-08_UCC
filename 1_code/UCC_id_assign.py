
import numpy as np
from string import ascii_lowercase


def main():
    """
    """


def assign_UCC_ids(final_DB):
    """
    UCC GXXX.X+YY.Y
    """
    lonlat = np.array([final_DB['GLON'], final_DB['GLAT']]).T
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
        if i > 0:
            print(ucc_id, final_DB['ID'][idx], final_DB['GLON'][idx], final_DB['GLAT'][idx])

        ucc_ids.append(ucc_id)
    breakpoint()

    return final_DB


def trunc(values, decs=1):
    return np.trunc(values*10**decs)/(10**decs)


if __name__ == '__main__':
    main()

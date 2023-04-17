
import pandas as pd


def main():
    """
    """
    df = pd.read_csv("../0_data/databases/BICA19.csv")

    all_names, names_orig = [], []
    for i, cl in enumerate(df['Name']):
        names = cl.split(',')
        for name in names:

            name = cluster_rename(name)

            all_names.append(name.lower().replace('_', '').replace(
                               ' ', '').replace('-', ''))
            names = cl.split(',')
            names = ",".join([_.strip() for _ in names])
            names_orig.append(names)

    for i, name in enumerate(all_names):
        for j, name_2 in enumerate(all_names[i + 1:]):
            if name == name_2:
                print(i, i + 1 + j, name, names_orig[i], "|", names_orig[i + 1 + j])
    print("Finished")


def cluster_rename(name):
    """
    Standardize the naming of these clusters watching for 

    FSR XXX w leading zeros
    FSR XXX w/o leading zeros
    FSR_XXX w leading zeros
    FSR_XXX w/o leading zeros

    ESO XXX-YY w leading zeros
    ESO XXX-YY w/o leading zeros
    ESO_XXX_YY w leading zeros
    ESO_XXX_YY w/o leading zeros
    ESO_XXX-YY w leading zeros
    ESO_XXX-YY w/o leading zeros
    """
    name = name.strip()

    if name.startswith("FSR"):
        if '_' in name:
            n2 = name.split('_')[1]
        else:
            n2 = name.split(' ')[1]
        n2 = int(n2)
        if n2 < 10:
            n2 = '000' + str(n2)
        elif n2 < 100:
            n2 = '00' + str(n2)
        elif n2 < 1000:
            n2 = '0' + str(n2)
        else:
            n2 = str(n2)
        name = "FSR_" + n2

    if name.startswith("ESO"):
        if ' ' in name[4:]:
            n1, n2 = name[4:].split(' ')
        elif '_' in name[4:]:
            n1, n2 = name[4:].split('_')
        elif '' in name[4:]:
            n1, n2 = name[4:].split('-')

        n1 = int(n1)
        if n1 < 10:
            n1 = '00' + str(n1)
        elif n1 < 100:
            n1 = '0' + str(n1)
        else:
            n1 = str(n1)
        n2 = int(n2)
        if n2 < 10:
            n2 = '0' + str(n2)
        else:
            n2 = str(n2)
        name = "ESO_" + n1 + '_' + n2

    return name


if __name__ == '__main__':
    # plt.style.use('science')
    main()

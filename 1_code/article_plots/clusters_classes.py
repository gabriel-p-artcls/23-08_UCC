
import matplotlib.pyplot as plt
import pandas as pd
import scienceplots
plt.style.use('science')

clpath = "/home/gabriel/Github/UCC/"
final_dbs_path = clpath + "add_New_DB/UCC_cat_20230702.csv"

print("Reading fastMP output file...")
fastMP_db = pd.read_csv(final_dbs_path)
fnames = [_.split(';')[0] for _ in fastMP_db['fnames']]


def main():
    """
    """
    fig = plt.figure(figsize=(12, 10))

    for q, cl_name in enumerate(('NGC 2516', 'Haffner 9', 'Czernik 10', 'LISC 3630')):

        fname = cl_name.lower().replace('_', '').replace(
            ' ', '').replace('-', '').replace('.', '').replace('+', 'p')

        # Read fastMP output catalogue
        i = fnames.index(fname)
        row = fastMP_db.iloc[i]
        # Search for the cluster file
        file_name = clpath + row['quad'] + '/datafiles/' + fname + '.parquet'
        fastMP_cl = pd.read_parquet(file_name)

        plot(q*4 + 1, fastMP_cl, cl_name)

    # plt.show()
    # plt.subplots_adjust(hspace=.025)
    fig.tight_layout()
    plt.savefig("clusters_classif.png", dpi=300)


def plot(q, clust, lbl):
    """
    """
    ax = plt.subplot(4, 4, q)
    plt.hist(clust['Plx'], alpha=.5, label=lbl)
    if q >= 13:
        plt.xlabel("Plx [mas]")
    plt.ylabel("N")
    plt.text(0.65, 0.82, lbl, transform=ax.transAxes)
    if q == 1:
        plt.text(0.71, 0.75, r"\textbf{(AA)}", transform=ax.transAxes)
    if q == 5:
        plt.text(0.71, 0.75, r"\textbf{(BB)}", transform=ax.transAxes)
    if q == 9:
        plt.text(0.71, 0.75, r"\textbf{(CC)}", transform=ax.transAxes)
    if q == 13:
        plt.text(0.71, 0.75, r"\textbf{(DD)}", transform=ax.transAxes)

    ax = plt.subplot(4, 4, q + 1)
    plt.scatter(clust['GLON'], clust['GLAT'], ec='w', lw=.2, alpha=.7, c=clust['probs'])
    if q >= 13:
        plt.xlabel("lon [deg]")
    plt.ylabel("lat [deg]")

    ax = plt.subplot(4, 4, q + 2)
    plt.scatter(clust['pmRA'], clust['pmDE'], ec='w', lw=.2, alpha=.7, c=clust['probs'])
    if q >= 13:
        plt.xlabel("pmRA [mas/yr]")
    plt.ylabel("pmDE [mas/yr]")

    ax = plt.subplot(4, 4, q + 3)
    plt.scatter(clust['BP-RP'], clust['Gmag'], ec='w', lw=.2, alpha=.7, c=clust['probs'])
    if q >= 13:
        plt.xlabel("BP-RP")
    plt.ylabel("G")
    plt.gca().invert_yaxis()
    plt.colorbar()


if __name__ == '__main__':
    main()

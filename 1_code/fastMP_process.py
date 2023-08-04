
import os
import sys

cluster_run = True
try:
    import call_fastMP
except ModuleNotFoundError:
    cluster_run = False

if cluster_run is False:
    # LOCAL RUN
    UCC_path = "/home/gabriel/Github/UCC/add_New_DB/"
    # insert at 1, 0 is the script path (or '' in REPL)
    # Path to modules
    sys.path.insert(1, UCC_path + "modules/")
    import call_fastMP
    import main_process_GDR3_query as G3Q
    # Path to fastMP
    sys.path.insert(1, '/home/gabriel/Github/fastmp/')
    from fastmp import fastMP
    GAIADR3_path = '/media/gabriel/backup/gabriel/GaiaDR3/'
    frames_path = GAIADR3_path + 'datafiles_parquet_G20/'
    frames_ranges = GAIADR3_path + 'files_G20/frame_ranges.txt'
    files = os.listdir(UCC_path)
    for i, _ in enumerate(files):
        if 'UCC_cat' in _ and _.endswith(".csv"):
            idx = i
    UCC_fname = files[idx]
    UCC_cat = UCC_path + UCC_fname
    GCs_cat = UCC_path + "databases/globulars.csv"
else:
    # CLUSTER RUN
    import call_fastMP
    import main_process_GDR3_query as G3Q
    from fastmp import fastMP
    # Path to the database with Gaia DR3 data
    GAIADR3_path = '/home/gperren/GaiaDR3/'
    frames_path = GAIADR3_path + 'datafiles_parquet_G20/'
    # File that contains the regions delimited by each frame (in this folder)
    frames_ranges = GAIADR3_path + 'files_G20/frame_ranges.txt'
    # Full database of clusters (in this folder)
    files = os.listdir()
    for i, _ in enumerate(files):
        if _.endswith("_in.csv"):
            idx = i
    UCC_cat = files[idx]
    GCs_cat = "globulars.csv"

out_path = "out/"

# Run ID
# ij = int(os.path.basename(__file__).split('.')[0].split('_')[-1]) - 1
try:
    ij = int(sys.argv[1]) - 1
except IndexError:
    ij = -1


def run(Nj=10):
    """
    """
    # Read data
    frames_data, df_UCC, df_gcs = call_fastMP.read_input(
        frames_ranges, UCC_cat, GCs_cat)

    # Split into N_j jobs
    N_r = int(len(df_UCC) / Nj + 2)
    if ij >= 0:
        clusters_list = df_UCC[ij*N_r:(ij + 1)*N_r]
    else:
        clusters_list = df_UCC

    call_fastMP.run(
        fastMP, G3Q, frames_path, frames_data, df_UCC, df_gcs, UCC_cat,
        out_path, clusters_list)


if __name__ == '__main__':
    run()


import pandas as pd


dbs_folder = '/media/gabriel/rest/Dropbox_nosync/Papers/2023/GDR3_members/0_data/databases/old/'


data = pd.read_csv(
    dbs_folder + "HAO21" + ".csv", sep=',', comment='#', index_col=False)

# Replace FSR clusters
names = []
for cl in data['Cluster']:
    if cl.startswith("FSR"):
        try:
            if cl[4] != "0":
                names.append(cl)
                continue
        except:
            print("here")
            breakpoint()
        cl_id = int(cl[4:])
        names.append("FSR_" + str(cl_id))
    else:
        names.append(cl)
data['Cluster'] = names
data.to_csv(dbs_folder + 'HAO21_2.csv', index=False)

# # Replace ESO clusters
# names = []
# for cl in data['name']:
#     if cl.startswith("ESO"):
#         cl_id = cl[4:].split('-')
#         cl_id = str(int(cl_id[0])) + '_' + str(int(cl_id[1]))
#         names.append("ESO_" + str(cl_id))
#         print(cl, "ESO_" + str(cl_id))
#     else:
#         names.append(cl)
# data['name'] = names
# data.to_csv(dbs_folder + 'MWSC_2.csv', index=False)


# data = pd.read_csv(
#     dbs_folder + "MWSC" + ".dat", sep=',', comment='#', index_col=False)

# # Replace FSR clusters
# names = []
# for cl in data['name']:
#     if cl.startswith("FSR"):
#         if cl[4] != "0":
#             names.append(cl)
#             continue
#         cl_id = int(cl[4:])
#         names.append("FSR_" + str(cl_id))
#     else:
#         names.append(cl)
# data['name'] = names
# data.to_csv(dbs_folder + 'MWSC_2.csv', index=False)

# # Replace ESO clusters
# names = []
# for cl in data['name']:
#     if cl.startswith("ESO"):
#         cl_id = cl[4:].split('-')
#         cl_id = str(int(cl_id[0])) + '_' + str(int(cl_id[1]))
#         names.append("ESO_" + str(cl_id))
#         print(cl, "ESO_" + str(cl_id))
#     else:
#         names.append(cl)
# data['name'] = names
# data.to_csv(dbs_folder + 'MWSC_2.csv', index=False)



# data = pd.read_csv(
#     dbs_folder + "CG20" + ".dat", sep=',', comment='#', index_col=False)

# # Replace FSR clusters
# names = []
# for cl in data['Name']:
#     if cl.startswith("FSR"):
#         if cl[4] != "0":
#             names.append(cl)
#             continue
#         cl_id = int(cl[4:])
#         names.append("FSR_" + str(cl_id))
#     else:
#         names.append(cl)
# data['Name'] = names
# data.to_csv(dbs_folder + 'CG20_2.csv', index=False)

# # Replace ESO clusters
# names = []
# for cl in data['Name']:
#     if cl.startswith("ESO"):
#         cl_id = cl[4:].split('_')
#         cl_id = str(int(cl_id[0])) + '_' + str(int(cl_id[1]))
#         names.append("ESO_" + str(cl_id))
#         print(cl, "ESO_" + str(cl_id))
#     else:
#         names.append(cl)
# breakpoint()
# data['Name'] = names
# data.to_csv(dbs_folder + 'CG20_2.csv', index=False)



# data = pd.read_csv(
#     dbs_folder + "DIAS21" + ".dat", sep=',', comment='#', index_col=False)

# # Replace FSR clusters
# names = []
# for cl in data['Cluster']:
#     if cl.startswith("FSR"):
#         if cl[4] != "0":
#             names.append(cl)
#             continue
#         cl_id = int(cl[4:])
#         names.append("FSR_" + str(cl_id))
#     else:
#         names.append(cl)
# data['Cluster'] = names
# data.to_csv(dbs_folder + 'DIAS21_2.csv', index=False)

# # Replace ESO clusters
# names = []
# for cl in data['Cluster']:
#     if cl.startswith("ESO"):
#         cl_id = cl[4:].split('_')
#         cl_id = str(int(cl_id[0])) + '_' + str(int(cl_id[1]))
#         names.append("ESO_" + str(cl_id))
#         print(cl, "ESO_" + str(cl_id))
#     else:
#         names.append(cl)
# data['Cluster'] = names
# data.to_csv(dbs_folder + 'DIAS21_2.csv', index=False)
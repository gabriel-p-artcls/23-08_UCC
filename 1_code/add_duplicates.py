
import csv
import pandas as pd
import sys
sys.path.insert(1, '/home/gabriel/Github/UCC/add_New_DB/modules/')
import duplicates_id

"""
Identify possible duplicates (and assign a probability) using
the positions estimated with the most likely members.
"""

# This is the catalogue *after* running the 'fastMP_process' script (that
# calls the 'call_fastMP' module), because it needs to have the positions
# estimated with the most likely members
UCC_cat = "/home/gabriel/Github/UCC/add_New_DB/UCC_cat_20230517.csv"
df_UCC = pd.read_csv(UCC_cat)

print("Finding final duplicates and their probabilities...")
dups_fnames, dups_probs = duplicates_id.run(df_UCC)
df_UCC['dups_fnames'] = dups_fnames  # This column is rewritten here
df_UCC['dups_probs'] = dups_probs
df_UCC.to_csv(
    UCC_cat, na_rep='nan', index=False, quoting=csv.QUOTE_NONNUMERIC)
print(f"File {UCC_cat} updated")

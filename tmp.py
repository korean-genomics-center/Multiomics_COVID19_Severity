# %%
import glob
import os

import pandas as pd

# %%
path_metatbl = "/BiO/Access/kyungwhan1998/Infectomics/Resources/Data/Methylation/MetaTable_Incl_EHR_LOO/metatable_*.cpg_table_file_path.Severe.tsv"
list_metatbl = glob.glob(path_metatbl, recursive=True)
for metatbl in list_metatbl: 
    df_metatbl = pd.read_csv(metatbl, sep="\t")
    num_samples = len(df_metatbl)
    print(os.path.basename(metatbl), num_samples)
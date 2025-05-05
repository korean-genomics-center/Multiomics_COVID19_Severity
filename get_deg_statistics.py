# %%
import os 
from pathlib import Path
import pandas as pd

# %%
WORKDIR = str(Path(os.path.abspath(__file__)).parents[3])
path_exp = f"{WORKDIR}/Results/5_deg/Visit1_Severe__Visit1_Mild_20240327.tsv.normcount"
path_list_gene = f"{WORKDIR}/Resources/Data/Methylation/DEG_LOO/list_genes_over_loo_overlap_threshold_7.txt"
file_deg = "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/Visit1_Severe__Visit1_Mild_20240327.tsv"

# %%
with open(path_list_gene, mode="r") as fr:
    list_deg_loo = list(map(lambda x: x.rstrip("\n"), fr.readlines()))

df_deg = pd.read_csv(file_deg, sep="\t")
df_deg_set_index = df_deg.set_index("ID")
df_deg_loo = df_deg_set_index.loc[list_deg_loo]
df_deg_loo

# %%

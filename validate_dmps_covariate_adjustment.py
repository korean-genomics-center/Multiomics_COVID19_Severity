# %%
import glob
import os
from pathlib import Path

import pandas as pd

# %%
WORKDIR = str(Path(os.path.abspath(__file__)).parents[3])
path_methyl_markers = f"{WORKDIR}/Resources/Scripts/Final/list_methyl_markers.txt"
with open(path_methyl_markers, mode="r") as fr:
    list_methyl_markers = list(map(lambda x: x.rstrip("\n"), fr.readlines()))

list_file_dmps = glob.glob(f"{WORKDIR}/Results/11_dmp/Covariate_Sex_Age_Smoking_Medication/DMPExtract/Methylation_DMP_Extract.Control_Mild.Case_Severe.tsv", recursive=True)
for file_dmp in list_file_dmps:
    df_dmp = pd.read_csv(file_dmp, sep="\t")
    df_dmp["Marker"] = df_dmp["chr"].astype(str) + ":" + df_dmp["start"].astype(str)
    df_dmp_set_idx = df_dmp.set_index("Marker")
    df_dmp_marker_filt = df_dmp_set_idx.loc[list_methyl_markers, :]
    df_dmp_marker_filt
    break
# %%

# %%
import os
from pathlib import Path

import pandas as pd

# # %%
# WORKDIR = str(Path(os.path.abspath(__file__)).parents[3])
# path_methykittable = f"{WORKDIR}/Resources/Data/Methylation/MethylKitTable/MethylCpGTable.Control.Mild.Case.Severe.tsv"
# df_methylkittable = pd.read_csv(path_methykittable, sep="\t")

# # %%
# path_dmp_pos = f"{WORKDIR}/Results/11_dmp/Covariate_Sex_Age_CellType_Smoking_Medication_Comorbidity/DMPExtract/Methylation_DMP_Extract.Control_Mild.Case_Severe.filtered.fdr_05.sorted.tsv"
# df_dmp_pos = pd.read_csv(path_dmp_pos, sep="\t")
# df_dmp_pos["chrpos"] = df_dmp_pos["chr"] + ":" + df_dmp_pos["start"].astype(str)
# list_dmp_pos = df_dmp_pos["chrpos"].to_list()

# # %%
# dir_methylkittable_filtered = f"{WORKDIR}/Resources/Data/Methylation/MethylKitTable_Incl_EHR"
# os.makedirs(dir_methylkittable_filtered, exist_ok=True)
# outfilename = "MethylCpGTable.Control.Mild.Case.Severe.filtered.Covariate_Sex_Age_CellType_Smoking_Medication_Comorbidity.tsv"
# path_methykittable_filtered = os.path.join(dir_methylkittable_filtered, outfilename)

# %%
def filter_by_pos(df_methylkittable, list_dmp_pos, path_methykittable_filtered, save=True):
    if "chrpos" not in list(df_methylkittable.columns):
        df_methylkittable["chrpos"] = df_methylkittable["chr"] + ":" + df_methylkittable["start"].astype(str)
        df_methylkittable_set_idx = df_methylkittable.set_index("chrpos")
    
    df_methylkittable_filtered = df_methylkittable_set_idx.loc[list_dmp_pos, :]
    
    if save:
        df_methylkittable_filtered.to_csv(path_methykittable_filtered, sep="\t", index=False)
    
    return df_methylkittable_filtered
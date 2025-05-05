# %%
import numpy as np
import pandas as pd

# %%
path_confirmrecov = "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.DEG_by_Sex.Control_M.Case_F.Cov_Age.20250304.rawcount.tsv"
path_healthy = "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.HealthyControl.20250304.rawcount.tsv"
path_concat = "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.Confirmed.Recovered.HealthyControl.20250304.rawcount.tsv"


# %%
# path_confirmrecov = "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.DEG_by_Sex.Control_M.Case_F.Cov_Age.20250304.normcount.vst.tsv"
# path_healthy = "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.HealthyControl.20250304.normcount.vst.tsv"
# path_concat = "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.Confirmed.Recovered.HealthyControl.20250304.normcount.vst.tsv"

# %
df_confirmrecov = pd.read_csv(path_confirmrecov, sep="\t").reset_index(drop=False).rename(columns={"index": "ID"})
df_healthy = pd.read_csv(path_healthy, sep="\t").set_index("ID")
df_healthy = df_healthy.applymap(int).reset_index(drop=False)
df_concat = pd.merge(df_confirmrecov, df_healthy, how="inner", on="ID")
df_concat.to_csv(path_concat, sep="\t", index=False)

# %%

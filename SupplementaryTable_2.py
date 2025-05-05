# %%
import os
import pandas as pd

# %%
sumtable = "/BiO/Access/kyungwhan1998/Infectomics/Results/11_dmp/Yes_LOO/Covariate_Sex_Age/Summary_Table.DMP_Covariate_Sex_Age_Yes_LOO_common_across_9folds_or_more.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_7.20240326.tsv"

# %%
df_sumtable = pd.read_csv(sumtable, sep="\t")

# %%
list_key_gene = ["ENSG00000115607.10_IL18RAP", "ENSG00000176928.7_GCNT4", "ENSG00000118689.15_FOXO3", "ENSG00000140564.13_FURIN", "ENSG00000271303.2_SRXN1", "ENSG00000096060.15_FKBP5"]

# %%
df_sumtable = df_sumtable.set_index("RNA").loc[list_key_gene]
# %%
df_sumtable
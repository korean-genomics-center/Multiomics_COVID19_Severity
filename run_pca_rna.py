# %%
import gzip
import pickle

import pandas as pd
from run_pca_methyl_rna import run_pca
from sklearn.decomposition import PCA

path_rna = "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.Confirmed.Recovered.HealthyControl.20250304.batch_adj.platform_corrected.severity_preserved.normcount.vst.tsv"
path_meta = "/BiO/Access/kyungwhan1998/Infectomics/Resources/Data/RNA/COVID19_master_table_added_CRF_20250306.txt"
path_pca_obj = "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.Confirmed.Recovered.HealthyControl.20250304.batch_adj.platform_corrected.severity_preserved.normcount.vst.16genes.pc12.obj"
path_pca = "/BiO/Access/kyungwhan1998/Infectomics/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.Confirmed.Recovered.HealthyControl.20250304.batch_adj.platform_corrected.severity_preserved.normcount.vst.16genes.pc12.tsv"
path_marker = "/BiO/Access/kyungwhan1998/Infectomics/Results/11_dmp/Yes_LOO/Covariate_Sex_Age/list_new_gene_markers.txt"

with open(path_marker, mode="r") as fr:
    list_marker = list(set(list(map(lambda x: x.rstrip("\n"), fr.readlines()))))

table_rna = pd.read_csv(path_rna, sep="\t", index_col=[0])
# list_ensembl = list(map(lambda x: "_".join(x.split("_")[1:]), list(table_rna.index)))
# table_rna.index = list_ensembl
table_rna_select = table_rna.loc[list_marker]

table_transformed, pca = run_pca(table_rna_select.T, n_comp = 2)
table_transformed = table_transformed.reset_index(drop=False).rename(columns={"index": "Project_ID_Alias"})
table_meta = pd.read_csv(path_meta, sep="\t")
table_merged = pd.merge(table_transformed, table_meta, how="inner", on="Project_ID_Alias")
table_merged.to_csv(path_pca, sep="\t", index=False)
with gzip.open(path_pca_obj, mode="wb") as fb:
    pickle.dump(pca, fb)
# %%

# %%
import gzip
import os
import pickle
from pathlib import Path
import numpy as np
import pandas as pd

# %%
def save_list_markers(path_out_marker, list_markers):
    with open(path_out_marker, mode="w") as fw:
        for marker in list_markers:
            fw.write(marker + "\n")

# %%
WORKDIR = str(Path(__file__).parents[3])
# list_covariates = ["Covariate_Sex_Age_CCIScore", "Covariate_Sex_Age_CellType_2samples_CellPropTotal_Not_100_Included_New", "Covariate_Sex_Age_CellType_2samples_CellPropTotal_Not_100_Included_Smoking", "Covariate_Sex_Age_CellType_2samples_CellPropTotal_Not_100_Included_Smoking_CCI", "Covariate_Sex_Age_Smoking_Status", "Covariate_Sex_Age_Smoking_Status_CCIScore"]
list_covariates = ["Covariate_Sex_Age"]
for covariate in list_covariates:
    file_dict_dmp_loo_common = f"{WORKDIR}/Results/11_dmp/Yes_LOO/{covariate}/dictionary_dmp_loo_common_9folds.pk.gz"
    file_dict_deg_loo_common = f"{WORKDIR}/Results/11_dmp/Yes_LOO/DEG_LOO/dictionary_deg_loo_common_7folds.pk.gz"
    file_deg_dmp_overlap = f"{WORKDIR}/Results/11_dmp/Yes_LOO/{covariate}/Methyl_RNA_Correlation.Filtered.DMP_{covariate}_Yes_LOO_common_across_9folds_or_more.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_7.20240326.tsv"
    file_stats = f"{WORKDIR}/Results/11_dmp/Yes_LOO/{covariate}/Summary_Table.DMP_{covariate}_Yes_LOO_common_across_9folds_or_more.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_7.20240326.tsv"
    path_out_methyl_marker = f"{WORKDIR}/Results/11_dmp/Yes_LOO/{covariate}/list_new_methyl_markers.txt"
    path_out_gene_marker = f"{WORKDIR}/Results/11_dmp/Yes_LOO/{covariate}/list_new_gene_markers.txt"

    df_deg_dmp_overlap = pd.read_csv(file_deg_dmp_overlap, sep="\t")
    df_deg_dmp_overlap = df_deg_dmp_overlap[["Methyl", "RNA", "corr_rho", "corr_fdr"]]
    df_deg_dmp_overlap["Methyl"] = df_deg_dmp_overlap["Methyl"].apply(lambda x: ":".join(x.split("_")[:-1]), df_deg_dmp_overlap["Methyl"].to_list())
    list_methyl_marker = df_deg_dmp_overlap["Methyl"].to_list()
    save_list_markers(path_out_methyl_marker, list_methyl_marker)
    list_gene_marker = list(set(df_deg_dmp_overlap["RNA"].to_list()))
    save_list_markers(path_out_gene_marker, list_gene_marker)

    with gzip.open(file_dict_dmp_loo_common) as frb:
        dict_dmp_loo_common = pickle.load(frb)

    dict_median_beta_dmp_loo_common = {marker: np.median(beta) for marker, beta in dict_dmp_loo_common.items()}
    dict_median_beta_marker = {marker: median_beta for marker, median_beta in dict_median_beta_dmp_loo_common.items() if marker in list_methyl_marker}
    df_median_beta_marker = pd.DataFrame.from_dict(dict_median_beta_marker, orient="index").reset_index(drop=False)
    df_median_beta_marker.columns = ["Methyl", "Median_Beta"]
    df_merged_methyl = pd.merge(df_deg_dmp_overlap, df_median_beta_marker, how="inner", on="Methyl")

    with gzip.open(file_dict_deg_loo_common) as frb:
        dict_deg_loo_common = pickle.load(frb)

    dict_median_log2fc_deg_loo_common = {marker: np.median(log2fc) for marker, log2fc in dict_deg_loo_common.items()}
    dict_median_log2fc_marker = {marker: median_log2fc for marker, median_log2fc in dict_median_log2fc_deg_loo_common.items() if marker in list_gene_marker}
    df_median_log2fc_marker = pd.DataFrame.from_dict(dict_median_log2fc_marker, orient="index").reset_index(drop=False)
    df_median_log2fc_marker.columns = ["RNA", "Median_log2FC"]
    df_merged_methyl_rna = pd.merge(df_merged_methyl, df_median_log2fc_marker, how="inner", on="RNA")

    df_merged_methyl_rna.to_csv(file_stats, sep="\t", index=False)

# %%

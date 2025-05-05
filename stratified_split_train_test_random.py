# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from methylation_severity import (get_dict_deg_dmp_overlap_markers,
                                  get_drop_samples, main)
from sklearn.model_selection import train_test_split

# %%
path_sev_info = "/BiO/Access/kyungwhan1998/Infectomics/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"
dir_methylcpgmin = "/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/MethylCpGMin"
infilenamehyper = f"/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/Epigenetic_changes/first/hyper/dictionary_marker_freqC_all_samples_20240220.pk.gz"
infilenamehypo = f"/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/Epigenetic_changes/first/hypo/dictionary_marker_freqC_all_samples_20240220.pk.gz"
file_deg_dmp_overlap = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.DMP_DEG_LOO_annot.20240402.tsv"
path_drop_samples = "/BiO/Access/kyungwhan1998/Infectomics/Resources/Scripts/Final/list_remove_samples.txt"
outdir = "/BiO/Access/kyungwhan1998/Infectomics/Results/InfectomicsPaper1/20240906"
os.makedirs(outdir, exist_ok=True)
path_marker_target = "/BiO/Access/kyungwhan1998/Infectomics/Resources/Scripts/Final/list_methyl_markers.txt"
list_drop_samples = get_drop_samples(path_drop_samples)

# %%
def plot_stratified_random_split(data, train_y, test_y, col_stratify):
    # Visualize stratification
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))

    # Distribution in the original dataset
    original_counts = data[col_stratify].value_counts().sort_index()
    ax[0].bar(original_counts.index, original_counts.values,
              color='blue', alpha=0.7, label='Original Data')
    ax[0].set_title("Original Data Distribution")
    ax[0].set_xlabel("Stratify Column Value")
    ax[0].set_ylabel("Count")
    ax[0].legend()

    # Train and Test set distributions
    train_counts = train_y.value_counts().sort_index()
    test_counts = test_y.value_counts().sort_index()

    # Ensure consistent indices and align for plotting
    all_indices = sorted(set(train_counts.index).union(test_counts.index))
    train_counts = train_counts.reindex(all_indices, fill_value=0)
    test_counts = test_counts.reindex(all_indices, fill_value=0)

    bar_width = 0.35

    ax[1].bar([i - bar_width / 2 for i in range(len(all_indices))], train_counts.values,
              width=bar_width, label='Train Set', color='green', alpha=0.7)
    ax[1].bar([i + bar_width / 2 for i in range(len(all_indices))], test_counts.values,
              width=bar_width, label='Test Set', color='orange', alpha=0.7)
    ax[1].set_title("Train/Test Set Distribution")
    ax[1].set_xlabel("Stratify Column Value")
    ax[1].set_ylabel("Count")
    ax[1].set_xticks(range(len(all_indices)))
    ax[1].set_xticklabels(all_indices)
    ax[1].legend()

    plt.tight_layout()
    plt.show()

# %%
df_beta_all_hyper, df_beta_all_hypo = main(path_sev_info, dir_methylcpgmin, infilenamehyper, infilenamehypo)
dict_markers_overlap = get_dict_deg_dmp_overlap_markers(file_deg_dmp_overlap)
marker_overlap_hyper = list(set(dict_markers_overlap.keys()).intersection(set(list(df_beta_all_hyper.columns))))
marker_overlap_hypo = list(set(dict_markers_overlap.keys()).intersection(set(list(df_beta_all_hypo.columns))))
dict_markers_hyper = {k: v for k, v in dict_markers_overlap.items() if k in marker_overlap_hyper}
dict_markers_hypo = {k: v for k, v in dict_markers_overlap.items() if k in marker_overlap_hypo}

list_markers_hyper = list(set(dict_markers_hyper.keys()))
list_markers_hypo = list(set(dict_markers_hypo.keys()))
df_beta_selec_hyper = df_beta_all_hyper[list_markers_hyper]
df_beta_selec_hypo = df_beta_all_hypo[list_markers_hypo] 
list_metainfo = list(df_beta_all_hyper.iloc[:, -8:].columns)
df_metainfo = df_beta_all_hyper[list_metainfo]
df_beta_selec_all = pd.concat([df_beta_selec_hyper, df_beta_selec_hypo, df_metainfo], axis=1)

list_id = list(df_beta_selec_all.index)
list_keep_id = list(filter(lambda x: not str(x).endswith("V3"), list_id))
list_keep_id = list(filter(lambda x: not str(x).endswith("V4"), list_keep_id))
list_keep_id = list(filter(lambda x: str(x).startswith("C19-C"), list_keep_id))
df_beta_selec_target = df_beta_selec_all.loc[list_keep_id, :]

list_markers_all = list_markers_hyper + list_markers_hypo
X = df_beta_selec_target[list_markers_all]
y = df_beta_selec_target[["Severity_group"]]

train_X, test_X, _, _ = train_test_split(X, y, test_size=0.3, random_state=42, stratify=y)
# plot_stratified_random_split(y, train_y, test_y, col_stratify="Severity_group")
list_train_X = list(train_X.index)
list_test_X = list(test_X.index)
y.loc[list_train_X, "Dataset"] = "Train"
y.loc[list_test_X, "Dataset"] = "Test"

# %%
table_values = X.T.reset_index(drop=False).rename(columns={"index": "CpG"})
table_meta = y.reset_index(drop=False)

table_values.to_csv(f"{outdir}/Methylation_Beta_32_Markers_COVID19.txt", sep="\t", index=False)
table_meta.to_csv(f"{outdir}/Metadata_Train_Test_Split_COVID19.txt", sep="\t", index=False)
# %%

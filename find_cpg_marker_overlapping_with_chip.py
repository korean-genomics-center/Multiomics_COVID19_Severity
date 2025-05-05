# %%
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
df_beta_all_hyper, df_beta_all_hypo = main(path_sev_info, dir_methylcpgmin, infilenamehyper, infilenamehypo)
dict_markers_overlap = get_dict_deg_dmp_overlap_markers(file_deg_dmp_overlap)
marker_overlap_hyper = list(set(dict_markers_overlap.keys()).intersection(set(list(df_beta_all_hyper.columns))))
marker_overlap_hypo = list(set(dict_markers_overlap.keys()).intersection(set(list(df_beta_all_hypo.columns))))
dict_markers_hyper = {k: v for k, v in dict_markers_overlap.items() if k in marker_overlap_hyper}
dict_markers_hypo = {k: v for k, v in dict_markers_overlap.items() if k in marker_overlap_hypo}

list_markers_hyper = list(set(dict_markers_hyper.keys()))
list_markers_hypo = list(set(dict_markers_hypo.keys()))
list_markers = list_markers_hyper + list_markers_hypo
# %%
markers_chip = "/BiO/Access/kyungwhan1998/Infectomics/Resources/Data/MethylationArray_EPIC_EPICplus_EPICv2_HM450k_HM27k_GRCh38_filter.tsv"
markers_overlapping = "/BiO/Access/kyungwhan1998/Infectomics/Resources/Data/Overlapping_with_Chip_7_Markers.tsv"

with open(markers_chip, mode="r") as fr, open(markers_overlapping, mode="w") as fw:
    header = fr.readline()
    fw.write(header)
    for line in fr:
        record = line.rstrip("\n").split("\t")
        chr_pos = record[2] + ":" + record[3]
        if chr_pos in list_markers:
            overlap = "\t".join(record) + "\n"
            fw.write(overlap)
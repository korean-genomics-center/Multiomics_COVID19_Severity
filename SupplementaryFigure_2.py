# %%
import math
import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from methylation_severity import (get_dict_deg_dmp_overlap_markers,
                                  get_dict_palette, get_dict_pos,
                                  get_dictionary_methyl_pvalsig,
                                  get_drop_samples, get_xticklabels, main,
                                  make_dataframe_stat_test)
from statannot import add_stat_annotation

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
with open(path_marker_target, mode="r") as fr:
    list_excl_methyl = list(map(lambda x: x.rstrip(), fr.readlines()))
    
df_beta_all_hyper, df_beta_all_hypo = main(path_sev_info, dir_methylcpgmin, infilenamehyper, infilenamehypo)
dict_markers_overlap = get_dict_deg_dmp_overlap_markers(file_deg_dmp_overlap)
marker_overlap_hyper = list(set(dict_markers_overlap.keys()).intersection(set(list(df_beta_all_hyper.columns))))
marker_overlap_hypo = list(set(dict_markers_overlap.keys()).intersection(set(list(df_beta_all_hypo.columns))))
dict_markers_hyper = {k: v for k, v in dict_markers_overlap.items() if k in marker_overlap_hyper}
dict_markers_hypo = {k: v for k, v in dict_markers_overlap.items() if k in marker_overlap_hypo}

list_markers_hyper = list(set(dict_markers_hyper.keys()).difference(set(list_excl_methyl)))
list_markers_hypo = list(set(dict_markers_hypo.keys()).difference(set(list_excl_methyl)))
df_beta_selec_hyper = df_beta_all_hyper[list_markers_hyper]
df_beta_selec_hypo = df_beta_all_hypo[list_markers_hypo] 
list_metainfo = list(df_beta_all_hyper.iloc[:, -8:].columns)
df_metainfo = df_beta_all_hyper[list_metainfo]
df_beta_selec_all = pd.concat([df_beta_selec_hyper, df_beta_selec_hypo, df_metainfo], axis=1)

list_id = list(df_beta_selec_all.index)
list_keep_id = list(filter(lambda x: not str(x).endswith("V3"), list_id))
list_keep_id = list(filter(lambda x: not str(x).endswith("V4"), list_keep_id))
df_beta_selec_target = df_beta_selec_all.loc[list_keep_id, :]

list_target_methyl = list_markers_hyper + list_markers_hypo
dict_marker_target = dict()
for marker, markers in dict_markers_overlap.items():
    if marker in list_target_methyl:
        target = "-".join([x.split("_")[-1] for x in markers])
        dict_marker_target[marker] = target
        
df_marker_stat = make_dataframe_stat_test(df_beta_selec_target, list_target_methyl)
dict_methyl_pvalsig = get_dictionary_methyl_pvalsig(df_marker_stat)

dict_pos = get_dict_pos(df_beta_selec_target, list_pos=[0, 1, 2, 3])
dict_palette = get_dict_palette(df_beta_selec_target)

categories = list(dict_pos.keys())
positions = list(dict_pos.values())
colors = [dict_palette[cat] for cat in categories]
order = ['Healthy', 'Mild_First', 'Severe_First', 'Convalescent']

# %%
plt.rcParams["font.size"] = 10
colsev = "Severity_visit"
colsample = "Sample_ID"
nrows = 6
ncols = int(np.ceil(len(list_target_methyl) / nrows))
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols, 6*nrows))
axes = axes.flatten()

for i, (marker, ax) in enumerate(zip(list_target_methyl, axes)):
    df_beta_selec_target_copy = df_beta_selec_target.copy()
    df_marker_stat_copy = df_marker_stat.copy()
    plot_data = {cat: [] for cat in categories}
    
    for category in categories:
        values = df_beta_selec_target_copy[df_beta_selec_target_copy[colsev] == category][marker].dropna().values
        plot_data[category] = values

    sns.boxplot(data=df_beta_selec_target_copy, x=colsev, y=marker, order=order, palette=colors, ax=ax, fliersize=0, width=0.3, zorder=2)
    
    for category, color in zip(categories, colors):
        x = np.full_like(plot_data[category], dict_pos[category], dtype=float)
        ax.scatter(x, plot_data[category], color=color, edgecolor='black', label=category, alpha=0.5, s=30, zorder=3)
    
    box_pairs = df_marker_stat_copy.loc[marker, "comp"]
    list_idx_selec = [idx for idx, x in zip(range(len(box_pairs)), box_pairs) if x[0].split("_")[-1]==x[1].split("_")[-1] or x[1]!="Healthy"]
    box_pairs_selec = [box_pairs[idx] for idx in list_idx_selec]
    pvalues = df_marker_stat_copy.loc[marker, "padj"]
    pvalues_selec = [pvalues[idx] for idx in list_idx_selec]

    add_stat_annotation(ax, 
                        data=df_beta_selec_target_copy, 
                        x=colsev,
                        y=marker,
                        order=order,
                        box_pairs=box_pairs_selec,
                        perform_stat_test=False,
                        pvalues=pvalues_selec,
                        text_offset=-0.1,
                        line_offset=0.05,
                        text_format='star',
                        fontsize=plt.rcParams["font.size"]+2,
                        loc='inside', 
                        verbose=0)
    
    ax.set_xlabel("")
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xticks(positions)
    list_xticklabels = get_xticklabels(df_beta_selec_target, colsample=colsample)
    ax.set_xticklabels(list_xticklabels, fontsize=plt.rcParams["font.size"]+2)
    ax.tick_params(axis='y', labelsize=plt.rcParams["font.size"])
    ax.set_ylabel("Proportion of\nMethylated CpGs (%)", fontsize=plt.rcParams["font.size"]+2)
    markername = f"{marker}\n$\it({dict_marker_target.get(marker, marker)})$\n"
    legend_elements = [
        mpatches.Patch(color='darkgrey', label='Healthy controls'),
        mpatches.Patch(color='forestgreen', label='Mild (Acute)'),
        mpatches.Patch(color='firebrick', label='Severe (Acute)'),
        mpatches.Patch(color='royalblue', label='Convalescent')
    ]
    ax.set_title(markername, fontsize=plt.rcParams["font.size"]+4)
    axes[-1].legend(handles=legend_elements, 
            loc='center',
            bbox_to_anchor=(0.4, 0.4), 
            title="Groups", 
            fontsize=plt.rcParams["font.size"]+4, 
            title_fontsize=plt.rcParams["font.size"]+6,
            frameon=False)

    axes[-1].axis("off")
    for j in range(len(list_target_methyl), len(axes)):
        axes[j].axis('off')
    
    ax.grid(axis="y")
    ax.set_axisbelow(True)

plt.subplots_adjust(hspace=0.5, wspace=0.4, left=0.1, right=0.9)
plt.savefig("/BiO/Access/kyungwhan1998/Infectomics/Results/InfectomicsPaper1/20240906/SupplementaryFigure2.pdf", bbox_inches="tight", dpi=300)
plt.show()
plt.close()

# %%

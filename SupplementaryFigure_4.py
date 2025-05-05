# %%
import os
from pathlib import Path
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from methylation_severity import (get_dict_deg_dmp_overlap_markers,
                                  get_dict_palette, get_dict_pos,
                                  get_drop_samples, get_xticklabels, main,
                                  make_dataframe_stat_test)
from statannot import add_stat_annotation

# %%
WORKDIR = str(Path(__file__).parents[3])
path_sev_info = f"{WORKDIR}/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"
dir_methylcpgmin = f"{WORKDIR}/Results/10_methyl/MethylCpGMin"
infilenamehyper = f"{WORKDIR}/Results/10_methyl/Epigenetic_changes/first/hyper/dictionary_marker_freqC_all_samples_20240220.pk.gz"
infilenamehypo = f"{WORKDIR}/Results/10_methyl/Epigenetic_changes/first/hypo/dictionary_marker_freqC_all_samples_20240220.pk.gz"
file_deg_dmp_overlap = f"{WORKDIR}/Results/11_dmp/Yes_LOO/Covariate_Sex_Age/Methyl_RNA_Correlation.Filtered.DMP_Covariate_Sex_Age_Yes_LOO_common_across_9folds_or_more.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_7.20240326.tsv"
path_drop_samples = f"{WORKDIR}/Resources/Scripts/Final/list_remove_samples.txt"
outdir = f"{WORKDIR}/Results/Paper"
os.makedirs(outdir, exist_ok=True)
list_drop_samples = get_drop_samples(path_drop_samples)

# %%

list_excl_methyl = ["chr2:102143316", "chr5:75320838", "chr6:108562564", "chr15:90100633", "chr20:654151"] + ["chr6:36697843"]

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

list_target_methyl = list_markers_hyper + list_markers_hypo
list_target_methyl = sorted(list_target_methyl, key=lambda x: int(x.split(":")[0].replace("chr", "")))
dict_marker_target = dict()
for marker, markers in dict_markers_overlap.items():
    if marker in list_target_methyl:
        target = "-".join([x.split("_")[-1] for x in markers])
        dict_marker_target[marker] = target

list_id = list(df_beta_selec_all.index)
list_id = list(set(list_id).difference(set(list_drop_samples)))
df_beta_selec_target = df_beta_selec_all.loc[list_id, :]
df_methyl_stat = make_dataframe_stat_test(df_beta_selec_target, list_target_methyl)

# %%
plt.rcParams["font.size"] = 10
colsev = "Severity_visit"
colsample = "Sample_ID"
nrows = 5
ncols = int(np.ceil(len(list_target_methyl) / nrows))
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(3*ncols, 6*nrows))
axes = axes.flatten()

dict_pos = get_dict_pos(df_beta_selec_target, list_pos=[0, 1, 2, 3, 4, 5])
dict_palette = get_dict_palette(df_beta_selec_target)

categories = list(dict_pos.keys())
positions = list(dict_pos.values())
colors = [dict_palette[cat] for cat in categories]
order = ['Healthy','Mild_First','Severe_First','Mild_Last','Severe_Last','Convalescent']

for i, (marker, ax) in enumerate(zip(list_target_methyl, axes)):
    df_beta_selec_target_copy = df_beta_selec_target.copy()
    df_methyl_stat_copy = df_methyl_stat.copy()
    plot_data = {cat: [] for cat in categories}
    
    for category in categories:
        values = df_beta_selec_target_copy[df_beta_selec_target_copy[colsev] == category][marker].dropna().values
        plot_data[category] = values

    g = sns.boxplot(data=df_beta_selec_target_copy, x=colsev, y=marker, order=order, palette=colors, ax=ax, fliersize=0, width=0.3, zorder=2)
    
    g.set(xlabel=None) 
    for category, color in zip(categories, colors):
        x = np.full_like(plot_data[category], dict_pos[category], dtype=float)
        ax.scatter(x, plot_data[category], color=color, edgecolor='black', label=category, alpha=0.5, s=30, zorder=3)
    
    box_pairs = df_methyl_stat_copy.loc[marker, "comp"]
    list_idx_selec = [idx for idx, x in zip(range(len(box_pairs)), box_pairs) if x[0].split("_")[-1]==x[1].split("_")[-1] or x[0] == "Healthy"]
    box_pairs_selec = [box_pairs[idx] for idx in list_idx_selec]
    pvalues = df_methyl_stat_copy.loc[marker, "pval"]
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
                        fontsize=plt.rcParams["font.size"]+4,
                        loc='inside', 
                        verbose=0)
    
    ax.set_xlabel("")
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xticks(positions)
    list_xticklabels = get_xticklabels(df_beta_selec_target, dict_name_vis={"First": " (Acute)", "Last": " (Recovery)"}, colsample=colsample)

    if (i // ncols) == (nrows - 1):
        ax.set_xticklabels(list_xticklabels, fontsize=plt.rcParams["font.size"]+2, rotation=45, ha="right", rotation_mode="anchor")
    else:
        ax.set_xticklabels([], fontsize=0, rotation=45, ha="right", rotation_mode="anchor")
    
    ax.tick_params(axis='y', labelsize=plt.rcParams["font.size"]+2)
    
    if i % ncols == 0:
        ax.set_ylabel("Proportion of\nMethylated CpGs (%)", fontsize=plt.rcParams["font.size"]+4)
    else:
        ax.set_ylabel("", fontsize=0)
    markername = f"{marker}\n$\it({dict_marker_target.get(marker, marker)})$\n"
    ax.set_title(markername, fontsize=plt.rcParams["font.size"]+5)

    for j in range(len(list_target_methyl), len(axes)):
        axes[j].axis('off')
    
    ax.grid(axis="y")
    ax.set_axisbelow(True)

plt.subplots_adjust(hspace=0.4, wspace=0.5, left=0.1, right=0.9)
plt.savefig(f"{WORKDIR}/Results/Paper/SupplementaryFigure4.pdf", bbox_inches="tight", dpi=300)
plt.savefig(f"{WORKDIR}/Results/Paper/SupplementaryFigure4.png", bbox_inches="tight", dpi=300)
plt.show()
plt.close()

# %%

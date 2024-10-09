# %%
import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from statannot import add_stat_annotation

from rna_severity import (get_dict_palette, get_dict_pos,
                          get_dictionary_gene_pvalsig, get_xticklabels, main,
                          make_dataframe_stat_test_rna, save_stat_test_result)

# %%
col_id = "ID"
colsev = "Severity_visit"
path_drop_samples = "/BiO/Access/kyungwhan1998/Infectomics/Resources/Scripts/Final/list_remove_samples.txt"
path_exp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/DEGExtract/DEGExtract.RNA_samples_with_Methyl.DEG_by_Sex.Control_M.Case_F.Cov_Age.mincount_1.20240402.tsv.normcount"
path_methyl_marker = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.DMP_DEG_LOO_annot.20240402.tsv"
path_meta = "/BiO/Access/kyungwhan1998/Infectomics/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"
outdir = f"/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/rna"
os.makedirs(outdir, exist_ok=True)
df_rna_count_meta = main(path_drop_samples, path_exp, path_methyl_marker, path_meta)

list_targeted_gene = list(filter(lambda x: str(x).startswith("ENSG"), list(df_rna_count_meta.columns)))

dict_pos = get_dict_pos(df_rna_count_meta, colsev=colsev, list_pos=[0, 1, 2, 3, 4])
dict_palette = get_dict_palette(df_rna_count_meta, colsev=colsev)

categories = list(dict_pos.keys())
positions = list(dict_pos.values())
colors = [dict_palette[cat] for cat in categories]
list_sev = sorted(df_rna_count_meta[colsev].unique(), key=lambda x: x.split("_")[-1])
order = list_sev[1:] + [list_sev[0]]

df_gene_diffexp_sorted = make_dataframe_stat_test_rna(df_rna_count_meta, list_targeted_gene)
df_gene_diffexp_sorted_melted = save_stat_test_result(df_gene_diffexp_sorted, outdir, "SupplementaryTable5.txt", omics="rna")
dict_gene_pvalsig = get_dictionary_gene_pvalsig(df_gene_diffexp_sorted)

# %%
plt.rcParams["font.size"] = 10
nrows = 4
ncols = int(np.ceil(len(list_targeted_gene) / nrows))
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols, 6*nrows))
axes = axes.flatten()
for i, (gene, ax) in enumerate(zip(list_targeted_gene, axes)):
    df_rna_count_meta_copy = df_rna_count_meta.copy()
    df_gene_diffexp_sorted_copy = df_gene_diffexp_sorted.copy()
    plot_data = {cat: [] for cat in categories}
    
    for category in categories:
        values = df_rna_count_meta_copy[df_rna_count_meta_copy[colsev] == category][gene].dropna().values
        plot_data[category] = values

    sns.boxplot(data=df_rna_count_meta_copy, x=colsev, y=gene, palette=colors, order=order, width=0.3, fliersize=0, ax=ax, zorder=2)

    # Plot scatter plots
    for category, color in zip(categories, colors):
        x = np.full_like(plot_data[category], dict_pos[category], dtype=float)
        ax.scatter(x, plot_data[category], color=color, edgecolor='black', label=category, alpha=0.5, s=30, zorder=3)

    from itertools import combinations
    box_pairs = df_gene_diffexp_sorted_copy.loc[gene, "comp"]
    list_idx_selec = [idx for idx, x in zip(range(len(box_pairs)), box_pairs) if x[0].startswith("C") or x[1].startswith("C") or x[0].split("_")[-1] == x[1].split("_")[-1]]
    box_pairs_selec = [box_pairs[i] for i in list_idx_selec]
    pvalues = df_gene_diffexp_sorted_copy.loc[gene, "padj"]
    pvalues_selec = [pvalues[i] for i in list_idx_selec]
    
    add_stat_annotation(ax, 
                        data=df_rna_count_meta_copy, 
                        x=colsev,
                        y=gene,
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
    list_xticklabels = get_xticklabels(df_rna_count_meta_copy)
    ax.set_xticklabels(list_xticklabels, fontsize=plt.rcParams["font.size"]+2)
    ax.tick_params(axis='y', labelsize=plt.rcParams["font.size"]+2)
    ax.set_ylabel("Expression", fontsize=plt.rcParams["font.size"]+4)
    ax.set_title(f"$\it{gene.split('_')[1]}$", fontsize=plt.rcParams["font.size"]+6)
    
    legend_elements = [
        mpatches.Patch(color='darkgrey', label='Healthy controls'),
        mpatches.Patch(color='forestgreen', label='Mild'),
        mpatches.Patch(color='firebrick', label='Severe'),
        mpatches.Patch(color='royalblue', label='Convalescent')
    ]
    axes[-1].legend(handles=legend_elements, 
            loc='center',
            bbox_to_anchor=(0.4, 0.4), 
            title="Groups", 
            fontsize=plt.rcParams["font.size"]+4, 
            title_fontsize=plt.rcParams["font.size"]+6,
            frameon=False)

    axes[-1].axis("off")
    for j in range(len(list_targeted_gene), len(axes)):
        axes[j].axis('off')
    
    ax.grid(axis="y")
    ax.set_axisbelow(True)

plt.subplots_adjust(hspace=0.3, wspace=0.3, left=0.1, right=0.9)
plt.savefig("/BiO/Access/kyungwhan1998/Infectomics/Results/InfectomicsPaper1/20240906/SupplementaryFigure5.png", bbox_inches="tight", dpi=600)
plt.show()
plt.close()

# %%

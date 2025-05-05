# %%
import os
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from rna_severity import (get_dict_palette,
                          get_dict_pos, 
                          get_xticklabels, 
                          main,
                          make_dataframe_stat_test_rna)
from statannot import add_stat_annotation

# %%
WORKDIR = str(Path(__file__).parents[3])
col_id = "Sample_ID"
colsev = "Severity_Binary"
path_drop_samples = f"{WORKDIR}/Resources/Scripts/Final/list_remove_samples.txt"
path_exp = f"{WORKDIR}/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.Confirmed.Recovered.HealthyControl.20250304.batch_adj.platform_corrected.severity_preserved.normcount.vst.tsv"
file_deg_dmp_overlap = f"{WORKDIR}/Results/11_dmp/Yes_LOO/Covariate_Sex_Age/Methyl_RNA_Correlation.Filtered.DMP_Covariate_Sex_Age_Yes_LOO_common_across_9folds_or_more.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_7.20240326.tsv"
path_meta = f"{WORKDIR}/Resources/Data/RNA/COVID19_master_table_added_CRF_20231007.txt"
outdir = f"{WORKDIR}/Results/Paper"
os.makedirs(outdir, exist_ok=True)
list_marker_excl = ['ENSG00000115607.10_IL18RAP',  'ENSG00000176928.7_GCNT4', 'ENSG00000118689.15_FOXO3',  'ENSG00000140564.13_FURIN', 'ENSG00000271303.2_SRXN1', 'ENSG00000096060.15_FKBP5']

# %%
df_rna_count_meta = main(path_drop_samples, path_exp, file_deg_dmp_overlap, path_meta)
df_rna_count_meta = df_rna_count_meta[~df_rna_count_meta[col_id].str.contains("V3")]
df_rna_count_meta = df_rna_count_meta[~df_rna_count_meta[col_id].str.contains("V4")]
list_genes = list(filter(lambda x: str(x).startswith("ENSG"), list(df_rna_count_meta.columns)))
list_targeted_gene = list(set(list_genes).difference(set(list_marker_excl)))

dict_pos = get_dict_pos(df_rna_count_meta, colsev=colsev, list_pos=[0, 1, 2, 3])
dict_palette = get_dict_palette(df_rna_count_meta, colsev=colsev)

categories = list(dict_pos.keys())
positions = list(dict_pos.values())
colors = [dict_palette[cat] for cat in categories]
list_sev = sorted(df_rna_count_meta[colsev].unique())
order = list_sev[1:] + [list_sev[0]]

df_gene_diffexp_sorted = make_dataframe_stat_test_rna(df_rna_count_meta, list_targeted_gene)
# %%
plt.rcParams["font.size"] = 10
nrows = 3
ncols = int(np.ceil(len(list_targeted_gene) / nrows))
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(3*ncols, 5*nrows))
axes = axes.flatten()

for i, (gene, ax) in enumerate(zip(list_targeted_gene, axes)):
    df_rna_count_meta_copy = df_rna_count_meta.copy()
    df_gene_diffexp_sorted_copy = df_gene_diffexp_sorted.copy()
    plot_data = {cat: [] for cat in categories}
    
    for category in categories:
        values = df_rna_count_meta_copy[df_rna_count_meta_copy[colsev] == category][gene].dropna().values
        plot_data[category] = values

    sns.boxplot(data=df_rna_count_meta_copy, x=colsev, y=gene, palette=colors, order=order, width=0.4, fliersize=0, ax=ax, zorder=3)

    for category, color in zip(categories, colors):
        x = np.full_like(plot_data[category], dict_pos[category], dtype=float)
        ax.scatter(x, plot_data[category], color=color, edgecolor='black', label=category, alpha=0.5, s=30, zorder=5)

    box_pairs = df_gene_diffexp_sorted_copy.loc[gene, "comp"]
    pvalues = df_gene_diffexp_sorted_copy.loc[gene, "pval"]
    
    add_stat_annotation(ax, 
                        data=df_rna_count_meta_copy, 
                        x=colsev,
                        y=gene,
                        order=order,
                        box_pairs=box_pairs,
                        perform_stat_test=False,
                        pvalues=pvalues,
                        text_offset=-0.1,
                        line_offset=0.05,
                        text_format='star',
                        fontsize=plt.rcParams["font.size"],
                        loc='inside', 
                        verbose=0)
    ax.set_xlabel("")
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xticks(positions)
    list_xticklabels = get_xticklabels(df_rna_count_meta_copy)
    if i // ncols == nrows - 1:
        ax.set_xticklabels(list_xticklabels, fontsize=plt.rcParams["font.size"]+2, rotation=45, ha="right", rotation_mode="anchor")
    else:
        ax.set_xticklabels([], fontsize=0, rotation=45, ha="right", rotation_mode="anchor")
    
    ax.tick_params(axis='y', labelsize=plt.rcParams["font.size"]+2)
    row_starts = [row * ncols for row in range(nrows) if row * ncols < len(list_targeted_gene)]

    if i in row_starts:
        ax.set_ylabel("VST count", fontsize=plt.rcParams["font.size"]+4)
    else:
        ax.set_ylabel("", fontsize=0)
    ax.set_title(f"$\it{gene.split('_')[1]}$", fontsize=plt.rcParams["font.size"]+6)
    ax.grid(axis="y")
    ax.set_axisbelow(True)

legend_elements = [
        mpatches.Patch(color='darkgrey', label='Healthy'),
        mpatches.Patch(color='forestgreen', label='Mild-Moderate (Acute)'),
        mpatches.Patch(color='firebrick', label='Severe-Critical (Acute)'),
        mpatches.Patch(color='royalblue', label='Convalescent')
]

axes[-1].legend(handles=legend_elements, 
           loc='center',
           bbox_to_anchor=(0.4, 0.4),
           title="Groups", 
           fontsize=plt.rcParams["font.size"]+4, 
           title_fontsize=plt.rcParams["font.size"]+6,
           frameon=False)

for j in range(len(list_targeted_gene), len(axes)):
    axes[j].axis('off')

plt.subplots_adjust(hspace=0.2, wspace=0.6, left=0.1, right=0.9)
plt.savefig(f"{WORKDIR}/Results/Paper/SupplementaryFigure4.pdf", bbox_inches="tight", dpi=300)
plt.savefig(f"{WORKDIR}/Results/Paper/SupplementaryFigure4.png", bbox_inches="tight", dpi=300)
plt.show()
plt.close()

# %%

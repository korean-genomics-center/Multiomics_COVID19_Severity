# %%
# %%
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from rna_severity import (
    get_dict_palette,
    get_dict_pos,
    get_xticklabels,
    main,
    make_dataframe_stat_test_rna
)
from statannot import add_stat_annotation

# %%
WORKDIR = str(Path(__file__).parents[3])
col_id = "Sample_ID"
colsev = "Severity_Binary"
colvis = "Visit_order"
dict_vis = {"Visit1": "First", "Visit3": "Last", "Visit4": "Last"}
path_drop_samples = f"{WORKDIR}/Resources/Scripts/Final/list_remove_samples.txt"
path_exp = f"{WORKDIR}/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.Confirmed.Recovered.HealthyControl.20250304.batch_adj.platform_corrected.severity_preserved.normcount.vst.tsv"
file_deg_dmp_overlap = f"{WORKDIR}/Results/11_dmp/Yes_LOO/Covariate_Sex_Age/Methyl_RNA_Correlation.Filtered.DMP_Covariate_Sex_Age_Yes_LOO_common_across_9folds_or_more.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_7.20240326.tsv"
path_meta = f"{WORKDIR}/Resources/Data/RNA/COVID19_master_table_added_CRF_20231007.txt"
outdir = f"{WORKDIR}/Results/Paper"
os.makedirs(outdir, exist_ok=True)

df_rna_count_meta = main(path_drop_samples, path_exp, file_deg_dmp_overlap, path_meta)
df_rna_count_meta["Visit"] = df_rna_count_meta[colvis].apply(lambda x: dict_vis.get(x, x))
df_rna_count_meta["Severity_visit"] = df_rna_count_meta.apply(
    lambda row: row[colsev] if row[colsev] in ["Healthy", "Convalescent"]
    else f"{row[colsev]}_{row['Visit']}",
    axis=1
)

list_targeted_gene = list(filter(lambda x: str(x).startswith("ENSG"), list(df_rna_count_meta.columns)))
dict_palette = get_dict_palette(df_rna_count_meta, colsev="Severity_visit")

order = ['Healthy', 'Mild_First', 'Severe_First', 'Mild_Last', 'Severe_Last', 'Convalescent']
colors = [dict_palette[x] for x in order]
dict_pos = {x: i for i, x in enumerate(order)}

df_gene_diffexp_sorted = make_dataframe_stat_test_rna(df_rna_count_meta, list_targeted_gene, colsev="Severity_visit")

# %%
plt.rcParams["font.size"] = 10
nrows = 4
ncols = int(np.ceil(len(list_targeted_gene) / nrows))
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(3*ncols, 5*nrows))
axes = axes.flatten()

for i, (gene, ax) in enumerate(zip(list_targeted_gene, axes)):
    df_rna_count_meta_copy = df_rna_count_meta.copy()
    df_gene_diffexp_sorted_copy = df_gene_diffexp_sorted.copy()

    plot_data = {cat: df_rna_count_meta_copy[df_rna_count_meta_copy["Severity_visit"] == cat][gene].dropna().values
                 for cat in order}

    sns.boxplot(data=df_rna_count_meta_copy, x="Severity_visit", y=gene,
                palette=dict_palette, order=order, width=0.5, fliersize=0, ax=ax, zorder=2)

    for category in order:
        values = plot_data[category]
        x = np.full_like(values, order.index(category), dtype=float)
        ax.scatter(x, values, color=dict_palette[category], edgecolor='black', alpha=0.5, s=30, zorder=3)

    box_pairs = df_gene_diffexp_sorted_copy.loc[gene, "comp"]
    pvalues = df_gene_diffexp_sorted_copy.loc[gene, "pval"]
    list_idx_selec = [
        idx for idx, (a, b) in enumerate(box_pairs)
        if a.startswith("C") or b.startswith("C") or a.split("_")[-1] == b.split("_")[-1]
    ]
    box_pairs_selec = [box_pairs[i] for i in list_idx_selec]
    pvalues_selec = [pvalues[i] for i in list_idx_selec]

    add_stat_annotation(ax,
                        data=df_rna_count_meta_copy,
                        x="Severity_visit", y=gene,
                        order=order, box_pairs=box_pairs_selec,
                        perform_stat_test=False,
                        pvalues=pvalues_selec,
                        text_offset=-0.1,
                        line_offset=0.05,
                        text_format='star',
                        fontsize=plt.rcParams["font.size"] + 2,
                        loc='inside',
                        verbose=0)

    ax.set_xlabel("")
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xticks(range(len(order)))
    list_xticklabels = get_xticklabels(df_rna_count_meta_copy,
                                       dict_name_vis={"First": " (Acute)", "Last": " (Recovery)"},
                                       colsample="Sample_ID", colsev="Severity_visit")
    if i // ncols == nrows - 1:
        ax.set_xticklabels(list_xticklabels, fontsize=plt.rcParams["font.size"]+2, rotation=45, ha="right", rotation_mode="anchor")
    else:
        ax.set_xticklabels([], fontsize=0)

    ax.tick_params(axis='y', labelsize=plt.rcParams["font.size"] + 2)
    if i % ncols == 0:
        ax.set_ylabel("VST count", fontsize=plt.rcParams["font.size"] + 4)
    else:
        ax.set_ylabel("")
    ax.set_title(f"$\it{gene.split('_')[1]}$", fontsize=plt.rcParams["font.size"] + 6)
    ax.grid(axis="y")
    ax.set_axisbelow(True)

plt.subplots_adjust(hspace=0.2, wspace=0.6, left=0.1, right=0.9)
plt.savefig(f"{outdir}/SupplementaryFigure5.png", bbox_inches="tight", dpi=300)
plt.savefig(f"{outdir}/SupplementaryFigure5.pdf", bbox_inches="tight", dpi=300)
plt.show()
plt.close()
# %%
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import os
import pandas as pd
from pathlib import Path
import matplotlib.patches as mpatches

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
path_methyl_target = f"{WORKDIR}/Results/11_dmp/Yes_LOO/Covariate_Sex_Age/list_new_methyl_markers.txt"

# %%
from methylation_severity import (get_dict_deg_dmp_overlap_markers,
                                  main,
                                  make_dataframe_normality_test,
                                  make_dataframe_stat_test,
                                  save_norm_test_result, 
                                  save_stat_test_result,
                                  plot_methyl_difference)

df_beta_all_hyper, df_beta_all_hypo = main(path_sev_info, dir_methylcpgmin, infilenamehyper, infilenamehypo)
list_metainfo = list(df_beta_all_hyper.iloc[:, -8:].columns)
list_marker_hyper = list(filter(lambda x: x not in list_metainfo, list(df_beta_all_hyper.columns)))
df_beta_all_hyper_hypo = pd.concat([df_beta_all_hyper[list_marker_hyper], df_beta_all_hypo], axis=1)
dict_markers_overlap = get_dict_deg_dmp_overlap_markers(file_deg_dmp_overlap)

marker_overlap_hyper = list(set(dict_markers_overlap.keys()).intersection(set(list(df_beta_all_hyper.columns))))
marker_overlap_hypo = list(set(dict_markers_overlap.keys()).intersection(set(list(df_beta_all_hypo.columns))))

with open(path_methyl_target, mode="r") as fr:
    list_target_methyl = list(map(lambda x: x.rstrip(), fr.readlines()))

with open(path_drop_samples, mode="r") as fr:
    list_drop_samples = list(map(lambda x: x.rstrip(), fr.readlines()))

dict_marker_target = dict()
for marker, markers in dict_markers_overlap.items():
    if marker in list_target_methyl:
        target = "-".join([x.split("_")[-1] for x in markers])
        dict_marker_target[marker] = target
     
df_beta_selec_target = df_beta_all_hyper_hypo[list_target_methyl + list_metainfo]
df_beta_selec_target = df_beta_selec_target[~(df_beta_selec_target.index.str.contains("V3"))]
df_beta_selec_target = df_beta_selec_target[~(df_beta_selec_target.index.str.contains("V4"))]
list_samples_filtered_in = list(filter(lambda x: x not in list_drop_samples, list(df_beta_selec_target.index)))
df_beta_selec_target = df_beta_selec_target.loc[list_samples_filtered_in]

df_methyl_marker_norm = make_dataframe_normality_test(df_beta_selec_target, list_target_methyl)
save_norm_test_result(df_methyl_marker_norm, outdir, "methylation_meanbeta_normality_within_severity_groups.txt")
df_beta_methyl_all = make_dataframe_stat_test(df_beta_selec_target, list_target_methyl)
df_stat_methyl = save_stat_test_result(df_beta_methyl_all, outdir, "methylation_meanbeta_differences_between_severity_groups.txt", dict_marker_target)

# %%
from rna_severity import (main, 
                          make_dataframe_normality_test_rna,
                          make_dataframe_stat_test_rna,
                          save_norm_test_result_rna, 
                          save_stat_test_result_rna,
                          plot_rna_difference)

col_id = "Sample_ID"
colsev = "Severity_Binary"
path_drop_samples = f"{WORKDIR}/Resources/Scripts/Final/list_remove_samples.txt"
path_exp = f"{WORKDIR}/Results/5_deg/DEGExtract/DEGExtract.RNA_samples_with_Methyl.Confirmed.Recovered.HealthyControl.20250304.batch_adj.platform_corrected.severity_preserved.normcount.vst.tsv"
file_deg_dmp_overlap = f"{WORKDIR}/Results/11_dmp/Yes_LOO/Covariate_Sex_Age/Methyl_RNA_Correlation.Filtered.DMP_Covariate_Sex_Age_Yes_LOO_common_across_9folds_or_more.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_7.20240326.tsv"
path_meta = f"{WORKDIR}/Resources/Data/RNA/COVID19_master_table_added_CRF_20231007.txt"
path_gene_target = f"{WORKDIR}/Results/11_dmp/Yes_LOO/Covariate_Sex_Age/list_new_gene_markers.txt"
outdir = f"{WORKDIR}/Results/Paper"
os.makedirs(outdir, exist_ok=True)

with open(path_gene_target, mode="r") as fr:
    list_target_gene = list(map(lambda x: x.rstrip(), fr.readlines()))

df_rna_count_meta = main(path_drop_samples, path_exp, file_deg_dmp_overlap, path_meta)
df_rna_count_meta = df_rna_count_meta[~df_rna_count_meta[col_id].str.contains("V3")]
df_rna_count_meta = df_rna_count_meta[~df_rna_count_meta[col_id].str.contains("V4")]

df_rna_marker_norm = make_dataframe_normality_test_rna(df_rna_count_meta, list_target_gene, colsev=colsev)
save_norm_test_result_rna(df_rna_marker_norm, outdir, "rna_meanexp_normality_within_severity_groups.txt")
df_marker_diffexp_all = make_dataframe_stat_test_rna(df_rna_count_meta, list_target_gene)
df_stat_rna = save_stat_test_result_rna(df_marker_diffexp_all, outdir, "rna_meanexp_differences_between_severity_groups.txt")

# %%
list_key_methyl = ["chr2:102143316", "chr5:75320838", "chr6:108562564", "chr6:36697843", "chr15:90100633", "chr20:654151"]
df_stat_key_methyl = save_stat_test_result(df_beta_methyl_all, outdir, "methylation_meanbeta_differences_between_severity_groups_key_methyl_only.txt", dict_marker_target, list_key_methyl)

list_key_gene = ["ENSG00000115607.10_IL18RAP", "ENSG00000176928.7_GCNT4", "ENSG00000118689.15_FOXO3", "ENSG00000140564.13_FURIN", "ENSG00000271303.2_SRXN1", "ENSG00000096060.15_FKBP5"]
df_stat_key_rna = save_stat_test_result_rna(df_marker_diffexp_all, outdir, "rna_meanexp_differences_between_severity_groups_key_genes_only.txt", list_key_gene)

# %%
plt.rcParams["font.size"] = 10
figsize = (12, 12)
plot_linewidth = 1

fig = plt.figure(figsize = figsize)
col = 100
row = 100
gsfig = gridspec.GridSpec(row, col, left = 0, right = 1, bottom = 0, top = 1, wspace = 1, hspace = 1)

gs1 = gsfig[0:45, 0:80]
ax1 = fig.add_subplot(gs1)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.axis("off")
ax1.text(-0.05, 1.15, "A", transform=ax1.transAxes,
        fontsize=plt.rcParams["font.size"]+10, fontweight='bold', va='top', ha='right')

gs2 = gsfig[0:45, 90:100]
ax2 = fig.add_subplot(gs2)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.axis("off")
ax2.text(-0.45, 1.15, "C", transform=ax2.transAxes,
        fontsize=plt.rcParams["font.size"]+10, fontweight='bold', va='top', ha='right')

gs3 = gsfig[60:100, 0:80]
ax3 = fig.add_subplot(gs3)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.axis("off")
ax3.text(-0.05, 1.10, "B", transform=ax3.transAxes,
        fontsize=plt.rcParams["font.size"]+10, fontweight='bold', va='top', ha='right')

gs4 = gsfig[60:100, 90:100]
ax4 = fig.add_subplot(gs4)
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.spines['left'].set_visible(False)
ax4.spines['bottom'].set_visible(False)
ax4.axis("off")
ax4.text(-0.45, 1.10, "D", transform=ax4.transAxes,
        fontsize=plt.rcParams["font.size"]+10, fontweight='bold', va='top', ha='right')

legend_elements = [
        mpatches.Patch(color='darkgrey', label='Healthy'),
        mpatches.Patch(color='forestgreen', label='Mild-Moderate (Acute)'),
        mpatches.Patch(color='firebrick', label='Severe-Critical (Acute)'),
        mpatches.Patch(color='royalblue', label='Convalescent')
]

list_plot_methyl_celltype_noadj = ["chr2:102143316", "chr5:75320838", "chr6:108562564", "chr15:90100633", "chr20:654151"]
list_plot_methyl_celltype_adj = ["chr6:36697843"]
plot_methyl_difference(df_beta_selec_target, df_beta_methyl_all, list_plot_methyl_celltype_noadj, dict_marker_target, gs1, fig, colsev="Severity_visit", pval="pval", nrows=1, list_pos=[0, 1, 2, 3])
axes_methyl = plot_methyl_difference(df_beta_selec_target, df_beta_methyl_all, list_plot_methyl_celltype_adj, dict_marker_target, gs2, fig, colsev="Severity_visit", pval="pval", nrows=1, list_pos=[0, 1, 2, 3])

axes_methyl[-1].legend(handles=legend_elements, 
        loc='center',
        bbox_to_anchor=(2.3, -0.12), 
        title="Groups", 
        fontsize=plt.rcParams["font.size"]+1, 
        title_fontsize=plt.rcParams["font.size"]+2,
        frameon=False)

list_plot_gene_celltype_noadj = ["ENSG00000115607.10_IL18RAP", "ENSG00000176928.7_GCNT4", "ENSG00000118689.15_FOXO3", "ENSG00000140564.13_FURIN", "ENSG00000271303.2_SRXN1"]
list_plot_gene_celltype_adj = ["ENSG00000096060.15_FKBP5"]
plot_rna_difference(df_rna_count_meta, df_marker_diffexp_all, list_plot_gene_celltype_noadj, gs3, fig, colsev="Severity_Binary", pval="pval", nrows=1, list_pos=[0, 1, 2, 3])
axes_rna = plot_rna_difference(df_rna_count_meta, df_marker_diffexp_all, list_plot_gene_celltype_adj, gs4, fig, colsev="Severity_Binary", pval="pval", nrows=1, list_pos=[0, 1, 2, 3])

# plt.savefig(f"{outdir}/Figure2.pdf", bbox_inches="tight", dpi=300)
# plt.savefig(f"{outdir}/Figure2.png", bbox_inches="tight", dpi=300)
plt.show()
plt.close()

# %%

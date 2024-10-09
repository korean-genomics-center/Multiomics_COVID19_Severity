# %%
import glob
import gzip
import math
import os
import pickle
from itertools import combinations

import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ranksums, ttest_ind
from statannot import add_stat_annotation
from statsmodels.stats.multitest import fdrcorrection

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

# %%
def get_drop_samples(path_drop_samples):
    with open(path_drop_samples, mode="r") as frmv:
        list_drop_samples = [x.rstrip() for x in frmv.readlines()]
    
    return list_drop_samples

def get_list_files_methylcpgmin(dir_methylcpgmin):
    list_files_methylcpgmin = glob.glob(f"{dir_methylcpgmin}/**/*pair_merged.methyl_cpg_min.tsv", recursive=True)

    return list_files_methylcpgmin

def get_list_methyl_sample(list_files_methylcpgmin):
    list_methyl_sample = list()
    for file_methylcpgmin in list_files_methylcpgmin:
        dir_methylcpgmin = os.path.basename(os.path.dirname(file_methylcpgmin))
        if dir_methylcpgmin == "HealthyControl":
            name_sample = os.path.basename(file_methylcpgmin).split(".")[0]
        else:
            name_sample = os.path.basename(dir_methylcpgmin)
        list_methyl_sample.append(name_sample)

    return list_methyl_sample

def load_pickle(loadfilename):
    with gzip.open(loadfilename,'rb') as fr:
        data = pickle.load(fr)
    
    return data

def get_dict_deg_dmp_overlap_markers(file_deg_dmp_overlap, col_methyl="Methyl", col_gene = "RNA"):
    from collections import defaultdict
    df_deg_dmp_overlap = pd.read_csv(file_deg_dmp_overlap, sep="\t")
    list_methyl = df_deg_dmp_overlap[col_methyl].to_list()
    list_gene = df_deg_dmp_overlap[col_gene].to_list()
    dict_all_markers = defaultdict(list)
    for marker, markersym in zip(list_methyl, list_gene):
        marker_edit = ":".join(marker.split("_")[:-1])
        dict_all_markers[marker_edit].append(markersym)

    return dict_all_markers

def get_dataframe_methylation_beta_samples(infilename, list_selected_samples):
    dict_cpgmin_all_samples = load_pickle(infilename)     
    dict_cpgmin_all_samples_marker_fixed = dict()
    for sampleid, dict_cpgs in dict_cpgmin_all_samples.items():
        dict_cpgs_marker_fixed = dict()
        for marker, freqC in dict_cpgs.items():
            fixed_marker = ":".join(marker)
            dict_cpgs_marker_fixed[fixed_marker] = freqC
        dict_cpgmin_all_samples_marker_fixed[sampleid] = dict_cpgs_marker_fixed

    df_beta = pd.DataFrame.from_dict(dict_cpgmin_all_samples_marker_fixed).astype(float)
    df_beta_transposed = df_beta.T.reset_index(drop=False).rename(columns={"index": "Sample_ID"})
    df_sev = pd.read_csv(path_sev_info, sep="\t")
    df_beta_sev = pd.merge(df_beta_transposed, df_sev, how="inner", on="Sample_ID")
    df_beta_sev["Visit"] = df_beta_sev["Visit"].fillna("Healthy")
    df_beta_sev["Visit_order"] = df_beta_sev["Visit_order"].fillna("Visit0")
    df_beta_sev["Severity_group"] = df_beta_sev["Severity_group"].fillna("Healthy")
    df_beta_sev["Subject_ID"] = df_beta_sev["Subject_ID"].fillna(df_beta_sev["Sample_ID"])
    df_beta_sev = df_beta_sev[df_beta_sev["Visit_order"]!="Visit6"]
    df_beta_sev["Severity_visit"] = df_beta_sev["Severity_visit"].apply(lambda x: "Healthy" if "Healthy" in x else x)
    df_beta_sev["Severity_visit"] = df_beta_sev["Severity_visit"].apply(lambda x: "Convalescent" if "Convalescent" in x else x)
    df_beta_set_idx = df_beta_sev.set_index("Sample_ID")
    # df_beta_set_idx_select = df_beta_set_idx.loc[list_selected_samples, :]
    df_beta_set_idx_select = df_beta_set_idx[df_beta_set_idx.index.isin(list_selected_samples)]
    
    return df_beta_set_idx_select

def find_meandiff(a, b):
    meandiff = np.mean(a) - np.mean(b)

    return meandiff

def make_dataframe_stat_test(df_marker, list_target_marker, colsev="Severity_visit"):
    dict_marker_stat = dict()
    for marker in list_target_marker:
        df_group_visit = df_marker.groupby(colsev)[marker].apply(np.array).reset_index(drop=False)
        list_comb_group = list(combinations(df_group_visit[colsev].values, r=2))
        list_comb_group = list(map(lambda x: x[::-1], list_comb_group))
        list_comb_marker_val = list(combinations(df_group_visit[marker].values, r=2))
        list_comb_marker_val = list(map(lambda x: x[::-1], list_comb_marker_val))
        list_meandiff = list(map(lambda x: find_meandiff(*x), list_comb_marker_val))
        list_stat = list(map(lambda x: ttest_ind(*x, equal_var=False)[0], list_comb_marker_val))
        list_pval = list(map(lambda x: ttest_ind(*x, equal_var=False)[1], list_comb_marker_val))
        list_sig, list_padj = fdrcorrection(list_pval)
        dict_stat_res = {"comp": list_comb_group, "diff": list_meandiff, "stat": list_stat, "pval": list_pval, "padj": list_padj, "is_sig": list_sig}
        dict_marker_stat[marker] = dict_stat_res
                            
    df_marker_stat = pd.DataFrame.from_dict(dict_marker_stat).T
    df_marker_stat_sorted = df_marker_stat.sort_values(by=["pval"], ascending=True)
    
    return df_marker_stat_sorted

def save_stat_test_result(df_marker_stat_sorted, outdir, filename, dict_marker_target):
    df_marker_stat_sorted_melted = pd.DataFrame(columns = ["CpG", "Gene", "Comp1", "Comp2", "diff", "stat", "pval", "padj", "is_sig"])

    ind_df = 0
    for cpgname, row in df_marker_stat_sorted.iterrows():
        comp_list = row["comp"]
        diff_list = row["diff"]
        stat_list = row["stat"]
        pval_list = row["pval"]
        padj_list = row["padj"]
        is_sig_list = row["is_sig"]
        
        for comp, diff, stat, pval, padj, is_sig in zip(comp_list, diff_list, stat_list, pval_list, padj_list, is_sig_list):
            comp1, comp2 = comp
            genename = dict_marker_target[cpgname]
            values_row = [cpgname, genename, comp1, comp2, diff, stat, pval, padj, is_sig]
            df_marker_stat_sorted_melted.loc[ind_df, :] = values_row
            ind_df += 1
    
    outfilepath = os.path.join(outdir, filename)
    first_col = list(df_marker_stat_sorted_melted.columns)[0]
    df_marker_stat_sorted_melted = df_marker_stat_sorted_melted.sort_values(by=[first_col], ascending=False)
    df_marker_stat_sorted_melted.to_csv(outfilepath, sep="\t", index=False)
    
    return df_marker_stat_sorted_melted


def get_sorted_severity(df_beta_methyl, 
                        dict_order_visit = {"Healthy": 0, "First": 1, "Last": 2, "Convalescent": 3},
                        dict_order_group = {"Healthy": 0, "Mild": 1, "Severe": 2, "Convalescent": 3},
                        colsev="Severity_visit"):
    list_severity = list(df_beta_methyl[colsev].unique())
    list_severity_sorted = sorted(list_severity, key=lambda x: dict_order_group[x.split("_")[0]])
    list_severity_sorted = sorted(list_severity_sorted, key=lambda x: dict_order_visit[x.split("_")[-1]])
    
    return list_severity_sorted

def get_dict_palette(df_beta_methyl, colsev="Severity_visit"):
    list_severity_sorted = get_sorted_severity(df_beta_methyl, colsev=colsev)
    
    list_colors = list()
    for sev_vis in list_severity_sorted:
        if sev_vis.split("_")[0]=="Mild":
            list_colors.append("forestgreen")
        elif sev_vis.split("_")[0]=="Severe":
            list_colors.append("firebrick") 
        elif sev_vis.split("_")[0]=="Healthy":
            list_colors.append("darkgrey")
        else: 
            list_colors.append("royalblue")
            
    dict_palette = dict(zip(list_severity_sorted, list_colors))

    return dict_palette

def get_dict_pos(df_beta_methyl, 
                 list_pos=[0, 1, 1.5, 2.5, 3, 4], 
                 colsev="Severity_visit"):
    list_severity_sorted = get_sorted_severity(df_beta_methyl, colsev=colsev)
    dict_pos = dict(zip(list_severity_sorted, list_pos))

    return dict_pos

def get_xticklabels(df_beta_methyl, 
                    colsample="ID", 
                    colsev="Severity_visit"):
    df_beta_methyl_reset_idx = df_beta_methyl.reset_index(drop=False).rename(columns={"index": colsample})
    df_num = df_beta_methyl_reset_idx.groupby(colsev)[colsample].unique().apply(len).reset_index(drop=False)
    dict_num = dict(zip(df_num[colsev], df_num[colsample]))
    
    list_severity_sorted = get_sorted_severity(df_beta_methyl, colsev=colsev)
    list_xticklabels = list()
    dict_name_sev = {"Severe": "Sev.", "Convalescent": "Conval."}
    dict_name_vis = {"First": "Acute", "Last": "Recov."}
    for key in list_severity_sorted:
        value = dict_num.get(key, None)
        if '_' in key:
            sev, vis = key.split('_')
            vis = dict_name_vis.get(vis, vis)
            vis = f"\n{vis}"
        else:
            sev = key
            vis = ''
        
        sev = dict_name_sev.get(sev, sev)
        xticklabel = f"{sev}{vis}\n({value})"
        list_xticklabels.append(xticklabel)
    
    return list_xticklabels

def get_list_comb_sev_compare(df_beta_methyl, colsev="Severity_visit"):
    list_severity_sorted = get_sorted_severity(df_beta_methyl, colsev=colsev)
    list_comb = list(combinations(list_severity_sorted, r=2))
    list_comb_select = list(filter(lambda x: (x[0].split("_")[-1] == x[1].split("_")[-1]), list_comb))

    return list_comb_select

# %%
list_files_methylcpgmin = get_list_files_methylcpgmin(dir_methylcpgmin)
list_methyl_sample = get_list_methyl_sample(list_files_methylcpgmin)
df_beta_all_hyper = get_dataframe_methylation_beta_samples(infilenamehyper, list_methyl_sample)
df_beta_all_hypo = get_dataframe_methylation_beta_samples(infilenamehypo, list_methyl_sample)
list_metainfo = list(df_beta_all_hyper.iloc[:, -8:].columns)
dict_markers_overlap = get_dict_deg_dmp_overlap_markers(file_deg_dmp_overlap)

marker_overlap_hyper = list(set(dict_markers_overlap.keys()).intersection(set(list(df_beta_all_hyper.columns))))
marker_overlap_hypo = list(set(dict_markers_overlap.keys()).intersection(set(list(df_beta_all_hypo.columns))))

with open(path_marker_target, mode="r") as fr:
    list_target_methyl = list(map(lambda x: x.rstrip(), fr.readlines()))

dict_marker_target = dict()
for marker, markers in dict_markers_overlap.items():
    if marker in list_target_methyl:
        target = "-".join([x.split("_")[-1] for x in markers])
        dict_marker_target[marker] = target
        
df_beta_selec_target = df_beta_all_hypo[list(dict_marker_target.keys()) + list_metainfo]
df_beta_selec_target = df_beta_selec_target[~df_beta_selec_target.index.str.contains("V3")]
df_beta_selec_target = df_beta_selec_target[~df_beta_selec_target.index.str.contains("V4")]

df_beta_methyl_all = make_dataframe_stat_test(df_beta_selec_target, list_target_methyl)
save_stat_test_result(df_beta_methyl_all, outdir, "methylation_meanbeta_differences_between_severity_groups.txt", dict_marker_target)

plt.rcParams["font.size"] = 10
cm = 1/2.54
width = 55*cm
height = 65*cm
figsize = (width, height)
plot_linewidth = 1

fig = plt.figure(figsize = figsize)
col = 100
row = 100
gsfig = gridspec.GridSpec(row, col, left = 0, right = 1, bottom = 0, top = 1, wspace = 1, hspace = 1)

gs1 = gsfig[0:60, 0:100]
ax1 = fig.add_subplot(gs1)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['bottom'].set_visible(False)
ax1.axis("off")
ax1.text(-0.05, 1.10, "A", transform=ax1.transAxes,
        fontsize=plt.rcParams["font.size"]+30, fontweight='bold', va='top', ha='right')

gs2 = gsfig[70:100, 0:100]
ax2 = fig.add_subplot(gs2)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.axis("off")
ax2.text(-0.05, 1.10, "B", transform=ax2.transAxes,
        fontsize=plt.rcParams["font.size"]+30, fontweight='bold', va='top', ha='right')

colsev = "Severity_visit"
nrows = 2
ncols = int(np.ceil(len(list_target_methyl) / nrows))
gs_subplots1 = gridspec.GridSpecFromSubplotSpec(nrows, ncols, subplot_spec=gs1, wspace=0.6, hspace=0.6)
axes1 = np.array([[fig.add_subplot(gs_subplots1[i, j]) for j in range(ncols)] for i in range(nrows)])
axes1 = axes1.flatten()

dict_pos = get_dict_pos(df_beta_selec_target, list_pos=[0, 1, 2, 3])
dict_palette = get_dict_palette(df_beta_selec_target)

categories = list(dict_pos.keys())
positions = list(dict_pos.values())
colors = [dict_palette[cat] for cat in categories]
list_sev_vis = list(df_beta_selec_target[colsev].unique())
order = [list_sev_vis[-1]] + list_sev_vis[:-1]

for i, (marker, ax) in enumerate(zip(list_target_methyl, axes1)):
    plot_data = {cat: [] for cat in categories}
    for category in categories:
        values = df_beta_selec_target[df_beta_selec_target[colsev] == category][marker].dropna().values
        plot_data[category] = values

    sns.boxplot(data=df_beta_selec_target, x=colsev, y=marker, order=order, palette=colors, ax=ax, fliersize=0, width=0.5, zorder=100)

    for category, color in zip(categories, colors):
        x = np.full_like(plot_data[category], dict_pos[category], dtype=float)
        ax.scatter(x, plot_data[category], color=color, edgecolor='black', label=category, alpha=0.5, s=100, zorder=999)
    
    dict_beta_sev_grp = df_beta_selec_target.groupby(colsev)[marker].apply(np.array).to_dict()
    box_pairs = df_beta_methyl_all.loc[marker, "comp"]
    list_idx_selec = [idx for idx, x in zip(range(len(box_pairs)), box_pairs) if x[1]!="Healthy"]
    box_pairs_selec = [box_pairs[idx] for idx in list_idx_selec]
    pvalues = df_beta_methyl_all.loc[marker, "pval"]
    pvalues_selec = [pvalues[idx] for idx in list_idx_selec]
    
    add_stat_annotation(ax,
                        data=df_beta_selec_target, 
                        x=colsev,
                        y=marker,
                        order=order,
                        box_pairs=box_pairs_selec,
                        perform_stat_test=False, 
                        pvalues=pvalues_selec,
                        text_offset=-0.1,
                        line_offset=0.05,
                        fontsize=plt.rcParams["font.size"]+10,
                        loc='outside', 
                        verbose=0)

    ax.set_xlabel("")    
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xticks(positions)
    list_xticklabels = get_xticklabels(df_beta_selec_target, colsample="Sample_ID")
    ax.set_xticklabels(list_xticklabels, fontsize=plt.rcParams["font.size"]+8)
    ax.tick_params(axis='y', labelsize=plt.rcParams["font.size"]+10)
    ax.set_ylabel("Proportion of\nMethylated CpGs (%)", fontsize=plt.rcParams["font.size"]+12)
    markername = f"{marker}\n$\it({dict_marker_target.get(marker, marker)})$\n\n"
    ax.set_title(markername, fontsize=plt.rcParams["font.size"]+15)
    ax.set_ylim(-10, 119)
    ax.grid(axis="y")
    ax.set_axisbelow(True)

legend_elements = [
    mpatches.Patch(color='darkgrey', label='Healthy controls'),
    mpatches.Patch(color='forestgreen', label='Mild (Acute)'),
    mpatches.Patch(color='firebrick', label='Severe (Acute)'),
    mpatches.Patch(color='royalblue', label='Convalescent')
]

axes1[-1].legend(handles=legend_elements, 
           loc='center',
           bbox_to_anchor=(0.4, 0.4), 
           title="Groups", 
           fontsize=plt.rcParams["font.size"]+15, 
           title_fontsize=plt.rcParams["font.size"]+17,
           frameon=False)

axes1[-1].axis("off")

from rna_severity import (get_dict_palette, get_dict_pos, get_xticklabels,
                          main, make_dataframe_stat_test_rna,
                          save_stat_test_result_rna)

col_id = "ID"
path_drop_samples = "/BiO/Access/kyungwhan1998/Infectomics/Resources/Scripts/Final/list_remove_samples.txt"
path_exp = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/DEGExtract/DEGExtract.RNA_samples_with_Methyl.DEG_by_Sex.Control_M.Case_F.Cov_Age.mincount_1.20240402.tsv.normcount"
path_methyl_marker = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.DMP_DEG_LOO_annot.20240402.tsv"
path_meta = "/BiO/Access/kyungwhan1998/Infectomics/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"
outdir = f"/BiO/Access/kyungwhan1998/Infectomics/Results/InfectomicsPaper1/20240906"
os.makedirs(outdir, exist_ok=True)
list_targeted_marker = ['ENSG00000115604.12_IL18R1', 'ENSG00000115607.10_IL18RAP',  'ENSG00000176928.7_GCNT4', 'ENSG00000118689.15_FOXO3',  'ENSG00000140564.13_FURIN', 'ENSG00000271303.2_SRXN1']
df_rna_count_meta = main(path_drop_samples, path_exp, path_methyl_marker, path_meta)
df_rna_count_meta = df_rna_count_meta[~df_rna_count_meta[col_id].str.contains("V3")]
df_rna_count_meta = df_rna_count_meta[~df_rna_count_meta[col_id].str.contains("V4")]
df_gene_diffexp_sorted = make_dataframe_stat_test_rna(df_rna_count_meta, list_targeted_marker)
save_stat_test_result_rna(df_gene_diffexp_sorted, outdir, "rna_meanexp_differences_between_severity_groups.txt")

nrows = 1
ncols = int(np.ceil(len(list_targeted_marker) / nrows))
gs_subplots2 = gridspec.GridSpecFromSubplotSpec(nrows, ncols, subplot_spec=gs2, wspace=0.7)
axes2 = np.array([[fig.add_subplot(gs_subplots2[i, j]) for j in range(ncols)] for i in range(nrows)])
axes2 = axes2.flatten()

dict_pos_rna = get_dict_pos(df_rna_count_meta, colsev=colsev, list_pos=[0, 1, 2])
dict_palette_rna = get_dict_palette(df_rna_count_meta, colsev=colsev)

categories = list(dict_pos_rna.keys())
positions = list(dict_pos_rna.values())
colors = [dict_palette_rna[cat] for cat in categories]
order = list(df_rna_count_meta[colsev].unique())

for i, (marker, ax2) in enumerate(zip(list_targeted_marker, axes2)):
    df_rna_count_meta_copy = df_rna_count_meta.copy()
    df_gene_diffexp_sorted_copy = df_gene_diffexp_sorted.copy()
    
    plot_data = {cat: [] for cat in categories}
    for category in categories:
        values = df_rna_count_meta_copy[df_rna_count_meta_copy[colsev] == category][marker].dropna().values
        plot_data[category] = values

    sns.boxplot(data=df_rna_count_meta_copy, x=colsev, y=marker, order=order, palette=colors, ax=ax2, fliersize=0, width=0.5, zorder=100)

    # Plot scatter plots
    for category, color in zip(categories, colors):
        x = np.full_like(plot_data[category], dict_pos_rna[category], dtype=float)
        ax2.scatter(x, plot_data[category], color=color, edgecolor='black', label=category, alpha=0.5, s=200, zorder=999)

    box_pairs = df_gene_diffexp_sorted_copy.loc[marker, "comp"]
    pvalues = df_gene_diffexp_sorted_copy.loc[marker, "padj"]
    
    add_stat_annotation(ax2, 
                        data=df_rna_count_meta_copy, 
                        x=colsev,
                        y=marker,
                        order=order,
                        box_pairs=box_pairs,
                        perform_stat_test=False,
                        pvalues=pvalues,
                        text_offset=-0.1,
                        line_offset=0.05,
                        text_format='star',
                        fontsize=plt.rcParams["font.size"]+10,
                        loc='inside', 
                        verbose=0)
    
    ax.set_xlabel("")      
    list_xticklabels = get_xticklabels(df_rna_count_meta_copy, colsample="ID")
    ax2.spines[['right', 'top']].set_visible(False)
    ax2.set_xticks(positions)
    ax2.set_xticklabels(list_xticklabels, fontsize=plt.rcParams["font.size"]+8)
    ax2.tick_params(axis='y', labelsize=plt.rcParams["font.size"]+10)
    ax2.set_xlabel("Severity", fontsize=0)
    ax2.set_ylabel("Expression", fontsize=plt.rcParams["font.size"]+15)
    ax2.set_title(f"$\it{marker.split('_')[1]}$\n", fontsize=plt.rcParams["font.size"]+15)
    ax2.grid(axis="y")
    ax2.set_axisbelow(True)

plt.savefig("/BiO/Access/kyungwhan1998/Infectomics/Results/InfectomicsPaper1/20240906/Figure2.png", bbox_inches="tight", dpi=600)
plt.show()
plt.close()

# %%

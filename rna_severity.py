# %%
import os
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import ranksums, shapiro, ttest_ind
from statsmodels.stats.multitest import fdrcorrection
from statannot import add_stat_annotation
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns

def get_drop_samples(path_drop_samples):
    with open(path_drop_samples, mode="r") as frmv:
        list_drop_samples = [x.rstrip() for x in frmv.readlines()]
    
    return list_drop_samples

def read_rna_count_matrix(path_exp, list_drop_samples):
    df_rna_count = pd.read_csv(path_exp, sep="\t")
    if list(df_rna_count.index)[0] == 0:
        df_rna_count = pd.read_csv(path_exp, sep="\t", index_col=[0])
    list_columns = list(df_rna_count.columns)
    list_drop_intsc = list(set(list_columns).intersection(set(list_drop_samples)))
    df_rna_count  = df_rna_count.drop(columns=list_drop_intsc)

    return df_rna_count

def select_rna_count_matrix(df_rna_count, select_pattern, delim_id="-", namepos=0):
    list_colname = list(df_rna_count.columns)
    list_colname_filt = list(filter(lambda x: x.split(delim_id)[namepos][0] != select_pattern if len(x.split("-")) > 1 else x, list_colname))
    df_rna_count = df_rna_count.loc[:, list_colname_filt]

    return df_rna_count

def get_list_gene_markers(file_deg_dmp_overlap):
    df_gene_marker = pd.read_csv(file_deg_dmp_overlap, sep="\t")
    list_gene_marker = df_gene_marker["RNA"].tolist()
    
    return list_gene_marker

def filter_rna_count_matrix(df_rna_count, list_gene):
    df_rna_count_filt = df_rna_count[df_rna_count.index.isin(list_gene)]
    
    return df_rna_count_filt

def transpose_rna_count_matrix(df_rna_count):
    df_rna_count_transposed = df_rna_count.T
    
    return df_rna_count_transposed

def reset_idx_rna_count_matrix(df_rna_count, colsample="Sample_ID"):
    df_rna_count_reset_idx = df_rna_count.reset_index(drop=False).rename(columns={"index": colsample})
    
    return df_rna_count_reset_idx
    
def merge_rna_count_matrix_meta(df_rna_count, path_meta, colsample="Sample_ID", colmeta="Project_ID_Alias", colsev="Severity_Binary", list_drop_samples=["C19-C014-V1", "C19-C059-V2"]):
    df_meta = pd.read_csv(path_meta, sep="\t")
    df_meta_rename_colsample = df_meta.rename(columns={colmeta: colsample})
    df_meta_rename_colsample_set_idx = df_meta_rename_colsample.set_index(colsample)
    list_drop_samples_select = list()
    for sample in set(df_meta_rename_colsample_set_idx.index):
        if sample in list_drop_samples:
            list_drop_samples_select.append(sample)
    df_meta_rename_colsample_drop_no_exp = df_meta_rename_colsample_set_idx.drop(index=list_drop_samples_select)
    df_meta_rename_colsample_drop_no_exp = df_meta_rename_colsample_drop_no_exp.reset_index(drop=False).rename(columns={"index": colsample})
    df_meta_rename_colsample_drop_no_exp_drop_by_visit = df_meta_rename_colsample_drop_no_exp.loc[~df_meta_rename_colsample_drop_no_exp["Visit_order"].isin(["-", "LongCOVID"])]
    df_meta_rename_colsample_drop_no_exp.loc[(df_meta_rename_colsample_drop_no_exp["Visit_order"] == "Recover"), colsev] = "Convalescent"
    df_meta_rename_colsample_drop_no_exp.loc[(df_meta_rename_colsample_drop_no_exp["HealthyControl"] == 1) | (df_meta_rename_colsample_drop_no_exp["Vaccination"] == 1), colsev] = "Healthy"
    
    df_merged = pd.merge(df_rna_count, df_meta_rename_colsample_drop_no_exp, how="inner", on=colsample)
    
    return df_merged

def find_meandiff(a, b):
    deltamean = np.mean(b) - np.mean(a)

    return deltamean

def make_dataframe_stat_test_rna(df_rna_count_meta, list_targeted_gene, colsev="Severity_Binary"): 
    dict_marker_diffexp_pval = dict()
    for marker in list_targeted_gene:
        df_group_visit = df_rna_count_meta.groupby(colsev)[marker].apply(np.array).reset_index(drop=False)
        list_comb_group = list(combinations(df_group_visit[colsev].values, r=2))
        list_comb_group = list(map(lambda x: x[::-1], list_comb_group))
        list_comb_marker_val = list(combinations(df_group_visit[marker].values, r=2))
        list_comb_marker_val = list(map(lambda x: x[::-1], list_comb_marker_val))
        list_meandiff = list(map(lambda x: find_meandiff(*x), list_comb_marker_val))
        list_stat = list(map(lambda x: ranksums(*x)[0], list_comb_marker_val))
        list_pval = list(map(lambda x: ranksums(*x)[1], list_comb_marker_val))
        _, list_padj = fdrcorrection(list_pval)
        list_meandiff = list(map(lambda x: round(x, 3), list_meandiff))
        list_stat = list(map(lambda x: round(x, 3), list_stat))
        dict_stat_res = {"comp": list_comb_group, "diff": list_meandiff, "stat": list_stat, "pval": list_pval, "padj": list_padj}
        dict_marker_diffexp_pval[marker] = dict_stat_res
        
    df_marker_diffexp = pd.DataFrame.from_dict(dict_marker_diffexp_pval).T
    
    return df_marker_diffexp

def save_stat_test_result_rna(df_marker_stat, outdir, filename, list_genename=[]):
    df_marker_stat_melted = pd.DataFrame(columns = ["Gene", "Comp1", "Comp2", "diff", "stat", "pval", "padj"])

    ind_df = 0
    for genename, row in df_marker_stat.iterrows():
        comp_list = row["comp"]
        diff_list = row["diff"]
        stat_list = row["stat"]
        pval_list = row["pval"]
        padj_list = row["padj"]
        pval_list = list(map(lambda x: round(x, 3) if x >= 0.001 else ('{:.3e}'.format(x) if int(x) == 0 else x), pval_list))
        padj_list = list(map(lambda x: round(x, 3) if x >= 0.001 else ('{:.3e}'.format(x) if int(x) == 0 else x), padj_list))
        for comp, diff, stat, pval, padj in zip(comp_list, diff_list, stat_list, pval_list, padj_list):
            comp1, comp2 = comp
            values_row = [genename, comp1, comp2, diff, stat, pval, padj]
            df_marker_stat_melted.loc[ind_df, :] = values_row
            ind_df += 1
    
    outfilepath = os.path.join(outdir, filename)
    first_col = list(df_marker_stat_melted.columns)[0]
    df_marker_stat_melted = df_marker_stat_melted.sort_values(by=[first_col], ascending=False)
    if list_genename != []:
        df_marker_stat_melted_set_index = df_marker_stat_melted.set_index("Gene")
        df_marker_stat_melted_selected = df_marker_stat_melted_set_index.loc[list_genename]
        df_marker_stat_melted = df_marker_stat_melted_selected.reset_index(drop=False)
    df_marker_stat_melted.to_csv(outfilepath, sep="\t", index=False)
    
    return df_marker_stat_melted

def make_dataframe_normality_test_rna(df_rna_count_meta, list_target_marker, colsev="Severity_Binary"):
    dict_marker_norm = dict()
    for marker in list_target_marker:
        df_group_visit = df_rna_count_meta.groupby(colsev)[marker].apply(np.array).reset_index(drop=False)
        list_marker_val = list(df_group_visit[marker].values)
        list_norm_pval = list(map(lambda x: shapiro(x)[1], list_marker_val))
        _, list_padj = fdrcorrection(list_norm_pval)
        dict_marker_norm[marker] = list_padj
        list_groups = list(df_group_visit[colsev])
    df_marker_norm = pd.DataFrame.from_dict(dict_marker_norm).T
    df_marker_norm= df_marker_norm.applymap(lambda x: round(x, 3) if round(x, 3) != 0.0 else '{:.3e}'.format(x))
    df_marker_norm.columns = list_groups
    
    return df_marker_norm

def save_norm_test_result_rna(df_marker_norm, outdir, filename):
    outfilepath = os.path.join(outdir, filename)
    df_marker_norm.to_csv(outfilepath, sep="\t", index=True)
    
    return df_marker_norm

def get_xticklabels(df_rna_count_meta, dict_name_vis = {"First": " Acute", "Last": " Recovery"}, colsample="Sample_ID", colsev="Severity_Binary"):
    df_num = df_rna_count_meta.groupby(colsev)[colsample].unique().apply(len).reset_index(drop=False)
    dict_num = dict(zip(df_num[colsev], df_num[colsample]))

    list_sev_vis = list(df_rna_count_meta[colsev].unique())
    list_sev_vis = sorted(
        list_sev_vis,
        key=lambda x: (
            0 if x == "Healthy" else 2 if x == "Convalescent" else 1,  # Healthy/Convalescent to top/bottom
            x.split("_")[1] if "_" in x else "",  # Visit type
            0 if x.startswith("Mild") else 1  # Mild before Severe
        )
    )
    
    list_xticklabels = list()
    dict_name_sev = {"Mild": "Mild-Moderate", "Severe": "Severe-Critical", "Convalescent": "Convalescent"}
    
    for key in list_sev_vis:
        value = dict_num.get(key, None)
        if '_' in key:
            sev, vis = key.split('_')
            vis = dict_name_vis.get(vis, vis)
        else:
            sev = key
            vis = ''
        
        sev = dict_name_sev.get(sev, sev)
        xticklabel = f"{sev}{vis}"
        list_xticklabels.append(xticklabel)
    
    return list_xticklabels

def get_dict_palette(df_rna_count_meta, colsev="Severity_Binary"):
    list_sev_vis = df_rna_count_meta[colsev].unique()
    list_sev_vis_sort = sorted(list_sev_vis, key=lambda x:x.split("_")[-1])
    list_sev_vis_shift = list_sev_vis_sort[1:] + [list_sev_vis_sort[0]]
    
    list_colors = list()
    for sev_vis in list_sev_vis_shift:
        if sev_vis.split("_")[0]=="Healthy":
            list_colors.append("darkgrey")
        elif sev_vis.split("_")[0]=="Mild":
            list_colors.append("forestgreen")
        elif sev_vis.split("_")[0]=="Severe":
            list_colors.append("firebrick") 
        else: 
            list_colors.append("royalblue")
            
    dict_palette = dict(zip(list_sev_vis_shift, list_colors))

    return dict_palette

def get_dict_pos(df_rna_count_meta, colsev="Severity_Binary", list_pos=[0, 0.5, 1.5, 2.0, 3.0]):
    list_sev_vis = df_rna_count_meta[colsev].unique()
    list_sev_vis_sort = sorted(list_sev_vis, key=lambda x:x.split("_")[-1])
    list_sev_vis_shift = list_sev_vis_sort[1:] + [list_sev_vis_sort[0]]
    list_num = list(range(len(list_sev_vis_shift)))
    dict_pos = dict(zip(list_sev_vis_shift, list_pos))
    
    return dict_pos

def get_dictionary_gene_pvalsig(df_gene_diffexp_first, target_pair = ("Mild_First", "Severe_First")):
    dict_gene_pvalsig = dict()
    list_all_pairs = list(map(lambda x: sorted(x), df_gene_diffexp_first["comp"][0]))
    idx_target = list_all_pairs.index(sorted(target_pair))
    genes = df_gene_diffexp_first["padj"].index
    for gene in genes:
        fdr = df_gene_diffexp_first["padj"].loc[gene][idx_target]
        if fdr < 0.05:
            dict_gene_pvalsig[gene] = True
        else:
            dict_gene_pvalsig[gene] = False
    
    return dict_gene_pvalsig

def main(path_drop_samples, path_exp, file_deg_dmp_overlap, path_meta):
    list_drop_samples = get_drop_samples(path_drop_samples)
    df_rna_count = read_rna_count_matrix(path_exp, list_drop_samples)
    list_gene_markers = get_list_gene_markers(file_deg_dmp_overlap)
    df_rna_count_gene_filtered = filter_rna_count_matrix(df_rna_count, list_gene_markers)
    df_rna_count_transposed = transpose_rna_count_matrix(df_rna_count_gene_filtered)
    df_rna_count_transposed_reset_idx = reset_idx_rna_count_matrix(df_rna_count_transposed)
    df_rna_count_meta = merge_rna_count_matrix_meta(df_rna_count_transposed_reset_idx, path_meta)
    
    return df_rna_count_meta

def plot_rna_difference(df_rna_count_meta, df_gene_diffexp, list_targeted_marker, gs, fig, colsample="Sample_ID", colsev="Severity_Binary", pval="pval", nrows=1, list_pos=[0, 1, 2, 3], order=['Healthy', 'Mild', 'Severe', 'Convalescent']):
    ncols = int(np.ceil(len(list_targeted_marker) / nrows))
    gs_subplots = gridspec.GridSpecFromSubplotSpec(nrows, ncols, subplot_spec=gs, wspace=0.6)
    axes = np.array([[fig.add_subplot(gs_subplots[i, j]) for j in range(ncols)] for i in range(nrows)])
    axes = axes.flatten()

    dict_pos_rna = get_dict_pos(df_rna_count_meta, colsev=colsev, list_pos=list_pos)
    dict_palette_rna = get_dict_palette(df_rna_count_meta, colsev=colsev)

    categories = list(dict_pos_rna.keys())
    positions = list(dict_pos_rna.values())
    colors = [dict_palette_rna[cat] for cat in categories]
    
    for i, (marker, ax) in enumerate(zip(list_targeted_marker, axes)):
        df_rna_count_meta_copy = df_rna_count_meta.copy()
        df_gene_diffexp_copy = df_gene_diffexp.copy()
        
        plot_data = {cat: [] for cat in categories}
        for category in categories:
            values = df_rna_count_meta_copy[df_rna_count_meta_copy[colsev] == category][marker].dropna().values
            plot_data[category] = values

        sns.boxplot(data=df_rna_count_meta_copy, x=colsev, y=marker, order=order, palette=colors, ax=ax, fliersize=0, width=0.5, zorder=100)

        # Plot scatter plots
        for category, color in zip(categories, colors):
            x = np.full_like(plot_data[category], dict_pos_rna[category], dtype=float)
            ax.scatter(x, plot_data[category], color=color, edgecolor='black', label=category, alpha=0.5, s=50, zorder=999)

        box_pairs = df_gene_diffexp_copy.loc[marker, "comp"]
        list_idx_selec = [idx for idx, x in zip(range(len(box_pairs)), box_pairs) if x[1]!="Healthy"]
        box_pairs_selec = [box_pairs[idx] for idx in list_idx_selec]
        pvalues = df_gene_diffexp_copy.loc[marker, pval]
        pvalues_selec = [pvalues[idx] for idx in list_idx_selec]
        
        add_stat_annotation(ax, 
                            data=df_rna_count_meta_copy, 
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

        ax.set_xlabel("", fontsize=0)      
        list_xticklabels = get_xticklabels(df_rna_count_meta_copy, colsample=colsample)
        ax.spines[['right', 'top']].set_visible(False)
        ax.set_xticks(positions)
        ax.set_xticklabels(list_xticklabels, fontsize=plt.rcParams["font.size"]+2, rotation = 45, ha="right", rotation_mode="anchor")
        ax.tick_params(axis='y', labelsize=plt.rcParams["font.size"]+2)
        ax.set_xlabel("", fontsize=0)
        
        if i == 0:
            ax.set_ylabel("VST count", fontsize=plt.rcParams["font.size"]+3)
        else:
            ax.set_ylabel("", fontsize=0)
        ax.set_title(f"$\it{marker.split('_')[1]}$\n", fontsize=plt.rcParams["font.size"]+5)
        ax.grid(axis="y")
        ax.set_axisbelow(True)
    
    return axes
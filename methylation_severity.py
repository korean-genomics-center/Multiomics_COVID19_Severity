# %%
import glob
import gzip
import os
import pickle
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import shapiro, ttest_ind
from statsmodels.stats.multitest import fdrcorrection

from statannot import add_stat_annotation
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import seaborn as sns

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

def get_dataframe_methylation_beta_samples(path_sev_info, infilename, list_selected_samples):
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

def make_dataframe_normality_test(df_marker, list_target_marker, colsev="Severity_visit"):
    dict_marker_norm = dict()
    for marker in list_target_marker:
        df_group_visit = df_marker.groupby(colsev)[marker].apply(np.array).reset_index(drop=False)
        list_marker_val = list(df_group_visit[marker].values)
        list_norm_pval = list(map(lambda x: shapiro(x)[1], list_marker_val))
        _, list_padj = fdrcorrection(list_norm_pval)
        dict_marker_norm[marker] = list_padj
        
    list_groups = list(df_group_visit[colsev])
    df_marker_norm = pd.DataFrame.from_dict(dict_marker_norm).T
    df_marker_norm.columns = list_groups
    df_marker_norm= df_marker_norm.applymap(lambda x: round(x, 3) if round(x, 3) != 0.0 else '{:.3e}'.format(x))
    
    return df_marker_norm

def save_norm_test_result(df_marker_norm, outdir, filename):
    outfilepath = os.path.join(outdir, filename)
    df_marker_norm.to_csv(outfilepath, sep="\t", index=True)
    
    return df_marker_norm

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
        _, list_padj = fdrcorrection(list_pval)
        list_meandiff = list(map(lambda x: round(x, 3), list_meandiff))
        list_stat = list(map(lambda x: round(x, 3), list_stat))
        dict_stat_res = {"comp": list_comb_group, "diff": list_meandiff, "stat": list_stat, "pval": list_pval, "padj": list_padj}
        dict_marker_stat[marker] = dict_stat_res
                            
    df_marker_stat = pd.DataFrame.from_dict(dict_marker_stat).T
    
    return df_marker_stat

def save_stat_test_result(df_marker_stat, outdir, filename, dict_marker_target, list_methyl=[]):
    df_marker_stat_melted = pd.DataFrame(columns = ["CpG", "Gene", "Comp1", "Comp2", "diff", "stat", "pval", "padj"])

    ind_df = 0
    for cpgname, row in df_marker_stat.iterrows():
        comp_list = row["comp"]
        diff_list = row["diff"]
        stat_list = row["stat"]
        pval_list = row["pval"]
        padj_list = row["padj"]
        pval_list = list(map(lambda x: round(x, 3) if x >= 0.001 else ('{:.3e}'.format(x) if int(x) == 0 else x), pval_list))
        padj_list = list(map(lambda x: round(x, 3) if x >= 0.001 else ('{:.3e}'.format(x) if int(x) == 0 else x), padj_list))
        for comp, diff, stat, pval, padj in zip(comp_list, diff_list, stat_list, pval_list, padj_list):
            comp1, comp2 = comp
            genename = dict_marker_target[cpgname]
            values_row = [cpgname, genename, comp1, comp2, diff, stat, pval, padj]
            df_marker_stat_melted.loc[ind_df, :] = values_row
            ind_df += 1
    
    outfilepath = os.path.join(outdir, filename)
    first_col = list(df_marker_stat_melted.columns)[0]
    df_marker_stat_melted = df_marker_stat_melted.sort_values(by=[first_col], ascending=False)
    if list_methyl != []:
        df_marker_stat_melted_set_index = df_marker_stat_melted.set_index("CpG")
        df_marker_stat_melted_selected = df_marker_stat_melted_set_index.loc[list_methyl]
        df_marker_stat_melted = df_marker_stat_melted_selected.reset_index(drop=False)
    df_marker_stat_melted.to_csv(outfilepath, sep="\t", index=False)
    
    return df_marker_stat_melted


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

def get_dict_pos(df_beta_methyl, list_pos=[0, 1, 1.5, 2.5, 3, 4], colsev="Severity_visit"):
    list_severity_sorted = get_sorted_severity(df_beta_methyl, colsev=colsev)
    dict_pos = dict(zip(list_severity_sorted, list_pos))

    return dict_pos

def get_xticklabels(df_beta_methyl, dict_name_vis={"First": "\nAcute", "Last": "\nRecovery"}, colsample="Sample_ID", colsev="Severity_visit"):
    if colsample not in list(df_beta_methyl.columns):
        df_beta_methyl_reset_idx = df_beta_methyl.reset_index(drop=False).rename(columns={"index": colsample})
        df_beta_methyl = df_beta_methyl_reset_idx
    df_num = df_beta_methyl_reset_idx.groupby(colsev)[colsample].unique().apply(len).reset_index(drop=False)
    dict_num = dict(zip(df_num[colsev], df_num[colsample]))
    
    list_severity_sorted = get_sorted_severity(df_beta_methyl, colsev=colsev)
    list_xticklabels = list()
    dict_name_sev = {"Mild": "Mild-Moderate", "Severe": "Severe-Critical", "Convalescent": "Convalescent"}
    dict_name_vis = dict_name_vis
    for key in list_severity_sorted:
        value = dict_num.get(key, None)
        if '_' in key:
            sev, vis = key.split('_')
            vis = dict_name_vis.get(vis, vis)
            # vis = f"\n{vis}"
        else:
            sev = key
            vis = ''
        
        sev = dict_name_sev.get(sev, sev)
        # xticklabel = f"{sev}{vis}\n({value})"
        xticklabel = f"{sev}{vis}"
        list_xticklabels.append(xticklabel)
    
    return list_xticklabels

def get_list_comb_sev_compare(df_beta_methyl, colsev="Severity_visit"):
    list_severity_sorted = get_sorted_severity(df_beta_methyl, colsev=colsev)
    list_comb = list(combinations(list_severity_sorted, r=2))
    list_comb_select = list(filter(lambda x: (x[0].split("_")[-1] == x[1].split("_")[-1]), list_comb))

    return list_comb_select

def main(path_sev_info, dir_methylcpgmin, infilenamehyper, infilenamehypo):
    list_files_methylcpgmin = get_list_files_methylcpgmin(dir_methylcpgmin)
    list_methyl_sample = get_list_methyl_sample(list_files_methylcpgmin)
    df_beta_all_hyper = get_dataframe_methylation_beta_samples(path_sev_info, infilenamehyper, list_methyl_sample)
    df_beta_all_hypo = get_dataframe_methylation_beta_samples(path_sev_info, infilenamehypo, list_methyl_sample)

    return df_beta_all_hyper, df_beta_all_hypo

def plot_methyl_difference(df_beta_selec_target, df_beta_methyl_all, list_plot_marker_target, dict_marker_target, gs_n, fig, colsev = "Severity_visit", pval="pval", nrows = 2, list_pos=[0, 1, 2, 3], dict_name_vis={"First": "", "Last": ""}):
    ncols = int(np.ceil(len(list_plot_marker_target) / nrows))
    gs_subplots = gridspec.GridSpecFromSubplotSpec(nrows, ncols, subplot_spec=gs_n, wspace=0.6)
    axes = np.array([[fig.add_subplot(gs_subplots[i, j]) for j in range(ncols)] for i in range(nrows)])
    axes = axes.flatten()

    dict_pos = get_dict_pos(df_beta_selec_target, list_pos=list_pos)
    dict_palette = get_dict_palette(df_beta_selec_target)

    categories = list(dict_pos.keys())
    positions = list(dict_pos.values())
    colors = [dict_palette[cat] for cat in categories]
    order = ['Healthy', 'Mild_First', 'Severe_First', 'Convalescent']

    for i, (marker, ax) in enumerate(zip(list_plot_marker_target, axes)):
        plot_data = {cat: [] for cat in categories}
        for category in categories:
            values = df_beta_selec_target[df_beta_selec_target[colsev] == category][marker].dropna().values
            plot_data[category] = values

        sns.boxplot(data=df_beta_selec_target, x=colsev, y=marker, order=order, palette=colors, ax=ax, fliersize=0, width=0.5, zorder=100)

        for category, color in zip(categories, colors):
            x = np.full_like(plot_data[category], dict_pos[category], dtype=float)
            ax.scatter(x, plot_data[category], color=color, edgecolor='black', label=category, alpha=0.5, s=50, zorder=999)
        
        # dict_beta_sev_grp = df_beta_selec_target.groupby(colsev)[marker].apply(np.array).to_dict()
        box_pairs = df_beta_methyl_all.loc[marker, "comp"]
        list_idx_selec = [idx for idx, x in zip(range(len(box_pairs)), box_pairs) if x[1]!="Healthy"]
        box_pairs_selec = [box_pairs[idx] for idx in list_idx_selec]
        pvalues = df_beta_methyl_all.loc[marker, pval]
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
                            fontsize=plt.rcParams["font.size"]+2,
                            loc='outside', 
                            verbose=0)

        ax.set_xlabel("", fontsize=0)    
        ax.spines[['right', 'top']].set_visible(False)
        ax.set_xticks(positions)
        list_xticklabels = get_xticklabels(df_beta_selec_target, dict_name_vis=dict_name_vis, colsample="Sample_ID")
        ax.set_xticklabels(list_xticklabels, fontsize=plt.rcParams["font.size"]+2, rotation = 45, ha="right", rotation_mode="anchor")
        ax.tick_params(axis='y', labelsize=plt.rcParams["font.size"]+2)
        if i == 0:
             ax.set_ylabel("Proportion of\nMethylated CpGs (%)", fontsize=plt.rcParams["font.size"]+3)
        else:
             ax.set_ylabel("", fontsize=0)
       
        markername = f"{marker}\n$\it({dict_marker_target.get(marker, marker)})$\n\n"
        ax.set_title(markername, fontsize=plt.rcParams["font.size"]+5)
        ax.set_ylim(-10, 119)
        ax.grid(axis="y")
        ax.set_axisbelow(True)

    return axes
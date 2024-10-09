import os
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import ranksums
from statsmodels.stats.multitest import fdrcorrection


def get_drop_samples(path_drop_samples):
    with open(path_drop_samples, mode="r") as frmv:
        list_drop_samples = [x.rstrip() for x in frmv.readlines()]
    
    return list_drop_samples

def read_rna_count_matrix(path_exp, list_drop_samples):
    df_rna_count = pd.read_csv(path_exp, sep="\t")
    list_columns = list(df_rna_count.columns)
    list_drop_intsc = list(set(list_columns).intersection(set(list_drop_samples)))
    df_rna_count  = df_rna_count.drop(columns=list_drop_intsc)

    return df_rna_count

def select_rna_count_matrix(df_rna_count, select_pattern, delim_id="-", namepos=0):
    list_colname = list(df_rna_count.columns)
    list_colname_filt = list(filter(lambda x: x.split(delim_id)[namepos][0] != select_pattern if len(x.split("-")) > 1 else x, list_colname))
    df_rna_count = df_rna_count.loc[:, list_colname_filt]

    return df_rna_count

def get_list_gene_markers(path_methyl_marker):
    df_gene_marker = pd.read_csv(path_methyl_marker, sep="\t")
    list_gene_marker = df_gene_marker["RNA"].tolist()
    
    return list_gene_marker

def filter_rna_count_matrix(df_rna_count, list_gene):
    df_rna_count_filt = df_rna_count[df_rna_count.index.isin(list_gene)]
    
    return df_rna_count_filt

def transpose_rna_count_matrix(df_rna_count):
    df_rna_count_transposed = df_rna_count.T
    
    return df_rna_count_transposed

def reset_idx_rna_count_matrix(df_rna_count, colsample="ID"):
    df_rna_count_reset_idx = df_rna_count.reset_index(drop=False).rename(columns={"index": colsample})
    
    return df_rna_count_reset_idx
    
def merge_rna_count_matrix_meta(df_rna_count, path_meta, colsample="ID", list_drop_samples=["C19-C014-V1", "C19-C059-V2"]):
    df_meta = pd.read_csv(path_meta, sep="\t")
    df_meta_rename_colsample = df_meta.rename(columns={"Sample_ID": colsample})
    df_meta_rename_colsample_set_idx = df_meta_rename_colsample.set_index(colsample)
    df_meta_rename_colsample_drop_no_exp = df_meta_rename_colsample_set_idx.drop(index=list_drop_samples)
    df_meta_rename_colsample_drop_no_exp = df_meta_rename_colsample_drop_no_exp.reset_index(drop=False).rename(columns={"index": colsample})
    df_meta_rename_colsample_drop_no_exp["Severity_visit"] = df_meta_rename_colsample_drop_no_exp["Severity_visit"].apply(lambda x: "Convalescent" if "Convalescent" in x else x)
    df_merged = pd.merge(df_rna_count, 
                         df_meta_rename_colsample_drop_no_exp, 
                         how="inner", 
                         on=colsample)
    
    return df_merged

def find_meandiff(a, b):
    deltamean = np.mean(b) - np.mean(a)

    return deltamean

def make_dataframe_stat_test_rna(df_rna_count_meta, list_targeted_gene, colsev="Severity_visit"):
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
        list_sig, list_padj = fdrcorrection(list_pval)
        dict_stat_res = {"comp": list_comb_group, "diff": list_meandiff, "stat": list_stat, "pval": list_pval, "padj": list_padj, "is_sig": list_sig}
        dict_marker_diffexp_pval[marker] = dict_stat_res
                            
    df_marker_diffexp = pd.DataFrame.from_dict(dict_marker_diffexp_pval).T
    df_marker_diffexp_sorted = df_marker_diffexp.sort_values(by=["pval"], ascending=True)
    
    return df_marker_diffexp_sorted

def save_stat_test_result_rna(df_marker_stat_sorted, outdir, filename):
    df_marker_stat_sorted_melted = pd.DataFrame(columns = ["Gene", "Comp1", "Comp2", "diff", "stat", "pval", "padj", "is_sig"])

    ind_df = 0
    for genename, row in df_marker_stat_sorted.iterrows():
        comp_list = row["comp"]
        diff_list = row["diff"]
        stat_list = row["stat"]
        pval_list = row["pval"]
        padj_list = row["padj"]
        is_sig_list = row["is_sig"]
        
        for comp, diff, stat, pval, padj, is_sig in zip(comp_list, diff_list, stat_list, pval_list, padj_list, is_sig_list):
            comp1, comp2 = comp
            values_row = [genename, comp1, comp2, diff, stat, pval, padj, is_sig]
            df_marker_stat_sorted_melted.loc[ind_df, :] = values_row
            ind_df += 1
    
    outfilepath = os.path.join(outdir, filename)
    first_col = list(df_marker_stat_sorted_melted.columns)[0]
    df_marker_stat_sorted_melted = df_marker_stat_sorted_melted.sort_values(by=[first_col], ascending=False)
    df_marker_stat_sorted_melted.to_csv(outfilepath, sep="\t", index=False)
    
    return df_marker_stat_sorted_melted

    
    outfilepath = os.path.join(outdir, filename)
    df_marker_diffexp_sorted_melted.to_csv(outfilepath, sep="\t", index=False)
    
    return df_marker_diffexp_sorted_melted

def get_xticklabels(df_rna_count_meta, colsample="ID", colsev="Severity_visit"):
    df_num = df_rna_count_meta.groupby(colsev)[colsample].unique().apply(len).reset_index(drop=False)
    dict_num = dict(zip(df_num["Severity_visit"], df_num["ID"]))

    list_sev_vis = df_rna_count_meta[colsev].unique()
    list_sev_vis_sort = sorted(list_sev_vis, key=lambda x:x.split("_")[-1])
    list_sev_vis_shift = list_sev_vis_sort[1:] + [list_sev_vis_sort[0]]
    
    list_xticklabels = list()
    dict_name_sev = {"Severe": "Sev.", "Convalescent": "Conval."}
    dict_name_vis = {"First": "Acute", "Last": "Recov."}
    for key in list_sev_vis_shift:
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

def get_dict_palette(df_rna_count_meta, colsev="Severity_visit"):
    list_sev_vis = df_rna_count_meta[colsev].unique()
    list_sev_vis_sort = sorted(list_sev_vis, key=lambda x:x.split("_")[-1])
    list_sev_vis_shift = list_sev_vis_sort[1:] + [list_sev_vis_sort[0]]
    
    list_colors = list()
    for sev_vis in list_sev_vis_shift:
        if sev_vis.split("_")[0]=="Mild":
            list_colors.append("forestgreen")
        elif sev_vis.split("_")[0]=="Severe":
            list_colors.append("firebrick") 
        else: 
            list_colors.append("royalblue")
            
    dict_palette = dict(zip(list_sev_vis_shift, list_colors))

    return dict_palette

def get_dict_pos(df_rna_count_meta, colsev="Severity_visit", list_pos=[0, 0.5, 1.5, 2.0, 3.0]):
    list_sev_vis = df_rna_count_meta[colsev].unique()
    list_sev_vis_sort = sorted(list_sev_vis, key=lambda x:x.split("_")[-1])
    list_sev_vis_shift = list_sev_vis_sort[1:] + [list_sev_vis_sort[0]]
    list_num = list(range(len(list_sev_vis_shift)))
    dict_pos = dict(zip(list_sev_vis_shift, list_pos))
    
    return dict_pos

def get_dictionary_gene_pvalsig(df_gene_diffexp_sorted_first, target_pair = ("Mild_First", "Severe_First")):
    dict_gene_pvalsig = dict()
    list_all_pairs = list(map(lambda x: sorted(x), df_gene_diffexp_sorted_first["comp"][0]))
    idx_target = list_all_pairs.index(sorted(target_pair))
    genes = df_gene_diffexp_sorted_first["padj"].index
    for gene in genes:
        fdr = df_gene_diffexp_sorted_first["padj"].loc[gene][idx_target]
        if fdr < 0.05:
            dict_gene_pvalsig[gene] = True
        else:
            dict_gene_pvalsig[gene] = False
    
    return dict_gene_pvalsig

def main(path_drop_samples, path_exp, path_methyl_marker, path_meta):
    list_drop_samples = get_drop_samples(path_drop_samples)
    df_rna_count = read_rna_count_matrix(path_exp, list_drop_samples)
    list_gene_markers = get_list_gene_markers(path_methyl_marker)
    df_rna_count_gene_filtered = filter_rna_count_matrix(df_rna_count, list_gene_markers)
    df_rna_count_transposed = transpose_rna_count_matrix(df_rna_count_gene_filtered)
    df_rna_count_transposed_reset_idx = reset_idx_rna_count_matrix(df_rna_count_transposed)
    df_rna_count_meta = merge_rna_count_matrix_meta(df_rna_count_transposed_reset_idx, path_meta)
    
    return df_rna_count_meta

# %%
# Create dummy scatter points for legend
# legend_elements = [
#     mpatches.Patch(color='forestgreen', label='Mild'),
#     mpatches.Patch(color='firebrick', label='Severe'), 
#     mpatches.Patch(color='royalblue', label='Convalescent')  
# ]
# # Add the common legend to the rightmost side of the figure
# plt.legend(handles=legend_elements, 
#            loc='lower left', 
#            bbox_to_anchor=(1.01, 1.0),
#            title="Severity", 
#            fontsize=plt.rcParams["font.size"]+5, 
#            title_fontsize=plt.rcParams["font.size"]+7,
#            frameon=False)
# # Hide unused subplots if ncols * nrows > n_genes
# for j in range(i + 1, len(axes)):
#     fig.delaxes(axes[j])
# def get_top_2_idx(df_rna_count_meta, list_select_gene, colsev="Severity_group"):
#     dict_top_2_idx = defaultdict(list)
#     for gene in list_select_gene:
#         # Calculate the mean expression value for "Mild" severity
#         mean_mild = df_rna_count_meta[df_rna_count_meta[colsev] == "Mild"][gene].mean()

#         # Calculate the mean expression value for non-Mild severity
#         mean_others = df_rna_count_meta[df_rna_count_meta[colsev] != "Mild"][gene].mean()

#         # Determine whether to use nlargest(2) or nsmallest(2)
#         if mean_mild > mean_others:
#             top_2_idx = list(df_rna_count_meta[gene].nsmallest(2).index)
#         else:
#             top_2_idx = list(df_rna_count_meta[gene].nlargest(2).index)
        
#         # Directly assign the top_2_idx list to the dictionary
#         dict_top_2_idx[gene] = top_2_idx
    
#     return dict_top_2_idx

# # Loop through the two selected samples and add labels
# top_2_idx = dict_top_2_idx[gene]
# for idx in top_2_idx:
#     # Ensure idx is a scalar
#     if isinstance(idx, pd.Series):
#         idx = idx.values[0]
    
#     # Check if idx is valid
#     if idx in df_rna_count_meta.index:
#         x_max = df_rna_count_meta.loc[idx, "x_pos"]
#         y_max = df_rna_count_meta.loc[idx, gene]
        
#         # Ensure idx is a scalar and convert to string if necessary
#         if np.isscalar(idx):
#             label = str(df_rna_count_meta.index[idx])
#         else:
#             label = str(idx)

#         # Label each sample with text on the right side of the scatterplot
#         ax.text(x_max + 0.20, 
#                 y_max, 
#                 s=label, 
#                 verticalalignment='center', 
#                 fontsize=16, 
#                 color='k',
#                 fontweight='bold')
#     else:
#         print(f"Warning: Index {idx} not found in df_rna_count_meta")
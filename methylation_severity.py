# %%
import glob
import gzip
import os
import pickle
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection


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

def get_dict_deg_dmp_overlap_markers(file_deg_dmp_overlap, col_marker="Methyl", col_gene = "RNA"):
    from collections import defaultdict
    df_deg_dmp_overlap = pd.read_csv(file_deg_dmp_overlap, sep="\t")
    list_marker = df_deg_dmp_overlap[col_marker].to_list()
    list_genesym = df_deg_dmp_overlap[col_gene].to_list()
    dict_all_markers = defaultdict(list)
    for marker, genesym in zip(list_marker, list_genesym):
        marker_edit = ":".join(marker.split("_")[:-1])
        dict_all_markers[marker_edit].append(genesym)

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
            list_colors.append("grey")
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
    df_beta_methyl_reset_idx = df_beta_methyl.reset_index(drop=False)
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

def find_meandiff(a, b):
    deltamean = np.mean(b) - np.mean(a)

    return deltamean

def make_dataframe_stat_test(df_marker, list_target_marker, colsev="Severity_visit"):
    dict_marker_stat = dict()
    for marker in list_target_marker:
        df_group_visit = df_marker.groupby(colsev)[marker].apply(np.array).reset_index(drop=False)
        list_groups = df_group_visit[colsev].values
        list_values = df_group_visit[marker].values
        list_comb_group = list(combinations(list_groups[::-1], r=2))
        list_comb_marker_val = list(combinations(list_values[::-1], r=2))   
        list_meandiff = list(map(lambda x: find_meandiff(*x), list_comb_marker_val))
        list_stat = list(map(lambda x: ttest_ind(*x, equal_var=False)[0], list_comb_marker_val))
        list_pval = list(map(lambda x: ttest_ind(*x, equal_var=False)[1], list_comb_marker_val))
        list_sig, list_padj = fdrcorrection(list_pval)
        dict_stat_res = {"comp": list_comb_group, "diff": list_meandiff, "stat": list_stat, "pval": list_pval, "padj": list_padj, "is_sig": list_sig}
        dict_marker_stat[marker] = dict_stat_res
    
    df_marker_stat = pd.DataFrame.from_dict(dict_marker_stat).T
    df_marker_stat_sorted = df_marker_stat.sort_values(by=["pval"], ascending=True)
    
    return df_marker_stat_sorted

def get_dictionary_methyl_pvalsig(df_marker_stat_sorted, target_pair = ("Mild_First", "Severe_First")):
    dict_gene_pvalsig = dict()
    list_all_pairs = list(map(lambda x: sorted(x), df_marker_stat_sorted["comp"][0]))
    idx_target = list_all_pairs.index(sorted(target_pair))
    genes = df_marker_stat_sorted["padj"].index
    for gene in genes:
        fdr = df_marker_stat_sorted["padj"].loc[gene][idx_target]
        if fdr < 0.05:
            dict_gene_pvalsig[gene] = True
        else:
            dict_gene_pvalsig[gene] = False
    
    return dict_gene_pvalsig

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

def main(path_sev_info, dir_methylcpgmin, infilenamehyper, infilenamehypo):
    list_files_methylcpgmin = get_list_files_methylcpgmin(dir_methylcpgmin)
    list_methyl_sample = get_list_methyl_sample(list_files_methylcpgmin)
    df_beta_all_hyper = get_dataframe_methylation_beta_samples(path_sev_info, infilenamehyper, list_methyl_sample)
    df_beta_all_hypo = get_dataframe_methylation_beta_samples(path_sev_info, infilenamehypo, list_methyl_sample)

    return df_beta_all_hyper, df_beta_all_hypo

# %%
import gzip
import math
import os
import pickle
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from joblib import Parallel, delayed
from matplotlib import pyplot as plt

# %%
def get_list_files_deg_loo(path_deg_loo_format):
    list_files_loo = glob(path_deg_loo_format)
    
    return list_files_loo

def get_dict_genes_deg_sig(list_files_loo, fc_cutoff, padj_cutoff):
    dict_genes_deg_significant_of_all_loo = dict()
    dict_genes_deg_significant_of_all_loo_sample_wo = dict()
    for deg_file in list_files_loo:
        table_deg_tmp = pd.read_csv(deg_file, sep = '\t')
        table_deg_tmp_sig = table_deg_tmp[np.logical_and(table_deg_tmp["log2FoldChange"].abs() > fc_cutoff, table_deg_tmp["padj"] < padj_cutoff)]
        for geneid, fc in zip(table_deg_tmp_sig["ID"].to_list(), table_deg_tmp_sig["log2FoldChange"].to_list()):
            if dict_genes_deg_significant_of_all_loo.get(geneid) == None:
                dict_genes_deg_significant_of_all_loo[geneid] = list()
                dict_genes_deg_significant_of_all_loo_sample_wo[geneid] = list()
            dict_genes_deg_significant_of_all_loo[geneid].append(fc)
            dict_genes_deg_significant_of_all_loo_sample_wo[geneid].append(deg_file.split('_')[-2])

    for geneid in list(dict_genes_deg_significant_of_all_loo.keys()):
        directions = list(map(lambda x : x < 0, dict_genes_deg_significant_of_all_loo[geneid]))
        is_consistent = sum(directions) == 0 or sum(directions) == len(directions)
        if not is_consistent:
            print(geneid)
            dict_genes_deg_significant_of_all_loo.pop(geneid)
    
    return dict_genes_deg_significant_of_all_loo, dict_genes_deg_significant_of_all_loo_sample_wo

def get_dmp_correlated_rna_fdr_sig(path_dmp_correlated_rna, dmp_name, fold, col_p, padj_cutoff):
    table_dmp_correlated_rna = pd.read_csv(path_dmp_correlated_rna.format(dmp_name, fold), sep = '\t')
    table_dmp_correlated_rna_fdr_sig = table_dmp_correlated_rna[table_dmp_correlated_rna[col_p] < padj_cutoff]
    
    return table_dmp_correlated_rna_fdr_sig

def filter_dict_genes_loo_overlap_threshold(dict_genes_deg_significant_of_all_loo, loo_overlap_cutoffs):
    count_deg_gene = dict(zip(dict_genes_deg_significant_of_all_loo.keys(), list(map(lambda g : len(dict_genes_deg_significant_of_all_loo[g]), dict_genes_deg_significant_of_all_loo.keys()))))

    dict_genes_over_loo_overlap_threshold = dict()
    for overlap_cutoff in loo_overlap_cutoffs:
        dict_genes_over_loo_overlap_threshold[overlap_cutoff] = list(filter(lambda gene : count_deg_gene[gene] >= overlap_cutoff, count_deg_gene.keys()))
    
    dict_genes_equal_loo_overlap_threshold = dict()
    for overlap_cutoff in loo_overlap_cutoffs:
        dict_genes_equal_loo_overlap_threshold[overlap_cutoff] = list(filter(lambda gene : count_deg_gene[gene] == overlap_cutoff, count_deg_gene.keys()))


    return dict_genes_over_loo_overlap_threshold, dict_genes_equal_loo_overlap_threshold

def filter_dict_genes_dmp_deg_corr_threshold(table_dmp_correlated_rna_fdr_sig, corr_cutoffs):
    # Methyl ~ RNA
    dict_genes_over_corr_threshold = dict()
    for corr_cutoff in corr_cutoffs:
        table_dmp_correlated_rna_fdr_sig_corr_sig = table_dmp_correlated_rna_fdr_sig[table_dmp_correlated_rna_fdr_sig["corr_rho"].abs() > corr_cutoff]
        dict_genes_over_corr_threshold[corr_cutoff] = list(table_dmp_correlated_rna_fdr_sig_corr_sig["RNA"].unique())
    
    return dict_genes_over_corr_threshold

def dump_pickle(data, outfile):
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    with gzip.open(outfile, 'wb') as f:
        pickle.dump(data, f)

def save_dict_deg_loo(dict_common_deg_loo, outdir, loop_overlap_cutoff):
    outfile = os.path.join(outdir, "DEG_LOO", f"dictionary_deg_loo_common_{loop_overlap_cutoff}folds.pk.gz")
    dict_deg_loo_overlap_threshold = {marker: log2fc for marker, log2fc in dict_common_deg_loo.items() if len(log2fc) >= loop_overlap_cutoff}
    dump_pickle(dict_deg_loo_overlap_threshold, outfile)
    
def save_list_genes_over_loo_overlap_threshold(dict_genes_over_loo_overlap_threshold, loo_overlap_threshold, outdir):
    outfile = os.path.join(outdir, "DEG_LOO", f"list_genes_over_loo_overlap_threshold_{str(loo_overlap_threshold)}.txt")
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    # if not os.path.exists(outfile):
    with open(outfile, mode="w") as fw:
        list_genes = dict_genes_over_loo_overlap_threshold.get(loo_overlap_threshold)
        print(f"Number of DEG-LOO Genes: {len(list_genes)}")
        for gene in list_genes:
            fw.write(gene+"\n")
        
def plot_deg_overlap_loo(dict_genes_over_loo_overlap_threshold, loo_overlap_cutoffs, outdir, outfigname):
    # LOO Stats
    outfile = os.path.join(outdir, f"plot_deg_overlap_loo_{outfigname}.png")
    # if not os.path.exists(outfile):
    plt.rcParams["font.size"] = 13
    plt.bar(loo_overlap_cutoffs, list(map(lambda cf : len(dict_genes_over_loo_overlap_threshold[cf]), loo_overlap_cutoffs)), color = "gray")
    for loo_cf in loo_overlap_cutoffs:
        cnt_genes = len(dict_genes_over_loo_overlap_threshold[loo_cf])
        plt.text(loo_cf, cnt_genes, cnt_genes, ha = "center", va = "bottom")
    plt.xlabel("Number of Overlaps on DEG-LOO\n(LOO on Severe,  samples)")
    plt.ylabel("Number of Genes")
    
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.show()
    plt.close()

def plot_sample_specific_deg_discovery(dict_genes_deg_significant_of_all_loo_sample_wo, list_files_loo, outdir, outfigname):
    outfile = os.path.join(outdir, f"plot_sample_specific_deg_discovery_{outfigname}.png")
    # if not os.path.exists(outfile):
    total_loo_samples = list(map(lambda x : x.split('_')[-2], list_files_loo))
    table_loo_sample_specific = pd.DataFrame(columns = total_loo_samples, index = list(dict_genes_deg_significant_of_all_loo_sample_wo.keys())).fillna(0)
    for gene, list_sample in dict_genes_deg_significant_of_all_loo_sample_wo.items():
        table_loo_sample_specific.loc[gene, list_sample] = 1
    table_loo_sample_specific["discovery"] = table_loo_sample_specific.apply(sum, axis = 1)
    table_loo_sample_specific = table_loo_sample_specific.sort_values(by = "discovery", ascending=False)
    dict_gene_first_cnt = dict()
    for loo_cnt in sorted(table_loo_sample_specific["discovery"].unique()):
        ind_first_loo_cnt = list(table_loo_sample_specific["discovery"]).index(loo_cnt)
        dict_gene_first_cnt[loo_cnt] = list(table_loo_sample_specific.index)[ind_first_loo_cnt]

    table_loo_sample_specific = table_loo_sample_specific.drop(columns = ["discovery"])
    plt.figure(figsize = (5, 12))
    sns.heatmap(table_loo_sample_specific,cbar_kws = {"label":"Whether marker discovered (True/False)"})
    plt.yticks([])
    for loo_cnt, gene_first in dict_gene_first_cnt.items():
        if loo_cnt == len(total_loo_samples):continue
        plt.axhline(list(table_loo_sample_specific.index).index(gene_first), linestyle = "--", linewidth = 1, color = 'r')
    plt.xlabel("Excluded Sample")
    plt.ylabel("DEG significant genes (Sorted by Number of LOO-Overlaps)")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.show()
    plt.close()

def plot_heatmap_dmp_deg_overlap(dict_genes_over_loo_overlap_threshold, dict_genes_over_corr_threshold, loo_overlap_cutoffs, corr_cutoffs, outdir, outfigname):
    # Gene Count Heatmap
    table_gene_overlap_counter_per_thresholds = pd.DataFrame(columns = list(loo_overlap_cutoffs), index = reversed(list(corr_cutoffs)))

    for loo_cf in loo_overlap_cutoffs:
        for corr_cf in corr_cutoffs:
            list_genes_over_loo_cf = dict_genes_over_loo_overlap_threshold[loo_cf]
            list_genes_over_corr_cf = dict_genes_over_corr_threshold[corr_cf]
            list_genes_overlap_between_two = list(set(list_genes_over_corr_cf) & set(list_genes_over_loo_cf))
            table_gene_overlap_counter_per_thresholds.loc[corr_cf, loo_cf] = len(list_genes_overlap_between_two)
            
    sns.heatmap(data = table_gene_overlap_counter_per_thresholds.astype(float), annot = True, fmt = ".0f")
    plt.xlabel("Number of DEG overlaps during 9 time of LOO")
    plt.ylabel("Correlation cutoffs between\nLOO-common DMSs and Gene expression")
    plt.title("Number of Genes overlap\n(DMSs-Gene Correlation & LOO DEG)")
    outfile = os.path.join(outdir, f"plot_heatmap_dmp_deg_overlap_{outfigname}.png")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.show()
    plt.close()

def save_table_dmp_deg_corr(table_dmp_correlated_rna_fdr_sig, dict_genes_over_corr_threshold, dict_genes_over_loo_overlap_threshold, outfileformat, path_out_marker, outfilename, fold, corr_cutoff, fc_cutoff, loo_cnt_cutoff):
    table_rna_methylcorr_sig = table_dmp_correlated_rna_fdr_sig[table_dmp_correlated_rna_fdr_sig["corr_rho"].abs() > corr_cutoff]
    table_rna_methylcorr_sig = table_rna_methylcorr_sig[table_rna_methylcorr_sig["RNA"].isin(set(dict_genes_over_corr_threshold[corr_cutoff]) & set(dict_genes_over_loo_overlap_threshold[loo_cnt_cutoff]))]
    table_rna_methylcorr_sig["-log10(fdr)"] = table_rna_methylcorr_sig["corr_fdr"].apply(lambda x : -math.log10(x+1))
    list_rna_methylcorr_sig = table_rna_methylcorr_sig["RNA"].to_list()
    with open(path_out_marker, mode="w") as fw:
        for rna in list_rna_methylcorr_sig:
            fw.write(rna + "\n")
            
    table_rna_methylcorr_sig.to_csv(outfileformat.format(outfilename, fold, corr_cutoff, fc_cutoff, loo_cnt_cutoff), sep = '\t', index = False)

def main(path_deg_loo_format, path_dmp_correlated_rna, fc_cutoff, padj_cutoff, col_p, corr_cutoffs, loo_overlap_cutoffs, outdir, outfileformat, path_out_marker, covariate, fold):
    list_files_loo = get_list_files_deg_loo(path_deg_loo_format)
    dict_genes_deg_significant_of_all_loo, dict_genes_deg_significant_of_all_loo_sample_wo = get_dict_genes_deg_sig(list_files_loo, fc_cutoff, padj_cutoff)
    for loop_overlap_cutoff in loo_overlap_cutoffs:
        save_dict_deg_loo(dict_genes_deg_significant_of_all_loo, outdir, loop_overlap_cutoff)
        
    dict_genes_over_loo_overlap_threshold, dict_genes_equal_loo_overlap_threshold = filter_dict_genes_loo_overlap_threshold(dict_genes_deg_significant_of_all_loo, loo_overlap_cutoffs)
    for loop_overlap_cutoff in loo_overlap_cutoffs:
        save_list_genes_over_loo_overlap_threshold(dict_genes_over_loo_overlap_threshold, loop_overlap_cutoff, outdir)
    plot_deg_overlap_loo(dict_genes_equal_loo_overlap_threshold, loo_overlap_cutoffs, outdir=outdir, outfigname="DEG_LOO")
    plot_sample_specific_deg_discovery(dict_genes_deg_significant_of_all_loo_sample_wo, list_files_loo, outdir=outdir, outfigname="DEG_LOO")
    
    table_dmp_correlated_rna_fdr_sig = get_dmp_correlated_rna_fdr_sig(path_dmp_correlated_rna, covariate, fold, col_p, padj_cutoff)
    dict_genes_over_corr_threshold = filter_dict_genes_dmp_deg_corr_threshold(table_dmp_correlated_rna_fdr_sig, corr_cutoffs)
    plot_heatmap_dmp_deg_overlap(dict_genes_over_loo_overlap_threshold, dict_genes_over_corr_threshold, loo_overlap_cutoffs, corr_cutoffs, outdir=outdir, outfigname=f"{covariate}_{fold}folds_or_more")
    save_table_dmp_deg_corr(table_dmp_correlated_rna_fdr_sig, dict_genes_over_corr_threshold, dict_genes_over_loo_overlap_threshold, outfileformat, path_out_marker.format(covariate), outfilename=covariate, fold=fold, corr_cutoff = 0.5, fc_cutoff = 1.3, loo_cnt_cutoff = 7)

# %%
fold = 9
col_p = "corr_fdr"
padj_cutoff = 0.05
fc_cutoff = 1.3
corr_cutoffs = np.linspace(0.10, 0.70, 7)
loo_overlap_cutoffs = range(1, 10, 1)
WORKDIR = str(Path(__file__).parents[3])
path_deg_loo_format = os.path.join(WORKDIR, "Results/5_deg/Visit1_Severe__Visit1_Mild_removed_loo_*_20240327.tsv")
path_dmp_correlated_rna = os.path.join(WORKDIR,"Results/11_dmp/Yes_LOO/Covariate_{0}/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_LOO_common_across_{1}folds_or_more.cis.distance_cf.1000000.20240326.sig_dmps.tsv")
outdir = os.path.join(WORKDIR, "Results/11_dmp/Yes_LOO")
os.makedirs(outdir, exist_ok=True)
outfileformat = os.path.join(outdir, "Covariate_{0}/Methyl_RNA_Correlation.Filtered.DMP_Covariate_{0}_Yes_LOO_common_across_{1}folds_or_more.cis.distance_cf.1000000.Overlap_DEG.corr_{2}.log2fc_{3}.loo_{4}.20240326.tsv")
path_out_marker =  os.path.join(outdir, "Covariate_{0}/list_new_gene_markers.txt")
list_covariates = ["Sex_Age"]
# list_covariates = ["Sex_Age_CCIScore", "Sex_Age_CellType_2samples_CellPropTotal_Not_100_Included_New", "Sex_Age_CellType_2samples_CellPropTotal_Not_100_Included_Smoking", "Sex_Age_CellType_2samples_CellPropTotal_Not_100_Included_Smoking_CCI", "Sex_Age_Smoking_Status", "Sex_Age_Smoking_Status_CCIScore"]

n_jobs = len(list_covariates)
with Parallel(n_jobs=n_jobs) as parallel:
    parallel(delayed(main)(path_deg_loo_format, path_dmp_correlated_rna, fc_cutoff, padj_cutoff, col_p, corr_cutoffs, loo_overlap_cutoffs, outdir, outfileformat, path_out_marker, covariate, fold) for covariate in list_covariates)

# %%

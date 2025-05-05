#%%
import math
from pathlib import Path

import pandas as pd


# %%
def main(path_corr, path_gene_info, path_save_cis_format, path_save_trans_format, list_distance_cutoff = list(), col_geneid = "gene_id", col_genename = "gene_name", col_chr = "chr", col_start = "start", col_end = "end"):
    table_corr = pd.read_csv(path_corr, sep = '\t')
    table_gene_info = pd.read_csv(path_gene_info, sep = '\t')
    
    table_corr["Methyl_chr"] = table_corr["Methyl"].apply(lambda x : x.split('_')[0])
    table_corr["Methyl_pos"] = table_corr["Methyl"].apply(lambda x : int(x.split('_')[1]))
    
    table_gene_info["RNA_name"] = table_gene_info.apply(lambda x : f"{x[col_geneid]}_{x[col_genename]}", axis = 1)
    
    table_corr_rna_info = get_gene_position_to_correlation_table(table_gene_info, table_corr, col_chr, col_start, col_end)
    dict_distance_cf_to_table = filter_trans_cis_table(table_corr_rna_info, list_distance_cutoff)
    for cutoff in dict_distance_cf_to_table.keys():
        dict_distance_cf_to_table[cutoff]["Cis"].to_csv(path_save_cis_format, sep = '\t', index = False)
        # dict_distance_cf_to_table[cutoff]["Trans"].to_csv(path_save_trans_format, sep = '\t', index = False)
   
def get_gene_position_to_correlation_table(table_gene_info, table_corr, col_chr, col_start, col_end):
    dict_gene_to_chr = dict(zip(table_gene_info["RNA_name"], table_gene_info[col_chr]))
    
    table_gene_info["TSS"] = table_gene_info.apply(lambda x : x[col_start] if x["strand"] == '+' else x[col_end], axis = 1)
    dict_gene_to_tss = dict(zip(table_gene_info["RNA_name"], table_gene_info["TSS"]))
    
    table_corr["RNA_chr"] = table_corr["RNA"].apply(dict_gene_to_chr.__getitem__)
    table_corr["RNA_TSS"] = table_corr["RNA"].apply(dict_gene_to_tss.__getitem__)
    
    return table_corr

def filter_trans_cis_table(table_corr, distance_cutoff):
    table_trans_diff_chr = table_corr[table_corr["RNA_chr"] != table_corr["Methyl_chr"]]
    table_same_chr = table_corr[table_corr["RNA_chr"] == table_corr["Methyl_chr"]]
    
    dict_cf_to_table = dict()
    for cutoff in distance_cutoff:
        if cutoff == None:
            cutoff = math.inf
        pos_diff = abs(table_corr["Methyl_pos"] - table_corr["RNA_TSS"])
        
        same_chr_is_in_distance = pos_diff <= cutoff
        same_chr_out_distance = pos_diff > cutoff
        table_same_chr_cis = table_same_chr[same_chr_is_in_distance]
        table_same_chr_trans = table_same_chr[same_chr_out_distance]
        dict_cf_to_table[cutoff] = {
            "Cis" : table_same_chr_cis,
            "Trans" : pd.concat([table_same_chr_trans, table_trans_diff_chr])
        }
    return dict_cf_to_table

#%%
if __name__ == "__main__":
    WORKDIR = str(Path(__file__).parents[3])
    path_corr = f"{WORKDIR}/Resources/Data/Methylation/MethylKitTable_Incl_EHR_Deconvoluton/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Covariate_Sex_Age_CellDeconv_Smoking_Medication_Comorbidity_Yes_LOO_common_across_4folds_or_more.tsv"
    path_gene_info = f"{WORKDIR}/Resources/References/GTF_Gene_Parsing/GTF.genecode.v42.chr_patch_hapl_scaff.gene_only.tsv"
    distance_cutoff = [1000000]
    path_save_cis_format = f"{WORKDIR}/Resources/Data/Methylation/MethylKitTable_Incl_EHR_Deconvoluton/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Covariate_Sex_Age_CellDeconv_Smoking_Medication_Comorbidity_Yes_LOO_common_across_4folds_or_more.cis.distance_cf.{distance_cutoff[0]}.tsv"
    path_save_trans_format = f"{WORKDIR}/Resources/Data/Methylation/MethylKitTable_Incl_EHR_Deconvoluton/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Covariate_Sex_Age_CellDeconv_Smoking_Medication_Comorbidity_Yes_LOO_common_across_4folds_or_more.trans.distance_cf.{distance_cutoff[0]}.tsv"

    
    # path_corr = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit3_4.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_VisitLast.LOO_common.20240326.tsv"
    # path_gene_info = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/GTF_Gene_Parsing/GTF.genecode.v42.chr_patch_hapl_scaff.gene_only.tsv"
    
    # distance_cutoff = [None, 1000000]
    # path_save_cis_format = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit3_4.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_VisitLast.LOO_common.cis.distance_cf.{0}.20240326.tsv"
    # path_save_trans_format = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit3_4.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_VisitLast.LOO_common.trans.distance_cf.{0}.20240326.tsv"
    
    main(path_corr, path_gene_info, path_save_cis_format, path_save_trans_format, list_distance_cutoff=distance_cutoff)

# %%
import matplotlib.pyplot as plt
import pandas as pd

WORKDIR = str(Path(__file__).parents[3])
table_corr = pd.read_csv(path_save_cis_format, sep="\t")
pos_diff = abs(table_corr["Methyl_pos"] - table_corr["RNA_TSS"])
plt.hist(pos_diff, color="grey")
plt.xlabel("positional difference (bp)")
plt.ylabel("Count")
plt.show()
plt.close()
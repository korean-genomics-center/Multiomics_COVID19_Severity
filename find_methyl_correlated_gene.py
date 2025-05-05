#%%
import math
import os
import warnings
from pathlib import Path

import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from joblib import Parallel, delayed
from scipy.stats import chi2, spearmanr
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import fdrcorrection

warnings.filterwarnings("ignore")


def main(path_methyl, path_exp, n_job, path_save, cols_methyl_feat = ["chr", "start", "end"], cols_exp_feat = ["index"]):
    table_methyl = pd.read_csv(path_methyl, sep = '\t')
    table_exp = pd.read_csv(path_exp, sep = '\t')
    
    table_methyl_indexed = make_feature_index(table_methyl, cols_methyl_feat)
    table_exp_indexed = make_feature_index(table_exp, cols_exp_feat)

    list_sample_methyl = list(table_methyl_indexed.columns)
    list_sample_exp = list(table_exp_indexed.columns)
    list_overlap_sample = list(set(list_sample_methyl) & set(list_sample_exp))
    print(len(list_overlap_sample), flush = True)
    print(list_overlap_sample, flush = True)
    
    list_methyl_indexes_split_by_job = split_data_by_n(list(table_methyl_indexed.index), n_job)
    list_table_methyl_indexed_by_job = list(map(lambda list_feat: table_methyl_indexed.loc[list_feat, :], list_methyl_indexes_split_by_job))
    
    with Parallel(n_jobs=n_job) as parallel:
        list_results = parallel(delayed(calculate_methyl_expression_relationship)(table_methyl_part, table_exp_indexed, list_overlap_sample) for table_methyl_part in list_table_methyl_indexed_by_job)
    table_tot_relation = pd.concat(list_results)
    
    table_tot_relation_corr_adj = add_fdr_corrected_pvalue(table_tot_relation, "corr_p", "corr_fdr")
    table_tot_relation_all_adj = add_fdr_corrected_pvalue(table_tot_relation_corr_adj, "lr_p", "lr_fdr")
    
    os.makedirs(os.path.dirname(path_save), exist_ok = True)
    table_tot_relation_all_adj.to_csv(path_save, sep = '\t', index = False)

def make_feature_index(table, cols_feat):
    if "index" in cols_feat:
        table["index"] = table.index
    table["featname"] = table[cols_feat].apply(lambda row : '_'.join(list(map(str, row))), axis = 1)
    table = table.drop(columns = cols_feat)
    table = table.set_index("featname", drop = True)
    return table

def split_data_by_n(list_values, n_split):
    num_data_per_split = math.floor(len(list_values)/n_split)
    num_leftover = len(list_values) % n_split
    list_values_split = list()
    ind_value_end = 0
    for i in range(n_split):
        n_values_append = num_data_per_split
        if i < num_leftover:
            n_values_append += 1
        list_values_part = list_values[ind_value_end:ind_value_end+n_values_append]
        ind_value_end += n_values_append
        list_values_split.append(list_values_part)
    return list_values_split

def calculate_methyl_expression_relationship(table_methyl, table_exp, list_samples_check):
    list_feat_methyl = list(table_methyl.index)
    list_feat_exp = list(table_exp.index)
    
    dict_relation = dict()
    
    for feat_methyl in list_feat_methyl:
        dict_relation[feat_methyl] = dict()
        dict_sample_to_methyl = table_methyl.loc[feat_methyl, :].to_dict()
        for feat_exp in list_feat_exp:
            dict_sample_to_exp = table_exp.loc[feat_exp, :].to_dict()
            
            corr_rho, corr_p = calculate_spearman_correlation_from_two_dictionaries(dict_sample_to_methyl, dict_sample_to_exp, list_samples_check)
            lr_beta, lr_p = fit_linear_regression_from_two_dictionaries(dict_sample_to_methyl, dict_sample_to_exp, list_samples_check)
            dict_relation[feat_methyl][feat_exp] = {
                "corr_rho":corr_rho,
                "corr_p":corr_p,
                "lr_beta":lr_beta,
                "lr_p":lr_p
            }
    table_relation = organize_methyl_expression_relation_result(dict_relation)
    return table_relation
            
def calculate_spearman_correlation_from_two_dictionaries(dict1, dict2, list_check):
    list_values1 = list(map(dict1.__getitem__, list_check))
    list_values2 = list(map(dict2.__getitem__, list_check))
    
    corr, p = spearmanr(list_values1, list_values2)
    return corr, p

def fit_linear_regression_from_two_dictionaries(dict_x, dict_y, list_check):
    values_x = list(map(dict_x.__getitem__, list_check))
    values_y = list(map(dict_y.__getitem__, list_check))
    
    try:
        lr_model = sm.GLM(values_y, sm.add_constant(values_x))
        fit_result = lr_model.fit()
        
        cons = fit_result.params[0]
        beta = fit_result.params[1]
        p = chi2.sf(fit_result.null_deviance - fit_result.deviance, 1)
    except Exception:
        beta = pd.NA
        p = pd.NA
    return beta, p

def organize_methyl_expression_relation_result(dict_relation):
    list_methyl_feat = list()
    list_rna_feat = list()
    list_corr_rho = list()
    list_corr_p = list()
    list_lr_beta = list()
    list_lr_p = list()
    for feat_methyl in dict_relation.keys():
        for feat_rna in dict_relation[feat_methyl].keys():
            corr_rho = dict_relation[feat_methyl][feat_rna]["corr_rho"]
            corr_p = dict_relation[feat_methyl][feat_rna]["corr_p"]
            lr_beta = dict_relation[feat_methyl][feat_rna]["lr_beta"]
            lr_p = dict_relation[feat_methyl][feat_rna]["lr_p"]
            list_methyl_feat.append(feat_methyl)
            list_rna_feat.append(feat_rna)
            list_corr_rho.append(corr_rho)
            list_corr_p.append(corr_p)
            list_lr_beta.append(lr_beta)
            list_lr_p.append(lr_p)
    table_res = pd.DataFrame({
        "Methyl":list_methyl_feat,
        "RNA" : list_rna_feat,
        "corr_rho" : list_corr_rho,
        "corr_p" : list_corr_p,
        "lr_beta" : list_lr_beta,
        "lr_p" : list_lr_p
    })
    return table_res

def add_fdr_corrected_pvalue(table, col_p, col_fdr_add):
    table_p_na = table[table[col_p].isna()]
    table_p_notna = table[table[col_p].notna()]
    
    _, fdr = fdrcorrection(table_p_notna[col_p].to_list())
    table_p_notna[col_fdr_add] = fdr
    table_p_na[col_fdr_add] = pd.NA
    
    return pd.concat([table_p_notna, table_p_na])

#%%
if __name__ == "__main__":
    fold = 9
    WORKDIR = str(Path(os.path.abspath(__file__)).parents[3])
    path_methyl = f"{WORKDIR}/Resources/Data/Methylation/MethylKitTable_Incl_EHR_Deconvoluton/MethylCpGTable.Control.Mild.Case.Severe.filtered.Covariate_Sex_Age_Yes_LOO_common_across_{fold}folds_or_more.tsv"
    path_exp = f"{WORKDIR}/Results/5_deg/Visit1_Severe__Visit1_Mild_20240327.tsv.normcount"
    path_save = f"{WORKDIR}/Resources/Data/Methylation/MethylKitTable_Incl_EHR_Deconvoluton/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Covariate_Sex_Age_Yes_LOO_common_across_{fold}folds_or_more.tsv"
    # path_methyl = f"{WORKDIR}/Resources/Data/Methylation/Infectomics.Copy_From_HjRyu/MethylCpGTable.Control.Mild.Case.Severe.Filtered.DMP.Hyper_Hypo.Sev_vs_Mild.VisitLast.LOO_common.tsv"
    # path_exp = "{WORKDIR}/Results/5_deg/Visit4__Visit3_DEG_20230925.tsv.normcount"
    # path_save = f"{WORKDIR}/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_VisitLast.LOO_common.20240326.tsv"
    n_job = 50
    
    main(path_methyl, path_exp, n_job, path_save)
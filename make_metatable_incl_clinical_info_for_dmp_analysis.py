# %%
import glob
import json
import os
import re
import unicodedata
import warnings
from collections import Counter
from itertools import chain
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

import clinical

warnings.filterwarnings("ignore")

# %%
def main(path_ehr, path_deconv, path_drop_samples, path_rename, path_comorbidity, dir_metadata, pattern_metadata, dir_metadata_incl_clin, flag_celltype, flag_deconv, flag_smoking, flag_medication, flag_comorbidity):
    list_drop_samples = get_drop_samples(path_drop_samples)
    parse_clin = clinical.Parse(path_ehr, list_drop_samples, path_rename)
    df_crf = parse_clin.give_parsed_dataframe()
    df_crf_inf = df_crf[df_crf["Sample_ID"].str.startswith("C19-C")]
    df_crf_inf_v1 = df_crf_inf[df_crf_inf["Sample_ID"].str.endswith("V1")]
    df_deconv = pd.read_csv(path_deconv, sep=",")
    df_crf_inf_v1_copy = df_crf_inf_v1.copy()
    covariates = list()
    if flag_celltype:
        list_celltype = ["Neutrophil", "Lymphocytes", "Monocytes", "Eosinophils", "Basophils"]
        covariates += list_celltype
        df_celltype_info_added = add_celltype_prop_df(df_crf_inf_v1_copy, list_celltype=list_celltype)
        df_crf_inf_v1_copy = df_celltype_info_added
    if flag_deconv:
        df_deconv_info_added = add_deconv_prop_df(df_crf_inf_v1_copy, df_deconv)
        list_celldeconv = list(filter(lambda x: str(x).startswith("Blood-"), list(df_deconv_info_added.columns)))
        covariates += list_celldeconv
        df_crf_inf_v1_copy = df_deconv_info_added
    if flag_smoking:
        covariates += ["SmokingStatus"]
        df_crf_smoke_info_added = add_smoking_info_df(df_crf_inf_v1_copy, col_smoke="SmokingStatus")
        df_crf_inf_v1_copy = df_crf_smoke_info_added
    if flag_medication:
        # covariates += ["MedicationUse"]
        # df_crf_med_info_added = add_mediciation_info_df(df_crf_inf_v1_copy, col_med="MedicationUse")
        covariates += ["SteroidUse", "AntiviralUse"]
        df_crf_med_info_added = add_mediciation_info_df(df_crf_inf_v1_copy)
        df_crf_inf_v1_copy = df_crf_med_info_added
    if flag_comorbidity:
        # dict_covariates = get_dict_covariates(path_comorbidity)
        # list_comorbidities = list(dict_covariates.keys())
        # covariates += list_comorbidities
        # df_crf_inf_v1_covariate_added = add_comorbidity_info_df(df_crf_inf_v1_copy, dict_covariates, col_comorbid="Comorbidity")
        # df_crf_inf_v1_copy = df_crf_inf_v1_covariate_added
        covariates += ["CCI_score"]
        
    df_crf_inf_v1_covariates = df_crf_inf_v1_copy[["Sample_ID", "Severity_group"] + covariates]

    if flag_deconv:
        list_celltype_ori = list(filter(lambda x: str(x).startswith("Blood-"), list(df_crf_inf_v1_covariates.columns)))
        list_celltype_new = list(map(lambda x: x.replace("-", "_").replace("+", "_"), list_celltype_ori))
        dict_celltype_rename = dict(zip(list_celltype_ori, list_celltype_new))
        df_crf_inf_v1_covariates = df_crf_inf_v1_covariates.rename(columns=dict_celltype_rename)
    
    df_crf_inf_v1_covariates.to_csv(os.path.join(dir_metadata_incl_clin, "COVID19_Covariate_Info.txt"), sep="\t", index=False)
    save_stats_proportion_covariates(df_crf_inf_v1_covariates, covariates)
    # plot_proportion_comparison(list_cov_prop, dir_metadata_incl_clin, list_mild_count, list_mild_nobs, list_sev_count, list_sev_nobs)

    
    df_crf_inf_v1_covariates_dropna = df_crf_inf_v1_covariates.dropna(axis=0)
    df_crf_inf_v1_covariates_dropna = df_crf_inf_v1_covariates_dropna.drop(columns=["Severity_group"])
    list_samples_ehr = df_crf_inf_v1_covariates_dropna["Sample_ID"].to_list()
    list_metadata = glob.glob(f"{dir_metadata}/**/{pattern_metadata}", recursive=True)
    list_samples_iter = list(map(lambda x: os.path.basename(x).split("_")[-1].replace(".tsv", "") + "-V1", list_metadata))
    for metadata, sample_iter in zip(list_metadata, list_samples_iter):
        if sample_iter in list_samples_ehr:
            filename = os.path.basename(metadata)
            df_metadata = pd.read_csv(metadata, sep="\t")
            df_merged = df_metadata.merge(df_crf_inf_v1_covariates_dropna, how="inner", on="Sample_ID")
            df_merged.to_csv(os.path.join(dir_metadata_incl_clin, filename), sep="\t", index=False)

def get_drop_samples(path_drop_samples):
    with open(path_drop_samples, mode="r") as frmv:
        list_drop_samples = [x.rstrip() for x in frmv.readlines()]
    
    return list_drop_samples

def add_celltype_prop_df(df_crf, list_celltype = ["Neutrophil", "Lymphocytes", "Monocytes", "Eosinophils", "Basophils"]):
    df_crf_nona_celltype_info = df_crf.copy()
    for celltype in list_celltype:
        df_crf_nona_celltype_info.loc[df_crf_nona_celltype_info[celltype].isna(),celltype] = np.nan
    # # hardcoding wrong value 102 --> 10.2 for the log
    df_crf_nona_celltype_info = df_crf_nona_celltype_info.set_index("Sample_ID")
    df_crf_nona_celltype_info.loc["C19-C005-V1", "Monocytes"] = 10.2
    df_crf_nona_celltype_info = df_crf_nona_celltype_info.reset_index(drop=False)
    df_crf_nona_celltype_info["CellPropTotal"] = df_crf_nona_celltype_info[list_celltype].sum(axis=1)

    # list_excl_sample_ind = list()
    # for ind, row in df_crf_nona_celltype_info[list_celltype + ["CellPropTotal"]].iterrows():
    #     dict_row = dict(row)
    #     if round(dict_row.get("CellPropTotal"), 1) != 100:
    #         list_excl_sample_ind.append(ind)

    # df_celltype_info_added = df_crf_nona_celltype_info.drop(list_excl_sample_ind)
    df_crf_nona_celltype_info.loc[:, list_celltype + ["CellPropTotal"]] /= 100
    df_crf_nona_celltype_info.loc[:, list_celltype + ["CellPropTotal"]] = df_crf_nona_celltype_info.loc[:, list_celltype + ["CellPropTotal"]].applymap(lambda x: round(x, 4))

    return df_crf_nona_celltype_info

def add_deconv_prop_df(df_crf, df_deconv, colsample="Sample_ID", coltype="CellType"):
    list_samples = list(map(lambda x: x.replace("_sorted", ""), list(df_deconv.columns)))
    df_deconv.columns = list_samples
    df_deconv_blood = df_deconv[df_deconv[coltype].str.contains("Blood")].set_index(coltype).T
    df_deconv_blood = df_deconv_blood.applymap(lambda x: round(float(x), 5))
    df_deconv_blood = df_deconv_blood.reset_index(drop=False).rename(columns={"index": colsample})
    df_deconv_blood = df_deconv_blood.rename_axis(None, axis=1)
    df_deconv_info_added = df_crf.merge(df_deconv_blood, how="inner", on=colsample)
    
    return df_deconv_info_added

def convert_smoking_info_binary(smoking):
    if smoking == 1:
        binary = 1
    else:
        binary = 0
    
    return binary

def add_smoking_info_df(df_crf, col_smoke):
    df_crf_smoke_info_added = df_crf.copy()
    df_crf_smoke_info_added.loc[df_crf_smoke_info_added[col_smoke].isna(), col_smoke] = np.nan
    df_crf_smoke_info_added[col_smoke] = df_crf_smoke_info_added[col_smoke].apply(convert_smoking_info_binary)
    
    return df_crf_smoke_info_added

# def add_mediciation_info_df(df_crf, col_med="MedicationUse"):
#     df_crf[col_med] = df_crf[col_med].apply(lambda x: 1 if x != "무" else 0)
#     df_crf_med_info_added = df_crf
    
#     return df_crf_med_info_added

def add_mediciation_info_df(df_crf, col_steroid="SteroidUse", col_antiviral = "AntiviralUse"):
    df_crf[col_antiviral] = df_crf[col_antiviral].apply(lambda x: 1 if x != "NA" else 0)
    df_crf_med_info_added = df_crf
    
    return df_crf_med_info_added

def get_dict_covariates(path_comorbidity):
    with open(path_comorbidity, mode="rb") as fr:
        dict_covariates = dict(json.load(fr))
    
    return dict_covariates

def add_comorbidity_info_df(df_crf, dict_covariates, col_comorbid="Comorbidity"):
    string_list = list(chain(*list(dict_covariates.values())))
    pattern = re.compile("|".join(map(re.escape, string_list)), re.IGNORECASE)
    df_crf[col_comorbid] = df_crf[col_comorbid].astype(str).apply(lambda x: unicodedata.normalize("NFC", x) if isinstance(x, str) else x)
    df_crf[col_comorbid] = df_crf[col_comorbid].apply(lambda x: str(x).replace(x, "diabetes_mellitus") if "DM" in str(x) else x)
    df_crf[col_comorbid] = df_crf[col_comorbid].str.replace(r"\s+", "", regex=True)
    df_crf["Matched_Comorbidities"] = df_crf[col_comorbid].apply(lambda x: set(pattern.findall(x)) if isinstance(x, str) else None)

    df_crf_copy = df_crf.copy()
    for covariate_name, list_vals in dict_covariates.items():
        df_crf_copy[covariate_name] = df_crf_copy["Matched_Comorbidities"].apply(lambda x: int(bool(set(list_vals) & x)))

    df_crf_inf_v1_covariate_added = df_crf_copy.drop(columns=[col_comorbid, "Matched_Comorbidities"])
    
    return df_crf_inf_v1_covariate_added

def two_proportion_z_test(x1, n1, x2, n2, alternative='smaller'):
    count = np.array([x1, x2])  # Number of successes
    nobs = np.array([n1, n2])   # Sample sizes

    _, p_value = sm.stats.proportions_ztest(count, nobs, alternative=alternative)
    
    return p_value

def get_count_and_nobs_group_comparison(df_crf_inf_v1_covariate_added, covariates, col_compare="Severity_group", col_group1="Mild", col_group2="Severe"):
    list_cov_prop = list()
    list_mild_count = list()
    list_mild_nobs = list()
    list_sev_count = list()
    list_sev_nobs = list()

    df_crf_inf_v1_covariate_added_copy = df_crf_inf_v1_covariate_added.copy()
    
    for covariate in covariates:
        dict_cnt_groupby_sev_covariate = df_crf_inf_v1_covariate_added_copy.groupby(col_compare)[covariate].apply(Counter).to_dict()
        if float(list(dict_cnt_groupby_sev_covariate.keys())[0][1]) in [0, 1]:
            dict_cnt_groupby_sev_covariate_modif= {k: int(v) if str(v) != "nan" else 0 for k, v in dict_cnt_groupby_sev_covariate.items()}
            dict_cnt_groupby_sev_covariate_modif = {"_".join([str(k[0]), str(int(k[1]))]): int(v) for k,v in dict_cnt_groupby_sev_covariate_modif.items() if str(k[1]) != "nan"}

            mild_no = dict_cnt_groupby_sev_covariate_modif.get(f"{col_group1}_0", 0)
            mild_yes = dict_cnt_groupby_sev_covariate_modif.get(f"{col_group1}_1", 0)
            sev_no = dict_cnt_groupby_sev_covariate_modif.get(f"{col_group2}_0", 0)
            sev_yes = dict_cnt_groupby_sev_covariate_modif.get(f"{col_group2}_1", 0)
            
            counts_mild = int(mild_yes)
            nobs_mild = (int(mild_yes) + int(mild_no))
            counts_sev = int(sev_yes)
            nobs_sev = (int(sev_yes) + int(sev_no))
            
            list_cov_prop.append(covariate)
            list_mild_count.append(counts_mild)
            list_mild_nobs.append(nobs_mild)
            list_sev_count.append(counts_sev)
            list_sev_nobs.append(nobs_sev)
        else:
            print(f"Skipping {covariate}: value not binary!")
            continue
    
    return list_cov_prop, list_mild_count, list_mild_nobs, list_sev_count, list_sev_nobs

def plot_proportion_comparison(categories, outdir, successes_group1, total_group1, successes_group2, total_group2):
    proportions_group1 = np.array(successes_group1) / np.array(total_group1)
    proportions_group2 = np.array(successes_group2) / np.array(total_group2)

    p_values = [two_proportion_z_test(x1, n1, x2, n2) for x1, n1, x2, n2 in zip(successes_group1, total_group1, successes_group2, total_group2)]

    sorted_indices = np.argsort(np.maximum(proportions_group1, proportions_group2))
    categories = np.array(categories)[sorted_indices]
    proportions_group1 = proportions_group1[sorted_indices]
    proportions_group2 = proportions_group2[sorted_indices]
    p_values = np.array(p_values)[sorted_indices]

    fig, ax = plt.subplots(figsize=(8, len(categories) * 0.6))
    
    y_pos = np.arange(len(categories))
    width = 0.4

    bar1 = ax.barh(y_pos + width/2, proportions_group2 * 100 + 0.5, height=width, label="Severe", color="tomato", alpha=1, zorder=3)
    bar2 = ax.barh(y_pos - width/2, proportions_group1 * 100 + 0.5, height=width, label="Mild", color="forestgreen", alpha=1, zorder=3)
    
    for i, p_val in enumerate(p_values):
        max_x = max(proportions_group1[i], proportions_group2[i]) * 100 + 2
        ax.text(max_x - 0.5, y_pos[i], f"$P$={p_val:.3f}", va='center', fontsize=10)
        ax.vlines(x=max_x - 0.7, ymin=(y_pos[i] - 0.2), ymax=(y_pos[i] + 0.2), color="k", linestyle="-", alpha=0.6)
        ax.hlines(y=(y_pos[i] - 0.2), xmin=max_x - 1.2, xmax=max_x - 0.7, color="k", linestyle="-", alpha=0.6)
        ax.hlines(y=(y_pos[i] + 0.2), xmin=max_x - 1.2, xmax=max_x - 0.7, color="k", linestyle="-", alpha=0.6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(categories)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Proportion (%)")
    ax.set_title("Proportions of clinical variables grouped by severity scale")
    ax.legend(loc="lower right", frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_zorder(6)
    ax.spines['bottom'].set_zorder(6)
        
    plt.grid(axis="x", zorder=-1)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "plot_proportion_comparison_clinical_values_severity.png"), dpi=300)
    plt.show()
    plt.close()
        
def save_stats_proportion_covariates(df_crf_inf_v1_covariates, covariates):
    list_cov_prop, list_mild_count, list_mild_nobs, list_sev_count, list_sev_nobs = get_count_and_nobs_group_comparison(df_crf_inf_v1_covariates, covariates, col_compare="Severity_group", col_group1="Mild", col_group2="Severe")
    proportions_mild = np.array(list_mild_count) / np.array(list_mild_nobs)
    proportions_sev = np.array(list_sev_count) / np.array(list_sev_nobs)
    p_values = [two_proportion_z_test(x1, n1, x2, n2) for x1, n1, x2, n2 in zip(list_mild_count, list_mild_nobs, list_sev_count, list_sev_nobs)]
    dict_prop_info = dict()
    for cov_prop, mild_count, mild_obs, mild_prop, sev_count, sev_obs, sev_prop, pval in zip(list_cov_prop, list_mild_count, list_mild_nobs, proportions_mild, list_sev_count,list_sev_nobs, proportions_sev, p_values):
        dict_prop_info[cov_prop] = [mild_count, mild_obs, round(float(mild_prop), 3), sev_count, sev_obs, round(float(sev_prop), 3), round(float(pval), 3)]
    
    df_prop_info = pd.DataFrame.from_dict(dict_prop_info, orient="index")
    df_prop_info.columns = ["Mild_Count", "Mild_Nobs", "Mild_Prop", "Sev_Count", "Sev_Nobs", "Sev_Prop", "Pvalue(one_tailed_two-sample_z_prop_test_sev_vs_mild)"]
    df_prop_info.to_csv(os.path.join(dir_metadata_incl_clin, "prop_info_covariates.txt"), sep="\t", index=True)
    
    return df_prop_info
    
def plot_proportion_covariates(df_crf_inf_v1_covariate_added, dict_covariates, outdir):
    covariates = list(dict_covariates.keys())
    list_cov_prop, list_mild_count, list_mild_nobs, list_sev_count, list_sev_nobs = get_count_and_nobs_group_comparison(df_crf_inf_v1_covariate_added, covariates, col_compare="Severity_group", col_group1="Mild", col_group2="Severe")
    plot_proportion_comparison(list_cov_prop, outdir, list_mild_count, list_mild_nobs, list_sev_count, list_sev_nobs)

# %%
if __name__ == "__main__":
    flag_celltype=False
    flag_deconv=False
    flag_smoking=True
    flag_medication=False
    flag_comorbidity=True
    WORKDIR = str(Path(os.path.abspath(__file__)).parents[3])
    path_comorbidity = f"{WORKDIR}/Resources/Scripts/Final/list_comorbidity.json"
    path_rename = f"{WORKDIR}/Resources/Scripts/Final/list_variable_name_change.json"
    path_drop_samples = f"{WORKDIR}/Resources/Scripts/Final/list_remove_samples.txt"
    # path_ehr = f"{WORKDIR}/Resources/Data/EHR/infectomics_CRF_20230410_edit.xlsx"
    path_ehr = f"{WORKDIR}/Resources/Data/EHR/Infectomics_CRF_UNIST_20250314.xlsx"
    path_deconv = f"{WORKDIR}/Results/13_uxm/COVID-19_V1_Conval_hg38_sortedByseverity.csv"
    dir_metadata = f"{WORKDIR}/Resources/Data/Methylation/MetaTable_LOO"
    pattern_metadata = "metatable_combined_all_firstVisit_WO_*.tsv"
    # pattern_metadata = "Methylseq_master_combined.tsv"
    # dir_metadata_incl_clin = f"{WORKDIR}/Resources/Data/Methylation/MetaTable_Incl_EHR_CellType_SmokingStatus_CCIScore"
    dir_metadata_incl_clin = f"{WORKDIR}/Resources/Data/Methylation/MetaTable_Incl_EHR_SmokingStatus_CCIScore"
    os.makedirs(dir_metadata_incl_clin, exist_ok=True)

    main(path_ehr, path_deconv, path_drop_samples, path_rename, path_comorbidity, dir_metadata, pattern_metadata, dir_metadata_incl_clin, flag_celltype, flag_deconv, flag_smoking, flag_medication, flag_comorbidity)

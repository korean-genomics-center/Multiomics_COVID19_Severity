# %%
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import clinical

warnings.filterwarnings("ignore")

# %%
def get_drop_samples(path_drop_samples):
    with open(path_drop_samples, mode="r") as frmv:
        list_drop_samples = [x.rstrip() for x in frmv.readlines()]
    
    return list_drop_samples

def get_clin_name(df_compare):
    clinical_name = list(df_compare.columns)[-1]
    
    return clinical_name

def calc_stats(df_compare, outdir):
    df_compare_copy = df_compare.copy()
    stat_clin = clinical.Stats(df_compare_copy, outdir)
    try:
        dict_stat_res = stat_clin.calculate_diff_stats()
        
        return dict_stat_res
    
    except Exception as e:
        print(e)

# %%
path_crf = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/infectomics_CRF_20230410_edit.xlsx"
path_rename = "/BiO/Access/kyungwhan1998/Infectomics/Resources/Scripts/Final/list_variable_name_change.json"
path_drop_samples = "/BiO/Access/kyungwhan1998/Infectomics/Resources/Scripts/Final/list_remove_samples.txt"
dict_palette = {"severe":"firebrick", 
                "mild":"forestgreen"}
outdir = "/BiO/Access/kyungwhan1998/Infectomics/Results/InfectomicsPaper1/20240906"
target = "Severity"
os.makedirs(outdir, exist_ok=True)

# %%
list_drop_samples = get_drop_samples(path_drop_samples)
parse_clin = clinical.Parse(path_crf, list_drop_samples, path_rename)
df_crf_raw = parse_clin.read_excel(engine="openpyxl", 
                                   sheet_name="21-22등록 대상자_modify", 
                                   skiprows=1)
df_crf = parse_clin.give_parsed_dataframe()
df_crf_inf = df_crf[df_crf["SampleID"].str.startswith("C19-C")]
df_crf_inf_v1 = df_crf_inf[df_crf_inf["SampleID"].str.endswith("V1")]
# lookup_visit = parse_clin.lookup_dict_sample_visit(df_crf_inf_v1)
list_df_compare = parse_clin.get_clinical_assoc_severity(df_crf_inf_v1, target)

list_clin_name = Parallel(n_jobs=len(list_df_compare))(delayed(get_clin_name)(df_compare) for df_compare in list_df_compare)
list_res_stat = Parallel(n_jobs=len(list_df_compare))(delayed(calc_stats)(df_compare, outdir) for df_compare in list_df_compare)
dict_clin_stat = {k: v for k, v in dict(zip(list_clin_name, list_res_stat)).items() if v is not None}
df_clin_stat = pd.DataFrame.from_dict(dict_clin_stat).T
# df_clin_2_sample_prop = df_clin_stat[df_clin_stat["test"]=="2_samp_prop"]
df_clin_mwu = df_clin_stat[df_clin_stat["test"]=="mwu"]
# dict_clin_stat_mwu = {k: v for k, v in dict_clin_stat.items() if v["test"]=="mwu"}
# list(df_clin_mwu[df_clin_mwu["pval"] < 0.05].index)
list_clinicals = ["SmokingYears", 
                    'CT(PCR)_E', 
                    'CT(PCR)_R', 
                    'CT(PCR)_N',
                    'Neutrophil',
                    'Monocytes',
                    'Glucose',
                    'LDH',
                    'CPK',
                    'FDP',
                    'Protein_S_Activity',
                    'CRP']

# %%
num_plots = len(list_clinicals)
cols = int(np.sqrt(num_plots))
rows = int(np.ceil(num_plots / cols))  # Ensure enough rows to accommodate all subplots
fig, axes = plt.subplots(rows, cols, figsize=(3 * cols, 4 * rows))  # Adjust width and height accordingly
axes = axes.flatten()  # Flatten in case of multiple rows and columns
i = 0  # Move counter outside the loop

for var_clin, df_compare in zip(list_clin_name, list_df_compare):
    if var_clin in list_clinicals:
        stat_clin = clinical.Stats(df_compare, outdir)
        dict_stat_res = stat_clin.calculate_diff_stats()
        plot_clin = clinical.Plot(dict_stat_res, dict_palette, outdir, axes[i])
        g, pvalues = plot_clin.draw_boxplot(df_compare, var_clin)
        g.yaxis.grid(True, zorder=-1)
        if float(pvalues[0]) < 0.05:
            axes[i].set_ylabel(ylabel=var_clin, color="firebrick", weight="bold", fontsize=20)
        
        i += 1

plt.subplots_adjust(hspace=0.3, wspace=0.8)
plt.savefig(f"{outdir}/SupplementaryFigure1.png", dpi=600)
plt.show()
plt.close()
# %%

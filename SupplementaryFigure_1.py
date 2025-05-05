# %%
import os
import warnings
from pathlib import Path
import clinical
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

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
WORKDIR = str(Path(__file__).parents[3])
path_crf = f"{WORKDIR}/Resources/Data/EHR/infectomics_CRF_20230410_edit.xlsx"
path_rename = f"{WORKDIR}/Resources/Scripts/Final/list_variable_name_change.json"
path_drop_samples = f"{WORKDIR}/Resources/Scripts/Final/list_remove_samples.txt"
dict_palette = {"Severe":"firebrick", 
                "Mild":"darkgreen"}
outdir = f"{WORKDIR}/Results/Paper"
colsev = "Severity_group"
list_clinicals = ['CT(PCR)_E', 
                'CT(PCR)_R', 
                'CT(PCR)_N',
                'RBC',
                'WBC',
                'Neutrophil',
                'Lymphocytes',
                'Monocytes',
                'Basophils',
                'Eosinophils',
                'Glucose',
                'CRP']
os.makedirs(outdir, exist_ok=True)

# %%
list_drop_samples = get_drop_samples(path_drop_samples)
parse_clin = clinical.Parse(path_crf, list_drop_samples, path_rename)
df_crf_raw = parse_clin.read_excel(engine="openpyxl", 
                                   sheet_name="21-22등록 대상자_modify", 
                                   skiprows=1)
df_crf = parse_clin.give_parsed_dataframe()
df_crf_inf = df_crf[df_crf["Sample_ID"].str.startswith("C19-C")]
df_crf_inf_v1 = df_crf_inf[df_crf_inf["Sample_ID"].str.endswith("V1")]
# lookup_visit = parse_clin.lookup_dict_sample_visit(df_crf_inf_v1)
list_df_compare = parse_clin.get_clinical_assoc_severity(df_crf_inf_v1, colsev)

# %%
list_clin_name = Parallel(n_jobs=len(list_df_compare))(delayed(get_clin_name)(df_compare) for df_compare in list_df_compare)
list_res_stat = Parallel(n_jobs=len(list_df_compare))(delayed(calc_stats)(df_compare, outdir) for df_compare in list_df_compare)
dict_clin_stat = {k: v for k, v in dict(zip(list_clin_name, list_res_stat)).items() if v is not None}
df_clin_stat = pd.DataFrame.from_dict(dict_clin_stat).T
df_clin_mwu = df_clin_stat[df_clin_stat["test"]=="ranksum"]

# %%
plt.rcParams["font.size"] = 10
num_plots = len(list_clinicals)
cols = 4
rows = int(np.ceil(num_plots / cols))  
fig, axes = plt.subplots(rows, cols, figsize=(3 * cols, 6 * rows))
axes = axes.flatten()
i = 0 
for var_clin, df_compare in zip(list_clin_name, list_df_compare):
    if var_clin in list_clinicals:
        stat_clin = clinical.Stats(df_compare, outdir)
        dict_stat_res = stat_clin.calculate_diff_stats()
        plot_clin = clinical.Plot(dict_stat_res, dict_palette, outdir, axes[i])
        g, pvalues = plot_clin.draw_boxplot(df_compare, var_clin)
        g.grid(axis="y", zorder=1)
        g.set_axisbelow(True)
        i += 1

plt.subplots_adjust(wspace=0.8, hspace=0.5)
plt.savefig(f"{outdir}/SupplementaryFigure1.png", bbox_inches="tight", dpi=300)
plt.savefig(f"{outdir}/SupplementaryFigure1.pdf", bbox_inches="tight", dpi=300)
plt.show()
plt.close()

# %%

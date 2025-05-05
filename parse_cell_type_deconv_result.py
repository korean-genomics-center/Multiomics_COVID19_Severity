# %%
import os

import numpy as np
import pandas as pd

# %%
path_deconv = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/cell_type_deconv/UXM/Results/run_250220/COVID-19_V1_hg38_sortedByseverity.csv"
df_deconv = pd.read_csv(path_deconv, sep=",")
list_cells = df_deconv["CellType"].to_list()
list_samples = list(map(lambda x: x.replace("_sorted", ""), list(df_deconv.columns)))
df_deconv.columns = list_samples
df_deconv_blood = df_deconv[df_deconv["CellType"].str.contains("Blood")].set_index("CellType").T
df_deconv_others = df_deconv[~df_deconv["CellType"].str.contains("Blood")].set_index("CellType").T
df_deconv_blood["Blood_Total"] = df_deconv_blood.sum(axis=1)
df_deconv_blood["Lymphocytes"] = df_deconv_blood[["Blood-B", "Blood-T", "Blood-NK"]].agg(sum, axis=1)
df_deconv_blood["Granulocytes"] = df_deconv_blood["Blood-Granul"]
df_deconv_blood["Monocytes"] = df_deconv_blood["Blood-Mono+Macro"]
df_deconv_others["Others_Total"] = df_deconv_others.agg(sum, axis=1)
df_deconv_blood["Others"] = df_deconv_others["Others_Total"]
list_col_fin = list(filter(lambda x: "Blood" not in x, list(df_deconv_blood.columns)))
df_deconv_blood = df_deconv_blood[list_col_fin]
df_deconv_blood.columns = list(map(lambda x: x + "_deconv", list_col_fin))
df_deconv_blood = df_deconv_blood.reset_index(drop=False).rename(columns={"index": "Sample_ID"})

# %%
path_ehr = "/BiO/Access/kyungwhan1998/Infectomics/Resources/Data/Methylation/MetaTable_Incl_EHR/Methylseq_master_combined_celltype_smoking_medication_comorbidity_info_added.tsv"
df_ehr = pd.read_csv(path_ehr, sep="\t")
list_celltype = ["Neutrophil", "Lymphocytes", "Monocytes", "Eosinophils", "Basophils"]
df_ehr_celltype = df_ehr[["Sample_ID"] + list_celltype]
df_ehr_celltype["Granulocytes"] = df_ehr_celltype[["Neutrophil", "Eosinophils", "Basophils"]].agg(sum, axis=1)
df_ehr_celltype = df_ehr_celltype.drop(columns = ["Neutrophil", "Eosinophils", "Basophils"])
df_ehr_celltype = df_ehr_celltype.set_index("Sample_ID")
list_col_ori = list(df_ehr_celltype.columns)
df_ehr_celltype.columns = list(map(lambda x: x + "_count", list_col_ori))
df_ehr_celltype = df_ehr_celltype.reset_index(drop=False)

# %%
df_merged = df_deconv_blood.merge(df_ehr_celltype, how="inner", on="Sample_ID")

# %%
import math
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr

list_comb_name = [(x, y) for x, y in combinations(df_merged.columns, 2) if x.split("_")[0] == y.split("_")[0] and x.split("_")[1] != y.split("_")[1]]
num_comb = len(list_comb_name)
list_colors = ["green", "orange", "blue"]
fig, axes = plt.subplots(ncols=4, nrows=1, figsize=(11, 2))
for ind, (name_1, name_2) in enumerate(list_comb_name):
    x_vals = df_merged[name_1].dropna()
    y_vals = df_merged[name_2].dropna()
    sns.regplot(x=x_vals, y=y_vals, ax=axes[ind], color=list_colors[ind], scatter_kws={"s":15, "alpha":0.8})
    sns.regplot(x=x_vals, y=y_vals, ax=axes[-1], color=list_colors[ind], scatter_kws={"s":15, "alpha":0.2})
    axes[-1].scatter([],[],color = list_colors[ind], label=name_1.split("_")[0])
    corr, pval = pearsonr(x=x_vals, y=y_vals)
    
    x_min, x_max = x_vals.min(), x_vals.max()
    y_min, y_max = y_vals.min(), y_vals.max()

    # x_pad = (x_max - x_min) * 0.2
    # y_pad = (y_max - y_min) * 0.2
    # axes[ind].set_xlim(x_min - x_pad, x_max + x_pad)
    # axes[ind].set_ylim(y_min - y_pad, y_max + y_pad)

    # tick_min = min(x_min, y_min) - min(x_pad, y_pad)
    # tick_max = max(x_max, y_max) + max(x_pad, y_pad)
    # tick_spacing = (tick_max - tick_min) / 5 
    # ticks = np.arange(tick_min, tick_max, tick_spacing)

    # axes[ind].set_xticks(ticks)
    # axes[ind].set_yticks(ticks)
    
    tick_min = min(x_min, y_min)
    tick_max = max(x_max, y_max)
    padding = 0.005
    
    tick_min_padded = tick_min - (tick_max - tick_min) * padding
    tick_max_padded = tick_max + (tick_max - tick_min) * padding
    
    # single_tick = 0.1 if (tick_max_padded - tick_min_padded) > 0.3 else 0.05
    if (tick_max_padded - tick_min_padded) < 0.5:
        adjust_param = 0.05
    else:
        adjust_param = 0.1
    single_tick = round(math.ceil(((tick_max_padded - tick_min_padded) / 4)/adjust_param) * adjust_param, 3)
    
    tick_min_padded_floor = round(math.floor(tick_min_padded / single_tick) * single_tick, 3)
    tick_max_padded_ceil = round(math.ceil(tick_max_padded / single_tick) * single_tick, 3)
    
    n_ticks_show = round((tick_max_padded_ceil - tick_min_padded_floor) / single_tick)
    
    print(tick_min_padded_floor, tick_max_padded_ceil, single_tick, n_ticks_show)
    
    axes[ind].set_xticks(np.linspace(tick_min_padded_floor, tick_max_padded_ceil, n_ticks_show+1))
    axes[ind].set_yticks(np.linspace(tick_min_padded_floor, tick_max_padded_ceil, n_ticks_show+1))
    
    textpadding = 0.03
    axes[ind].text(x=textpadding, y=1-textpadding, s=f"corr:{corr:.2f}\npval:{pval:.1e}", va="top", transform = axes[ind].transAxes)
    
    axes[ind].set_xlabel("Deconvoluted")
    axes[ind].set_ylabel("Measured")
    axes[ind].set_title(name_1.split("_")[0], weight="bold", fontsize=11)

axes[-1].set_xlim(0, 1)
axes[-1].set_ylim(0, 1)
axes[-1].set_xticks(np.arange(0, 1.1, 0.2))
axes[-1].set_yticks(np.arange(0, 1.1, 0.2))
axes[-1].set_title("Combined", weight="bold", fontsize=11)
axes[-1].set_xlabel("Deconvoluted")
axes[-1].set_ylabel("Measured")

fig.subplots_adjust(wspace=0.4)
axes[-1].legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
plt.show()
plt.close()
# %%

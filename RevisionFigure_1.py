# %%
import os
import pandas as pd
from pathlib import Path

# %%
WORKDIR = str(Path(__file__).parents[3])
path_uxm = f"{WORKDIR}/Results/13_uxm/COVID-19_V1_hg38_sortedByseverity.csv"
path_meta = f"{WORKDIR}/Resources/Data/RNA/COVID19_master_table_added_CRF_20231007.txt"
df_uxm = pd.read_csv(path_uxm, sep=",")
df_uxm_set_index = df_uxm.set_index("CellType")
list_celltype_blood = list(filter(lambda x: x.startswith("Blood"), list(df_uxm_set_index.index)))
list_celltype_other = list(set(list(df_uxm_set_index.index)).difference(set(list_celltype_blood)))
df_uxm_celltype_blood = df_uxm_set_index.loc[list_celltype_blood]
df_uxm_celltype_other = df_uxm_set_index.loc[list_celltype_other]
dict_uxm_celltype_other_sum = dict(df_uxm_celltype_other.sum(axis=0))
df_uxm_celltype_blood.loc["Others"] = dict_uxm_celltype_other_sum
df_uxm_celltype_blood_transposed = df_uxm_celltype_blood.T
df_uxm_celltype_blood_transposed.index = list(map(lambda x: x.replace("_sorted", ""), list(df_uxm_celltype_blood_transposed.index)))

df_meta = pd.read_csv(path_meta, sep="\t")
df_meta_set_index = df_meta.set_index("Project_ID_Alias")
df_meta_covid_v1 = df_meta_set_index.loc[list(df_uxm_celltype_blood_transposed.index)]
dict_sev = dict(zip(df_meta_covid_v1.index, df_meta_covid_v1["Severity_Binary"]))

dict_uxm_celltype_blood_sev_added = dict()
for ind, row in df_uxm_celltype_blood_transposed.iterrows():
    dict_row = dict(row)
    if ind in dict_sev:
        sev = dict_sev.get(ind, None)
        dict_row["Severity_group"] = sev 
    
    dict_uxm_celltype_blood_sev_added[ind] = (dict_row)

df_uxm_celltype_blood_sev_added = pd.DataFrame.from_dict(dict_uxm_celltype_blood_sev_added, orient="index")

# %%
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from statannotations.Annotator import Annotator

nrows = 1
ncols = -(-len(list_celltype_blood) // nrows)
list_celltype_blood = ["Blood-Granul", "Blood-T", "Blood-B", "Blood-NK", "Blood-Mono+Macro", "Others"]
dict_title = {"Blood-Granul": "Granul.", "Blood-T": "T", "Blood-B": "B", "Blood-NK": "NK", "Blood-Mono+Macro": "Mono.&\nMacro."}
fig, axes = plt.subplots(nrows, ncols, figsize=(6, 5))
axes = axes.flatten()

pairs = [("Mild", "Severe")]
palette = { "Mild": "darkgreen", "Severe": "firebrick"}
order = ["Mild", "Severe"]

for i, (celltype, ax) in enumerate(zip(list_celltype_blood, axes)):
    ax.set_axisbelow(True)
    sns.boxplot(data=df_uxm_celltype_blood_sev_added, x="Severity_group", y=celltype, width=0.5, showfliers=False, ax=ax, palette=palette, order=order, zorder=1)
    sns.stripplot(data=df_uxm_celltype_blood_sev_added, x="Severity_group", y=celltype, ax=ax, marker='o', size=4, facecolor=None, edgecolor='black', linewidth=0.5, order=order, jitter=False, alpha=0.6, zorder=2, palette=palette, label=None)
    annotator = Annotator(ax, pairs, data=df_uxm_celltype_blood_sev_added, x="Severity_group", y=celltype)
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=0)
    annotator.apply_and_annotate()
    ax.set_title(dict_title.get(celltype, celltype))
    ax.set_xticklabels(["Mild-Moderate", "Severe-Critical"], rotation=45, ha="right", rotation_mode="anchor")
    ax.set_xlabel("")
    ax.set_ylim(-0.09, 1.09)
    ax.set_yticks(np.arange(0, 1.1, 0.1))
    if i == 0:
        ax.set_ylabel("Cell-type Proportion", fontsize=13)
        ax.set_yticklabels(range(0, 110, 10))
    else:
        ax.set_ylabel("")
        ax.set_yticklabels([])
    
    ax.grid(axis="y", zorder=0)

# Remove unused axes
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.subplots_adjust(wspace=0.3)
plt.show()

# %%

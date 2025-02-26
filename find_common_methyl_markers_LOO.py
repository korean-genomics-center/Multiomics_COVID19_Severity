# %%
import glob
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

import math_trick
import matplotlib.pyplot as plt
import pandas as pd
from filter_MethylKitTable_by_position import filter_by_pos


def get_dict_dmp_loo(dir_dmp_loo, pattern):
    dict_dmp_loo = defaultdict(list)
    
    list_dmp_loo = glob.glob(f"{dir_dmp_loo}/{pattern}", recursive=True)
    for dmp_loo in list_dmp_loo:
        df_dmp_loo = pd.read_csv(dmp_loo, sep="\t")
        df_dmp_loo["chrpos"] = df_dmp_loo["chr"] + ":" + df_dmp_loo["start"].astype(str)
        list_dmp_loo = df_dmp_loo["chrpos"].to_list()
        list_dmp_methdiff = df_dmp_loo["meth.diff"].to_list()
        for dmp, methdiff in zip(list_dmp_loo, list_dmp_methdiff):
            dict_dmp_loo[dmp].append(methdiff)
    
    return dict_dmp_loo

def plot_common_methyl_markers_loo(dict_dmp_loo, outdir, covariate):
    num_common_loo = list(map(lambda x: len(x), list(dict_dmp_loo.values())))        
    dict_cnt = dict(Counter(num_common_loo))
    plt.bar(dict_cnt.keys(), dict_cnt.values(), color="grey")
    for x, y in dict_cnt.items():
        plt.text(x, y + 10, str(y), ha='center', va='bottom', fontsize=12)
    plt.xticks(range(1, max(dict_cnt.keys())+1), fontsize=12)
    plt.yticks(range(0, math_trick.round_up(max(dict_cnt.values()))+1, 1000), fontsize=12)
    plt.xlabel("Number of overlaps between LOO rounds", fontsize=16)
    plt.ylabel("Count", fontsize=16)
    plt.title(covariate, fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"plot_common_methyl_markers_loo_{covariate}.png"), dpi=300)
    plt.show()
    plt.close()

def filter_methylkittable_dmp_loo_by_pos(path_methykittable, dict_dmp_loo, covariate, fold, outdir):
    dict_common_dmp_loo = {k: v for k, v in dict_dmp_loo.items() if len(v) >= fold}
    list_common_dmp_loo = list(dict_common_dmp_loo.keys())

    outfilename = f"MethylCpGTable.Control.Mild.Case.Severe.filtered.Covariate_{covariate}_Yes_LOO_common_across_{fold}folds_or_more.tsv"
    path_methykittable_filtered = os.path.join(outdir, outfilename)

    df_methylkittable = pd.read_csv(path_methykittable, sep="\t")
    df_methylkittable_filtered = filter_by_pos(df_methylkittable, list_common_dmp_loo, path_methykittable_filtered, save=True)
    
    return df_methylkittable_filtered

def main(dir_dmp_loo, path_methykittable, covariate, fold, outdir, pattern = "*.fdr_05.sorted.tsv"):
    print(f"running: {covariate}...")
    dict_dmp_loo = get_dict_dmp_loo(dir_dmp_loo.format(covariate), pattern)
    plot_common_methyl_markers_loo(dict_dmp_loo, outdir, covariate)
    filter_methylkittable_dmp_loo_by_pos(path_methykittable, dict_dmp_loo, covariate, fold, outdir)
    
    print("Done!")
    
    return None

# %%
from joblib import Parallel, delayed

WORKDIR = str(Path(os.path.abspath(__file__)).parents[3])
dir_dmp = f"{WORKDIR}/Results/11_dmp"
dir_dmp_loo = os.path.join(dir_dmp, "Yes_LOO/CellDeconvAdjusted/Covariate_{0}/DMPExtract")
path_methykittable = f"{WORKDIR}/Resources/Data/Methylation/MethylKitTable/MethylCpGTable.Control.Mild.Case.Severe.tsv"
outdir = f"{WORKDIR}/Resources/Data/Methylation/MethylKitTable_Incl_EHR_Deconvoluton"
os.makedirs(outdir, exist_ok=True)
list_covariates = ["Sex_Age", "Sex_Age_Smoking_Medication_Comorbidity", "Sex_Age_CellDeconv"]
# list_covariates = ["Sex_Age_CellDeconv_Smoking_Medication_Comorbidity"]
fold = 9

# %%
n_jobs = len(list_covariates)
with Parallel(n_jobs=n_jobs) as parallel:
    parallel(delayed(main)(dir_dmp_loo, path_methykittable, covariate, fold, outdir) for covariate in list_covariates)
# %%

# %%
import os
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed


# %%
def main(path_methylkittable, path_dmp_deg_corr, path_sig_corr_only, covariate, fold):
    df_methylkittable = pd.read_csv(path_methylkittable.format(covariate, fold), sep="\t")
    df_methylkittable["chrpos"] = df_methylkittable["chr"] + "_" + df_methylkittable["start"].astype(str) + "_" + df_methylkittable["end"].astype(str)
    list_common_dmp_loo = df_methylkittable["chrpos"].to_list()
    
    with open(path_sig_corr_only.format(covariate, fold), mode="w") as fw:
        with open(path_dmp_deg_corr, mode="r") as fr:
            header = fr.readline().rstrip("\n").split("\t")
            new_header = ["Methyl", "RNA" , "corr_rho", "corr_fdr", "methyl_chr", "rna_chr"]
            fw.write("\t".join(new_header) + "\n")
            ind_methyl = header.index("Methyl")
            ind_rna = header.index("RNA")
            ind_fdr = header.index("corr_fdr")
            ind_rho = header.index("corr_rho")
            ind_chr = header.index("RNA_chr")
            for line in fr:
                if line is None:
                    break
                record = line.rstrip("\n").split("\t")
                methyl_marker = record[ind_methyl].strip()
                gene_marker = record[ind_rna].strip()
                methyl_fdr = record[ind_fdr].strip()
                methyl_rho = record[ind_rho].strip()
                methyl_chr = methyl_marker.split("_")[0]
                rna_chr = record[ind_chr].strip()
                try:
                    methyl_fdr = round(float(methyl_fdr), 5)
                    methyl_rho = round(float(methyl_rho), 5)
                    if (methyl_marker in set(list_common_dmp_loo)) and (methyl_fdr < 0.05) and (methyl_chr == rna_chr):
                        fw.write("\t".join([str(methyl_marker), str(gene_marker), str(methyl_rho), str(methyl_fdr)]) + "\n")
                        
                except ValueError:
                    print(f"Skipping invalid values: marker={methyl_marker}, fdr={methyl_fdr}, rho={methyl_rho}")

WORKDIR = str(Path(__file__).parents[3])
path_methylkittable = os.path.join(WORKDIR, "Resources/Data/Methylation/MethylKitTable_Incl_EHR_Deconvoluton/MethylCpGTable.Control.Mild.Case.Severe.filtered.Covariate_{0}_Yes_LOO_common_across_{1}folds_or_more.tsv")
path_dmp_deg_corr = os.path.join(WORKDIR, "Results/12_cis_eQTL/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.20240326.tsv")
path_sig_corr_only = os.path.join(WORKDIR, "Results/11_dmp/Yes_LOO/CellDeconvAdjusted/Covariate_{0}/Methylation_DMP_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common_across_{1}.cis.distance_cf.1000000.tsv")
list_covariates = ["Sex_Age", "Sex_Age_Smoking_Medication_Comorbidity", "Sex_Age_CellDeconv"]
# list_covariates = ["Sex_Age_CellDeconv_Smoking_Medication_Comorbidity"]
fold = 4

n_jobs = len(list_covariates)
with Parallel(n_jobs=n_jobs) as parallel:
    parallel(delayed(main)(path_methylkittable, path_dmp_deg_corr, path_sig_corr_only, covariate, fold) for covariate in list_covariates)
# %%

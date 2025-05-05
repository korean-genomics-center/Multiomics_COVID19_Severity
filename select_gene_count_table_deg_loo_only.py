# %%
import os
from pathlib import Path

# %%
WORKDIR = str(Path(os.path.abspath(__file__)).parents[3])
path_exp = f"{WORKDIR}/Results/5_deg/Visit1_Severe__Visit1_Mild_20240327.tsv.normcount"
path_list_gene = f"{WORKDIR}/Resources/Data/Methylation/DEG_LOO/list_genes_over_loo_overlap_threshold_7.txt"
path_exp_gene_filt = f"{WORKDIR}/Results/5_deg/Visit1_Severe__Visit1_Mild_20240327_genes_over_loo_overlap_threshold_7_filt.normcount"

# %%
with open(path_list_gene, mode='r') as frg:
    list_genes = list(map(lambda x: x.rstrip("\n"), frg.readlines()))
    
with open(path_exp, mode='r') as fr, open(path_exp_gene_filt, mode='w') as fw:
    sample_id = fr.readline().rstrip("\n").split("\t")
    header = sample_id
    fw.write("\t".join(header) + "\n")
    for line in fr:
        record = line.rstrip("\n").split("\t")
        gene_id = record[0]
        if gene_id in set(list_genes):
            fw.write("\t".join(record) + "\n")
# %%

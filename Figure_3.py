# %%
import glob
import gzip
import os
import pickle
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import gridspec

# %%
VISIT_CONV = {
    "Visit0" : "Healthy",
    "Visit1" : "Acute",
    "Visit3" : "Recovery",
    "Visit4" : "Recovery",
    "Visit5" : "Convalescent"
}
SEV_CONV = {
    "0": "Healthy",
    "1" : "Mild",
    "2" : "Mild",
    "3" : "Severe",
    "4" : "Severe"
}
PALETTE = {
    "Healthy controls": (0.6627450980392157, 0.6627450980392157, 0.6627450980392157, 1.0),
    "Convalescent" : (0.2549019607843137, 0.4117647058823529, 0.8823529411764706, 1.0),
    "Mild (Recovery)" : (0.13333333333333333, 0.5450980392156862, 0.13333333333333333, 1.0),
    "Mild (Acute)" : (0.13333333333333333, 0.5450980392156862, 0.13333333333333333, 1.0),
    "Severe (Recovery)" : (0.6980392156862745, 0.13333333333333333, 0.13333333333333333, 1.0),
    "Severe (Acute)" : (0.6980392156862745, 0.13333333333333333, 0.13333333333333333, 1.0)
}
MARKER = {
    "Healthy controls" : "D",
    "Convalescent" : "D",
    "Mild (Recovery)" : "^",
    "Mild (Acute)" : "o",
    "Severe (Recovery)" : "^",
    "Severe (Acute)" : "o"
}

# %%
def load_pickle(loadfilename):
    with gzip.open(loadfilename,'rb') as fr:
        data = pickle.load(fr)
    
    return data

def get_list_files_methylcpgmin(dirmethylcpgmin):
    list_files_methylcpgmin = glob.glob(f"{dirmethylcpgmin}/**/*pair_merged.methyl_cpg_min.tsv", recursive=True)

    return list_files_methylcpgmin

def get_list_methyl_sample(list_files_methylcpgmin, list_drop_samples):
    list_methyl_sample = list()
    for file_methylcpgmin in list_files_methylcpgmin:
        dirmethylcpgmin = os.path.basename(os.path.dirname(file_methylcpgmin))
        if dirmethylcpgmin == "HealthyControl":
            name_sample = os.path.basename(file_methylcpgmin).split(".")[0]
        else:
            name_sample = os.path.basename(dirmethylcpgmin)
        list_methyl_sample.append(name_sample)

    list_methyl_sample = list(filter(lambda x: "C19-C" in x, list_methyl_sample))
    list_methyl_sample = list(filter(lambda x: x not in list_drop_samples , list_methyl_sample))
    list_methyl_sample = list(filter(lambda x: "L1" not in x, list_methyl_sample))

    return list_methyl_sample

def get_dict_sample_severity(pathsevinfo, list_methyl_sample):
    dict_sample_severity = dict()
    with open (pathsevinfo, mode="r") as fr:
        list_header = fr.readline().rstrip("\n").split("\t")
        idx_id = list_header.index("Sample_ID")
        idx_sev_visit = list_header.index("Severity_visit")
        for line in fr:
            record = line.rstrip("\n").split("\t")
            sampleid = record[idx_id]
            sev_vis = record[idx_sev_visit]
            if sampleid in list_methyl_sample:
                sev, vis = sev_vis.split("_")
                vis_conv = VISIT_CONV.get(vis, vis)
                sev_vis = sev + "_" + vis_conv
                dict_sample_severity[sampleid] = sev_vis

    return dict_sample_severity  

def sort_dict_sample_severity_visit(dict_sample_severity):
    dict_sample_severity_sort = dict(sorted(dict_sample_severity.items(), key=lambda x: (x[1].split("_")[1], x[1].split("_")[0])))
    
    return dict_sample_severity_sort

def get_dmp_deg_stats(filetargetgene, colmethdiff="dmp_median_methdiff", colfc="deg_median_log2FC"):
    colfcabs = "abs" + colfc
    df_dmp_deg = pd.read_csv(filetargetgene, sep="\t")
    df_dmp_deg[colfcabs] = df_dmp_deg[colfc].apply(abs)
    df_dmp_deg_sorted = df_dmp_deg.sort_values(by=[colmethdiff, colfcabs], ascending=False)
    
    return df_dmp_deg_sorted

def get_dict_target_marker_gene(df_dmp_deg_sorted, colmeth="Methyl", colrna="RNA"):
    dict_marker_gene = defaultdict(set)
    for _, rows in df_dmp_deg_sorted.iterrows():
        dict_rows = dict(rows)
        gene_sym = dict_rows[colrna]
        dmp_markers = ":".join(dict_rows[colmeth].split("_")[:-1])
        dict_marker_gene[dmp_markers].add(gene_sym)
    
    return dict(dict_marker_gene)

def get_dict_target_marker_direction(df_dmp_deg_sorted, colmeth="Methyl", colmethdiff="dmp_median_methdiff"):
    dict_marker_gene = dict()
    for _, rows in df_dmp_deg_sorted.iterrows():
        dict_rows = dict(rows)
        direction = dict_rows[colmethdiff]
        if np.sign(direction) > 0:
            direction = "hyper"
        else:
            direction = "hypo"
        dmp_markers = ":".join(dict_rows[colmeth].split("_")[:-1])
        dict_marker_gene[dmp_markers] = direction
    
    return dict(dict_marker_gene)

def get_samplewise_genewise_cpg(outfilename, dict_markers, dict_sample_severity, dict_marker_gene):
    with open(outfilename, mode='w') as fw:
        list_header = ["ID", "Severity", "Marker", "GeneSymbol", "CpGbeta"]
        fw.write("\t".join(list_header) + "\n")
        for sample_id, dict_marker in dict_markers.items():
            if sample_id in dict_sample_severity.keys():
                severity_visit = dict_sample_severity[sample_id]
            else:
                continue
            for marker_name, cpg in dict_marker.items():
                pos = marker_name[0] + ":" + marker_name[1]
                if pos in dict_marker_gene.keys():
                    list_gene = dict_marker_gene[pos]
                    for gene in list_gene:
                        list_content = [sample_id, severity_visit, pos, gene, cpg]
                        fw.write("\t".join(list_content) + "\n")

def modify_samplewise_genewise_cpg_table(outfile):
    df_cpg_sample_all = pd.read_csv(outfile, sep="\t")
    df_cpg_sample_all["GeneSymbol_plot"] = df_cpg_sample_all["GeneSymbol"].apply(lambda x: "_".join(x.split("_")[1:]))
    df_cpg_sample_all["Marker_plot"] = df_cpg_sample_all["GeneSymbol_plot"] + " (" + df_cpg_sample_all["Marker"] + ")"
    df_cpg_sample_all["Severity_Group"] = df_cpg_sample_all["Severity"].str.split("_").str[0]
    df_cpg_sample_all["Visit"] = df_cpg_sample_all["Severity"].str.split("_").str[1]

    return df_cpg_sample_all

def get_list_cluster_markers(df_cpg_sample_gene_pivot, n_clust):
    import numpy as np
    from sklearn.cluster import AgglomerativeClustering
    data = df_cpg_sample_gene_pivot.values
    cluster = AgglomerativeClustering(n_clusters=n_clust, affinity='euclidean', linkage='average')
    list_clust = cluster.fit_predict(data)

    return list_clust

def get_drop_samples(path_drop_samples):
    with open(path_drop_samples, mode="r") as frmv:
        list_drop_samples = [x.rstrip() for x in frmv.readlines()]
    
    return list_drop_samples

def format_ticklabel(label):
    # Split text into parts: words outside and inside brackets
    parts = re.split(r'(\(.*?\))', label.get_text())  # This splits by brackets but keeps brackets in the list
    formatted_label = ""
    
    for part in parts:
        if part.startswith('(') and part.endswith(')'):
            # Do not italicize the part inside brackets
            formatted_label += part
        else:
            # Italicize the part outside the brackets
            formatted_label += f"$\\it{{{part}}}$ "  # Use LaTeX-style italicization
            
    return formatted_label

# %%
visit = "first"
pathsevinfo = "/BiO/Access/kyungwhan1998/Infectomics/Results/9_clinical/Infectomics_Severity_Information_20240926.tsv"
dirmethylcpgmin = "/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/MethylCpGMin"
# filetargetgene = f"/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/DMPDEG/{visit}/table_deg_dmp_overlap_abslog2fc_1.3_qval_0.05_sorted.tsv"
filetargetgene = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.DMP_DEG_LOO_annot.20240402.tsv"
dictmarkershyper = f"/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/Epigenetic_changes/{visit}/hyper/dictionary_marker_freqC_all_samples_20240220.pk.gz"
dictmarkershypo = f"/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/Epigenetic_changes/{visit}/hypo/dictionary_marker_freqC_all_samples_20240220.pk.gz"
outfilehyper = f"/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/Epigenetic_changes/{visit}/hyper/samplewise_genewise_cpg_dmpdeg_overlap_20240402.tsv"
outfilehypo = f"/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/Epigenetic_changes/{visit}/hypo/samplewise_genewise_cpg_dmpdeg_overlap_20240402.tsv"
# fileexp = "/BiO/Access/kyungwhan1998/Infectomics/Results/4_expmtx/ConfirmedRecovered/expression_matrix_genes.results_TPM.tsv" 
outfig = "/BiO/Access/kyungwhan1998/Infectomics/Results/InfectomicsPaper1/methylation_heatmap_severity_and_infectionphase.png"
list_drop_samples = ["C19-C009-V1",
                    "C19-C011-V1",
                    "C19-C016-V1",
                    "C19-C021-V1",
                    "C19-C022-V1",
                    "C19-C045-V2",
                    "C19-C045-V3",
                    "C19-C047-V2",
                    "C19-C047-V3",
                    "C19-C050-V2",
                    "C19-C050-V3",
                    "C19-C051-V2",
                    "C19-C051-V3",
                    "C19-C052-V2",
                    "C19-C052-V3",
                    "C19-C053-V2",
                    "C19-C053-V3",
                    "C19-C055-V2",
                    "C19-C055-V3",
                    "C19-C056-V2",
                    'C19-C056-V3',
                    'C19-C060-V2',
                    'C19-C060-V3',
                    'C19-C061-V2',
                    'C19-C061-V3']

# %%
list_files_methylcpgmin = get_list_files_methylcpgmin(dirmethylcpgmin)
list_methyl_sample = get_list_methyl_sample(list_files_methylcpgmin, list_drop_samples)
dict_markers_hyper = load_pickle(dictmarkershyper)
dict_markers_hypo = load_pickle(dictmarkershypo)
dict_sample_severity = get_dict_sample_severity(pathsevinfo, list_methyl_sample)
dict_sample_severity_sorted = sort_dict_sample_severity_visit(dict_sample_severity)
df_dmp_deg_sorted = get_dmp_deg_stats(filetargetgene)
dict_marker_gene = get_dict_target_marker_gene(df_dmp_deg_sorted)
dict_marker_direction = get_dict_target_marker_direction(df_dmp_deg_sorted)
get_samplewise_genewise_cpg(outfilehyper, dict_markers_hyper, dict_sample_severity_sorted, dict_marker_gene)
get_samplewise_genewise_cpg(outfilehypo, dict_markers_hypo, dict_sample_severity_sorted, dict_marker_gene)
df_cpg_sample_hyper = modify_samplewise_genewise_cpg_table(outfilehyper)
df_cpg_sample_hypo = modify_samplewise_genewise_cpg_table(outfilehypo)
df_cpg_sample_all = pd.concat([df_cpg_sample_hypo, df_cpg_sample_hyper], axis=0)

# %%
df_cpg_sample_gene_pivot_hyper = df_cpg_sample_hyper.pivot(index="Marker_plot", columns="ID", values="CpGbeta")
df_cpg_sample_gene_pivot_hyper["Cluster"] = get_list_cluster_markers(df_cpg_sample_gene_pivot_hyper, n_clust=3)
df_cpg_sample_gene_pivot_hyper_sorted = df_cpg_sample_gene_pivot_hyper.sort_values(by=["Cluster"], ascending=True)
df_cpg_sample_gene_pivot_hyper_sorted = df_cpg_sample_gene_pivot_hyper_sorted.drop(columns=["Cluster"])
df_cpg_sample_gene_pivot_hypo = df_cpg_sample_hypo.pivot(index="Marker_plot", columns="ID", values="CpGbeta")
df_cpg_sample_gene_pivot_hypo["Cluster"] = get_list_cluster_markers(df_cpg_sample_gene_pivot_hypo, n_clust=5)
df_cpg_sample_gene_pivot_hypo_sorted = df_cpg_sample_gene_pivot_hypo.sort_values(by=["Cluster"], ascending=True)
df_cpg_sample_gene_pivot_hypo_sorted = df_cpg_sample_gene_pivot_hypo_sorted.drop(columns=["Cluster"])
df_cpg_sample_gene_pivot = pd.concat([df_cpg_sample_gene_pivot_hypo_sorted, df_cpg_sample_gene_pivot_hyper_sorted], axis=0)

list_sample_sort = dict_sample_severity_sorted.keys()
df_cpg_sample_gene_pivot_sorted = df_cpg_sample_gene_pivot[list_sample_sort]
df_cpg_sample_gene_pivot_sorted["Pos"] = df_cpg_sample_gene_pivot_sorted.index.str.split("(").str[-1].str[:-1]
df_cpg_sample_gene_pivot_sorted.to_csv("/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/DMPDEG/first/cpgbeta_overlap.tsv", sep="\t", index=False)

# %%
list_sampleid = list(df_cpg_sample_all.groupby("ID")["Severity_Group"].apply(list).index)
list_sev = df_cpg_sample_all.groupby("ID")["Severity_Group"].apply(list).str[0].to_list()
list_vis = df_cpg_sample_all.groupby("ID")["Visit"].apply(list).str[0].to_list()
list_genes = df_cpg_sample_gene_pivot_sorted.index.to_list()
list_markers = df_cpg_sample_gene_pivot_sorted["Pos"].to_list()
df_cpg_sample_gene_pivot_sorted = df_cpg_sample_gene_pivot_sorted.drop(columns=["Pos"])
list_directions = list(map(lambda x: dict_marker_direction[x], list_markers))
dict_color_sev = {"Mild": "forestgreen", "Severe": "firebrick"}
dict_color_vis = {"Acute": "pink", "Recovery": "powderblue"}
dict_color_direction = {"hypo": "tab:blue", "hyper": "tab:red"}
list_sev_conv = list(map(lambda x: SEV_CONV.get(str(x), str(x)), list_sev))
list_vis_conv = list(map(lambda x: VISIT_CONV.get(str(x), str(x)), list_vis))
list_sev_color = list(map(dict_color_sev.__getitem__, list_sev_conv))
list_vis_color = list(map(dict_color_vis.__getitem__, list_vis_conv))
list_direction_color = list(map(dict_color_direction.__getitem__, list_directions))
color_df = pd.DataFrame({"Infection Phase": list_vis_color, "Severity Group": list_sev_color}, index=list_sampleid)
row_df = pd.DataFrame({"Methyl.\nStatus": list_direction_color}, index=list_genes)

# %%
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch

# Set figure fontdicterties
plt.rcParams["font.size"] = 10
cm = 1 / 2.54
width = 45 * cm
height = 65 * cm
figsize = (width, height)

g = sns.clustermap(df_cpg_sample_gene_pivot_sorted,
                   row_colors=row_df,
                   col_colors=color_df,
                   colors_ratio=0.02,
                   cmap="RdYlBu_r",
                   row_cluster=False,
                   col_cluster=False,
                   xticklabels=False,
                   yticklabels=True,
                   square=True,
                   vmin=0, 
                   vmax=100,
                   cbar_kws={"ticks": list(range(0, 110, 25)), 
                             "orientation": 'horizontal'})

yticklabels = g.ax_heatmap.get_yticklabels()
g.ax_heatmap.set_yticklabels([format_ticklabel(label) for label in yticklabels], fontdict={'fontsize': plt.rcParams["font.size"]+1}, rotation=0)
g.ax_cbar.remove()
g.ax_heatmap.set_xlabel("COVID-19 Patients", size=plt.rcParams["font.size"]+6)
g.ax_heatmap.set_ylabel("Methyation Markers", size=0)
g.ax_heatmap.axvline(x=37, color="white", linewidth=3, linestyle="solid")
g.ax_heatmap.axvline(x=46, color="white", linewidth=3, linestyle="solid")
g.ax_heatmap.axvline(x=83, color="white", linewidth=3, linestyle="solid")
g.ax_heatmap.axhline(y=0, color="white", linewidth=5, linestyle="solid")
g.ax_heatmap.axhline(y=28, color="white", linewidth=3, linestyle="solid")

g.gs.update(left=-0.1, right=1.1, top=1.2, bottom=0.55)

gs2 = gridspec.GridSpec(4, 1, left=-0.02, right=0.01, top=1.0, bottom=0.60, hspace=0.5)
ax_legend1 = g.fig.add_subplot(gs2[0])
ax_legend2 = g.fig.add_subplot(gs2[1])
ax_legend3 = g.fig.add_subplot(gs2[2])
ax_cbar1 = g.fig.add_subplot(gs2[3])

# Add legends for infection phase, severity, and methylation status
legend_elements1 = [Patch(facecolor='pink', edgecolor='k', label='Acute'),
                    Patch(facecolor='powderblue', edgecolor='k', label='Recovery')]
legend_elements2 = [Patch(facecolor='firebrick', edgecolor='k', label='Severe'),
                    Patch(facecolor='forestgreen', edgecolor='k', label='Mild')]
legend_elements3 = [Patch(facecolor='tab:blue', edgecolor='k', label='Hypomethylation'),
                    Patch(facecolor='tab:red', edgecolor='k', label='Hypermethylation')]

ax_legend1.text(-3.0, 2.3, "A", transform=ax_legend1.transAxes, fontsize=20, weight='bold', va='top', ha='left')
ax_legend1.legend(handles=legend_elements1, loc="center", frameon=False, fontsize=plt.rcParams["font.size"]+4)
ax_legend1.set_title("Infection Phase", fontdict={"size": plt.rcParams["font.size"]+6})
ax_legend1.axis('off')  # Hide the axis for the legend

ax_legend2.legend(handles=legend_elements2, loc="center", frameon=False, fontsize=plt.rcParams["font.size"]+4)
ax_legend2.set_title("Severity Group", fontdict={"size": plt.rcParams["font.size"]+6})
ax_legend2.axis('off')

ax_legend3.legend(handles=legend_elements3, loc="center", frameon=False, fontsize=plt.rcParams["font.size"]+4)
ax_legend3.set_title("Methylation Status", fontdict={"size": plt.rcParams["font.size"]+6})
ax_legend3.axis('off')

plt.colorbar(g.ax_heatmap.collections[0], cax=ax_cbar1, orientation='horizontal')
ax_cbar1.set_title('$\\beta$-value', fontdict={"fontsize":plt.rcParams["font.size"]+6})
ax_cbar1.tick_params(axis='x', length=0.5)
ax_cbar1.set_position([-0.08, 0.65, 0.15, 0.02])

from run_pca_methyl_rna import plot_pca_methyl, plot_pca_rna

path_meta_methyl = "/BiO/Access/kyungwhan1998/Infectomics/Results/9_clinical/Infectomics_Severity_Information_20240926.tsv"
col_meta_id = "Sample_ID"
col_meta_visit = "Visit_order"
col_meta_sev = "Severity"

path_methyl = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/MethylCpGTable/Infectomics.Copy_From_HjRyu/MethylCpGTable.Control.Mild.Case.Severe.Filtered.DMP.Hyper_Hypo.Sev_vs_Mild.Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.tsv"
cols_feat_methyl = ["chr", "start", "end"]

path_meta_rna = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/COVID19_master_table_20231007.Methyl_Overlap.with_Severity.20240402.txt"
path_rna_table = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/RNAExpressionTable/RNAExpression.COVID19.RNA_samples_with_Methyl.filter_Genes.Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.tsv"
col_rna_sampleid = "Project_ID_Alias"
col_rna_geneid = "Gene_ID"

gs5 = gridspec.GridSpec(1, 1, left=-0.025, right=0.40, top=0.4, bottom=0.01)
ax5 = g.fig.add_subplot(gs5[0])
ax_methyl = plot_pca_methyl(path_meta_methyl,
                            col_meta_id, 
                            col_meta_visit, 
                            col_meta_sev, 
                            path_methyl, 
                            cols_feat_methyl, 
                            ax5, 
                            list_pc = [1,2], 
                            fig_letter="B", 
                            letter_pos=(-0.2, 1.04))
ax_methyl.spines['right'].set_visible(False)
ax_methyl.spines['top'].set_visible(False)
ax_methyl.set_title(label="DNA methylation dynamics of COVID-19\npatients across infection phases", fontsize=plt.rcParams["font.size"]+6)
ax_methyl.tick_params(axis='y', labelsize=plt.rcParams["font.size"]+8)
ax_methyl.tick_params(axis='x', labelsize=plt.rcParams["font.size"]+8)
ax_methyl.legend().set_visible(False)

gs6 = gridspec.GridSpec(1, 1, left=0.60, right=1.0, top=0.4, bottom=0.01)
ax6 = g.fig.add_subplot(gs6[0])
ax_rna = plot_pca_rna(path_meta_rna,
                      col_rna_sampleid, 
                      col_meta_visit, 
                      col_meta_sev, 
                      ax6, 
                      list_pc = [1,2], 
                      fig_letter="C", 
                      letter_pos=(-0.2, 1.04))
ax_rna.spines['right'].set_visible(False)
ax_rna.spines['top'].set_visible(False)
ax_rna.set_title(label="Gene expression dynamics of COVID-19\npatients across infection phases", fontsize=plt.rcParams["font.size"]+6)
ax_rna.tick_params(axis='y', labelsize=plt.rcParams["font.size"]+8)
ax_rna.tick_params(axis='x', labelsize=plt.rcParams["font.size"]+8)
ax_rna.legend().set_visible(False)

import matplotlib.lines as mlines
import matplotlib.patches as mpatches

gs7 = gridspec.GridSpec(1, 1, left=1.1, right=1.3, top=0.4, bottom=0.01)
ax7 = g.fig.add_subplot(gs7[0])
legend_handles = []
for label, color in PALETTE.items():
    marker = MARKER[label]
    handle = mlines.Line2D([], [], color=color, marker=marker, linestyle='None', markersize=10, label=label)
    legend_handles.append(handle)

ax7.legend(handles=legend_handles, title="Groups", loc="center", fontsize=14, title_fontsize=16, frameon=False)
ax7.spines['right'].set_visible(False)
ax7.spines['top'].set_visible(False)
ax7.spines['left'].set_visible(False)
ax7.spines['bottom'].set_visible(False)
ax7.axis('off')

# Show the plot with the legend
plt.tight_layout()
plt.savefig("/BiO/Access/kyungwhan1998/Infectomics/Results/InfectomicsPaper1/20240906/Figure3.png", bbox_inches="tight", dpi=600)
plt.show()
plt.close()
# %%

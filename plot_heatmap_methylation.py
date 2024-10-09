# %%
import glob
import gzip
import os
import pickle
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import StandardScaler


# %%
def main(pathsevinfo, dirmethylcpgmin, filetargetgene, dictmarkershyper, dictmarkershypo, outfilehyper, outfilehypo, outfig, list_drop_samples):
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

    list_sampleid = list(df_cpg_sample_all.groupby("ID")["Severity_Group"].apply(list).index)
    list_sev = df_cpg_sample_all.groupby("ID")["Severity_Group"].apply(list).str[0].to_list()
    list_vis = df_cpg_sample_all.groupby("ID")["Visit"].apply(list).str[0].to_list()
    list_genes = df_cpg_sample_gene_pivot_sorted.index.to_list()
    list_markers = df_cpg_sample_gene_pivot_sorted["Pos"].to_list()
    df_cpg_sample_gene_pivot_sorted = df_cpg_sample_gene_pivot_sorted.drop(columns=["Pos"])
    list_directions = list(map(lambda x: dict_marker_direction[x], list_markers))
    dict_color_sev = {"Mild": "forestgreen", "Severe": "firebrick"}
    dict_color_vis = {"First": "pink", "Last": "powderblue"}
    dict_color_direction = {"hypo": "grey", "hyper": "black"}
    list_sev_color = list(map(dict_color_sev.__getitem__, list_sev))
    list_vis_color = list(map(dict_color_vis.__getitem__, list_vis))
    list_direction_color = list(map(dict_color_direction.__getitem__, list_directions))
    color_df = pd.DataFrame({"Infection Phase": list_vis_color, "Severity Group": list_sev_color}, index=list_sampleid)
    row_df = pd.DataFrame({"Methyl.\nStat.": list_direction_color}, index=list_genes)

    g = sns.clustermap(df_cpg_sample_gene_pivot_sorted,
                    row_colors=row_df,
                    col_colors=color_df,
                    colors_ratio=0.02,
                    cmap="RdYlBu_r", 
                    row_cluster=False, 
                    col_cluster=False,
                    xticklabels=False,  # Changed from True to False
                    yticklabels=True, 
                    square=True,
                    vmin=0, 
                    vmax=100, 
                    figsize=(20, 10), 
                    cbar_kws={"ticks":list(range(0, 110, 25)), "orientation":'horizontal'})  


    g.ax_col_colors.tick_params(axis='x', length=0)

    # Increase fontsize and bold labels for heatmap x-axis and y-axis
    g.ax_heatmap.tick_params(axis='y', labelsize=14)

    for label in g.ax_heatmap.get_yticklabels():
        label.set_fontweight('bold')

    # Position and style colorbar
    x0, _y0, _w, _h = g.cbar_pos
    g.ax_cbar.set_position([x0+0.7, -0.1, g.ax_row_dendrogram.get_position().width-0.1, 0.02])
    g.ax_cbar.set_title('$\\beta$-value', fontdict={"fontweight":"bold","fontsize":16})
    g.ax_cbar.tick_params(axis='x', length=0.5)

    # Axis labels and lines
    g.ax_heatmap.set_xlabel("COVID-19 Patients", size=20)
    g.ax_heatmap.set_ylabel("", size=20)
    g.ax_heatmap.axvline(x=0, color="white", linewidth=5, linestyle="solid")
    g.ax_heatmap.axvline(x=37, color="white", linewidth=3, linestyle="solid")
    g.ax_heatmap.axvline(x=46, color="white", linewidth=3, linestyle="solid")
    g.ax_heatmap.axvline(x=83, color="white", linewidth=3, linestyle="solid")
    g.ax_heatmap.axhline(y=0, color="white", linewidth=5, linestyle="solid")
    g.ax_heatmap.axhline(y=28, color="white", linewidth=3, linestyle="solid")


    from matplotlib.patches import Patch
    legend_elements1 = [Patch(facecolor='pink', edgecolor='k',label='Acute'),
                        Patch(facecolor='powderblue', edgecolor='k',label='Recovery')]
    legend_elements2 = [Patch(facecolor='firebrick', edgecolor='k',label='Severe'),
                        Patch(facecolor='forestgreen', edgecolor='k',label='Mild')]
    legend_elements3 = [Patch(facecolor='grey', edgecolor='k',label='Hypomethylation'),
                        Patch(facecolor='black', edgecolor='k',label='Hypermethylation')]

    legend1 = plt.legend(handles=legend_elements1, loc="center left", bbox_to_anchor=(-7.0, -0.1), frameon=False, fontsize=14)
    legend1.set_title("Infection Phase", prop={"weight":"bold", "size":16})  
    legend2 = plt.legend(handles=legend_elements2, loc="center left", bbox_to_anchor=(-5.0, -0.1), frameon=False, fontsize=14)
    legend2.set_title("Severity Group", prop={"weight":"bold", "size":16})  
    legend3 = plt.legend(handles=legend_elements3, loc="center left", bbox_to_anchor=(-3.0, -0.1), frameon=False, fontsize=14)
    legend3.set_title("Methylation Status", prop={"weight":"bold", "size":16})  
    legend3.get_frame().set_alpha(None)
    plt.gca().add_artist(legend1)
    plt.gca().add_artist(legend2)
    # plt.gca().axes.xaxis.set_ticklabels([])
    # plt.rcParams['xtick.bottom'] = False
    # plt.rcParams['xtick.labelbottom'] = False
    plt.tight_layout()
    # plt.savefig(outfig, dpi=600, bbox_inches="tight")
    plt.show()
    plt.close()
    
    # ------

    # list_target_gene_sort = df_dmp_deg_sorted["Gene_Symbol"].to_list()
    # list_idx_gene_sort = list(range(len(list_target_gene_sort)))
    # dict_gene_sort = dict(zip(list_target_gene_sort, list_idx_gene_sort))
    
    # df_exp = pd.read_csv(fileexp, sep="\t")
    # df_exp["C19-C014-V1"] = np.nan
    # df_exp_idxed = df_exp.set_index("ID")
    # list_gene_id = df_exp_idxed.index
    # list_genesymbol = list(map(lambda x: "_".join(x.split("_")[1:]), list_gene_id))
    # list_all_samples = df_exp_idxed.columns
    # mtx_exp_scaled = StandardScaler().fit_transform(df_exp_idxed)
    # df_exp_scaled = pd.DataFrame(mtx_exp_scaled, index=list_gene_id, columns=list_all_samples)
    # list_expsamples = list(df_exp_scaled.columns)
    # list_colsamples = list(filter(lambda x: x in list_expsamples, list_sample_sort))
    # df_exp_selec = df_exp_scaled[list_colsamples]
    # df_exp_selec["gene_symbol"] = list_genesymbol
    # df_exp_selec_filt = df_exp_selec[df_exp_selec["gene_symbol"].isin(list_target_gene_sort)]
    # df_exp_selec_filt.index = list(map(lambda x: "_".join(x.split("_")[1:]) + "(" + x.split("_")[0] + ")", df_exp_selec_filt.index))
    # df_exp_selec_filt["Order"] = df_exp_selec_filt["gene_symbol"].map(dict_gene_sort)
    # df_exp_selec_filt_order = df_exp_selec_filt.sort_values(by=["Order"], ascending=True)
    # df_exp_selec_filt_order_rmv_dup = df_exp_selec_filt_order[~df_exp_selec_filt_order["gene_symbol"].duplicated(keep="first")]
    # df_exp_selec_filt_drop_genesym = df_exp_selec_filt_order_rmv_dup.drop(columns=["gene_symbol", "Order"])

    # sns.clustermap(df_exp_selec_filt_drop_genesym, cmap="rocket_r", row_cluster=False, col_cluster=False, yticklabels=True, xticklabels=True, linewidths=.5, square=True, robust=True, mask=df_exp_selec_filt_drop_genesym.isnull(), figsize=(20, 8))
    # plt.show()
    # plt.close()


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

# %%
if __name__ == "__main__":
    visit = "first"
    pathsevinfo = "/BiO/Access/kyungwhan1998/Infectomics/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"
    dirmethylcpgmin = "/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/MethylCpGMin"
    # filetargetgene = f"/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/DMPDEG/{visit}/table_deg_dmp_overlap_abslog2fc_1.3_qval_0.05_sorted.tsv"
    filetargetgene = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/Methyl_RNA_Relation/Methyl_RNA_Correlation.Infectomics.Visit1_only.Mild_Severe.Methyl_perc.RNA_DESeqNormcount.Methyl_Filtered.DMP_Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.DMP_DEG_LOO_annot.20240402.tsv"
    dictmarkershyper = f"/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/Epigenetic_changes/{visit}/hyper/dictionary_marker_freqC_all_samples_20240220.pk.gz"
    dictmarkershypo = f"/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/Epigenetic_changes/{visit}/hypo/dictionary_marker_freqC_all_samples_20240220.pk.gz"
    outfilehyper = f"/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/Epigenetic_changes/{visit}/hyper/samplewise_genewise_cpg_dmpdeg_overlap_20240402.tsv"
    outfilehypo = f"/BiO/Access/kyungwhan1998/Infectomics/Results/10_methyl/Epigenetic_changes/{visit}/hypo/samplewise_genewise_cpg_dmpdeg_overlap_20240402.tsv"
    # fileexp = "/BiO/Access/kyungwhan1998/Infectomics/Results/4_expmtx/ConfirmedRecovered/expression_matrix_genes.results_TPM.tsv" 
    outfig = "/BiO/Access/kyungwhan1998/Infectomics/Results/InfectomicsPaper1/methylation_heatmap_severity_and_infectionphase.png"
    list_drop_samples = [
                        "C19-C009-V1",
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
                        'C19-C061-V3'
                        ]
        
    main(pathsevinfo, dirmethylcpgmin, filetargetgene, dictmarkershyper, dictmarkershypo, outfilehyper, outfilehypo, outfig, list_drop_samples)

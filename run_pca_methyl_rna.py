#%%
import itertools
from collections import defaultdict
from functools import partial

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Ellipse
from scipy.stats import sem
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

VISIT_CONV = {
    "Visit0" : "Healthy",
    "Visit1" : "Acute",
    "Visit3" : "Recovery",
    "Visit4" : "Recovery",
    "Visit5" : "Convalescent"
}
SEV_CONV = {
    "0" : "Healthy",
    "1" : "Mild",
    "2" : "Mild",
    "3" : "Severe",
    "4" : "Severe"
}
PALETTE = {
    "Healthy": (0.6627450980392157, 0.6627450980392157, 0.6627450980392157, 1.0),
    "Convalescent" : (0.2549019607843137, 0.4117647058823529, 0.8823529411764706, 1.0),
    "Mild(Recovery)" : (0.13333333333333333, 0.5450980392156862, 0.13333333333333333, 1.0),
    "Mild(Acute)" : (0.13333333333333333, 0.5450980392156862, 0.13333333333333333, 1.0),
    "Severe(Recovery)" : (0.6980392156862745, 0.13333333333333333, 0.13333333333333333, 1.0),
    "Severe(Acute)" : (0.6980392156862745, 0.13333333333333333, 0.13333333333333333, 1.0)
}
MARKER = {
    "Healthy" : "D",
    "Convalescent" : "D",
    "Mild(Recovery)" : "^",
    "Mild(Acute)" : "o",
    "Severe(Recovery)" : "^",
    "Severe(Acute)" : "o"
}

# %%
class pca_rna:
    explained_variance_ratio_ = [0.36382185, 0.17454983, 0.09582511, 0.06211652, 0.04747781]

def reshape_rna_table(table_rna, col_sampleid, col_geneid, col_exp):
    table_rna_sorted = table_rna.sort_values(by = col_geneid)
    
    list_samples = list(table_rna[col_sampleid].unique())
    list_genes = table_rna_sorted[table_rna_sorted[col_sampleid] == list_samples[0]][col_geneid].to_list()
    table_rna_sorted = table_rna_sorted.set_index(col_geneid, drop = False)
    
    table_rna_reshaped = pd.DataFrame(index = list_genes, columns = list_samples)
    for sample in list_samples:
        table_rna_sample = table_rna_sorted[table_rna_sorted[col_sampleid] == sample]
        table_rna_reshaped.loc[:, sample] = table_rna_sample[col_exp]
    table_rna_reshaped_dropna = table_rna_reshaped.dropna()
    return table_rna_reshaped_dropna.T

def run_pca(data, n_comp = 2):
    pca = PCA(n_components=n_comp).fit(data)
    transformed_data = pca.fit_transform(data)
    table_transformed = pd.DataFrame(transformed_data)
    table_transformed.index = data.index
    table_transformed.columns = list(map(lambda x : f"PC{x}", range(1, n_comp+1)))
    return table_transformed, pca

def add_meta_info_for_transformed_table(table_transformed, table_meta, col_meta_id, cols_add):
    table_transformed["Sample_ID"] = table_transformed.index
    for col_add in cols_add:
        dict_sample_to_col = dict(zip(table_meta[col_meta_id], table_meta[col_add]))
        table_transformed[col_add] = table_transformed["Sample_ID"].apply(dict_sample_to_col.__getitem__)
    return table_transformed

def combine_visit_order_and_severity_info(visit, sev):
    visit_conv = VISIT_CONV.get(str(visit), str(visit))
    sev_conv = SEV_CONV.get(str(sev), str(sev))
    
    if visit_conv == "Convalescent":
        return visit_conv
    
    elif visit_conv == "Healthy":
        return visit_conv
    
    else:
        return f"{sev_conv}({visit_conv})"    

def get_indexes_of_list_from_other_list(input, other = list()):
    return list(map(lambda val : other.index(val), input.to_list()))

from matplotlib.patches import FancyArrowPatch


def plot_pca(table_transformed, pca_obj, ax, list_pc=[1, 2], hue = "Severity(Phase)"):
    partial_get_indexes_of_list_from_hue_order = partial(get_indexes_of_list_from_other_list, other = list(PALETTE.keys()))
    table_transformed = table_transformed.sort_values(by = hue, key = partial_get_indexes_of_list_from_hue_order)
    pcs_sorted = sorted(list_pc)
    dict_pc_mean = defaultdict(dict)
    for ind_pc1, pc1 in enumerate(pcs_sorted[:-1]):
        for pc2 in pcs_sorted[ind_pc1+1:]:
            sns.scatterplot(data = table_transformed, x = f"PC{pc1}", y = f"PC{pc2}", hue = hue, ax = ax, palette = PALETTE, style = hue, markers = MARKER, s = 0, alpha = 1, edgecolor = "k")
            ax.set_xlabel(f"principal component {pc1}\n{format(pca_obj.explained_variance_ratio_[pc1-1]*100, '.2f')}%")
            ax.set_ylabel(f"principal component {pc2}\n{format(pca_obj.explained_variance_ratio_[pc2-1]*100, '.2f')}%")
            for hue_name in table_transformed[hue].unique():
                table_transformed_hue = table_transformed[table_transformed[hue] == hue_name]
                list_x = table_transformed_hue[f"PC{pc1}"].to_numpy()
                list_y = table_transformed_hue[f"PC{pc2}"].to_numpy()
                color = PALETTE[hue_name]
                confidence_ellipse(list_x, list_y, ax, edgecolor = color, linewidth=2, zorder=-1) 
                center_x = np.mean(list_x)
                center_y = np.mean(list_y)
                marker = MARKER[hue_name]
                ax.scatter(x=center_x, y=center_y, s=150, facecolor=color, marker=marker, zorder=-1)
                if "Severe" in hue_name or "Conv" in hue_name:
                    dict_pc_mean["Severe"].update({hue_name: (center_x, center_y)})
                if "Mild" in hue_name or "Conv" in hue_name:
                    dict_pc_mean["Mild"].update({hue_name: (center_x, center_y)})
    
    S = list(dict_pc_mean["Severe"].values()) 
    list_iter_S = list(itertools.combinations(S, r=2))
    list_iter_Sx = [list(zip(*itS))[0] for itS in list_iter_S]
    list_iter_Sx_select = list([list_iter_Sx[2]])
    list_iter_Sy = [list(zip(*itS))[1] for itS in list_iter_S]
    list_iter_Sy_select = list([list_iter_Sy[2]])
    for sx, sy in zip(list_iter_Sx_select, list_iter_Sy_select):
        arrow_s = FancyArrowPatch((sx[1], sy[1]), (sx[0], sy[0]),
                                arrowstyle="Simple,head_length=8,head_width=8,tail_width=2",
                                # connectionstyle="arc3,rad=.5",
                                color="k", linewidth=0.1, zorder=999)
        ax.add_patch(arrow_s)

    M = list(dict_pc_mean["Mild"].values()) 
    list_iter_M = list(itertools.combinations(M, r=2))
    list_iter_Mx = [list(zip(*itM))[0] for itM in list_iter_M]
    list_iter_Mx_select = list([list_iter_Mx[2]])
    list_iter_My = [list(zip(*itM))[1] for itM in list_iter_M]
    list_iter_My_select = list([list_iter_My[2]])
    for mx, my in zip(list_iter_Mx_select, list_iter_My_select):
        arrow_m = FancyArrowPatch((mx[1], my[1]), (mx[0], my[0]),
                                arrowstyle="Simple,head_length=8,head_width=8,tail_width=2",
                                # connectionstyle="arc3,rad=.5",
                                color="k", linewidth=0.1, zorder=999)
        ax.add_patch(arrow_m)
        
    return ax

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)

    return ax.add_patch(ellipse)

def convert_rgb_to_cmyk(r, g, b, CMYK_SCALE = 100, RGB_SCALE = 255):
    if (r, g, b) == (0, 0, 0):
        # black
        return 0, 0, 0, CMYK_SCALE

    # rgb [0,255] -> cmy [0,1]
    c = 1 - r / RGB_SCALE
    m = 1 - g / RGB_SCALE
    y = 1 - b / RGB_SCALE

    # extract out k [0, 1]
    min_cmy = min(c, m, y)
    c = (c - min_cmy) / (1 - min_cmy)
    m = (m - min_cmy) / (1 - min_cmy)
    y = (y - min_cmy) / (1 - min_cmy)
    k = min_cmy

    # rescale to the range [0,CMYK_SCALE]
    return c * CMYK_SCALE, m * CMYK_SCALE, y * CMYK_SCALE, k * CMYK_SCALE

def convert_cmyk_to_rgb(c, m, y, k, CMYK_SCALE = 100, RGB_SCALE=255):
    r = RGB_SCALE * (1.0 - c / float(CMYK_SCALE)) * (1.0 - k / float(CMYK_SCALE))
    g = RGB_SCALE * (1.0 - m / float(CMYK_SCALE)) * (1.0 - k / float(CMYK_SCALE))
    b = RGB_SCALE * (1.0 - y / float(CMYK_SCALE)) * (1.0 - k / float(CMYK_SCALE))
    
    return r, g, b

def plot_pca_rna(path_meta, col_meta_id, col_meta_visit, col_meta_sev, ax, list_pc = [1,2], fig_letter="A", letter_pos=(0, 0)):
    table_meta = pd.read_csv(path_meta, sep = '\t')
    table_rna_transformed = pd.read_csv("/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/DEGExtract/DEGExtract.RNA_samples_with_Methyl.DEG_by_Sex.Control_M.Case_F.Cov_Age.mincount_1.18_genes_only.20240402.pca.pc12345.tsv", sep = '\t')
    table_rna_transformed = table_rna_transformed.set_index("name", drop = True)
    table_rna_transformed_annotated = add_meta_info_for_transformed_table(table_rna_transformed, table_meta, col_meta_id, [col_meta_visit, col_meta_sev])
    table_rna_transformed_annotated["Visit_order"] = table_rna_transformed_annotated["Visit_order"].apply(lambda x: "Visit5" if x == "Recover" else x)
    table_rna_transformed_annotated["Severity(Phase)"] = table_rna_transformed_annotated.apply(lambda row : combine_visit_order_and_severity_info(row[col_meta_visit], row[col_meta_sev]), axis = 1)
    ax_rna = plot_pca(table_rna_transformed_annotated, pca_rna(), ax, list_pc)
    ax_rna.set_xlim(-12, 7)
    ax_rna.set_ylim(-5.1, 5.1)
    ax_rna.set_xlabel(ax_rna.get_xlabel(), fontsize=16)
    ax_rna.set_ylabel(ax_rna.get_ylabel(), fontsize=16)
    ax_rna.set_title("Gene Expression", fontsize=16)
    ax_rna.annotate(fig_letter,
                    xy=letter_pos,  # Adjusted further left to avoid overlap
                    xycoords='axes fraction',
                    xytext=(0, 0),
                    textcoords='offset points',
                    size=plt.rcParams["font.size"]+10, 
                    ha='left', 
                    va='center',
                    fontweight="bold",
                    color="black")

    return ax_rna

def plot_pca_methyl(path_meta, col_meta_id, col_meta_visit, col_meta_sev, path_methyl, cols_feat_methyl, ax, list_pc = [1,2], fig_letter="B", letter_pos=(0, 0)):
    table_meta = pd.read_csv(path_meta, sep = '\t')
    table_methyl = pd.read_csv(path_methyl, sep = '\t')
    table_methyl["Feat"] = table_methyl[cols_feat_methyl].apply(lambda row : '_'.join(list(map(str, row))), axis = 1)
    table_methyl = table_methyl.drop(columns = cols_feat_methyl)
    table_methyl_reshaped = table_methyl.set_index("Feat", drop = True).T
    list_methyl_sample = list(table_methyl_reshaped.index)
    list_meta_sample = table_meta[col_meta_id].to_list()
    list_intsc = list(set(list_methyl_sample).intersection(set(list_meta_sample)))
    table_methyl_reshaped = table_methyl_reshaped.loc[list_intsc, :]
    scaler = StandardScaler().fit(table_methyl_reshaped)
    list_samples = table_methyl_reshaped.index
    list_cpgs = table_methyl_reshaped.columns
    table_methyl_reshaped = scaler.transform(table_methyl_reshaped)
    table_methyl_reshaped = pd.DataFrame(table_methyl_reshaped)
    table_methyl_reshaped.index = list_samples
    table_methyl_reshaped.columns = list_cpgs
    table_methyl_transformed, pca_methyl = run_pca(table_methyl_reshaped, max(list_pc))
    table_methyl_transformed_annotated = add_meta_info_for_transformed_table(table_methyl_transformed, table_meta, col_meta_id, [col_meta_visit, col_meta_sev])
    table_methyl_transformed_annotated["Severity(Phase)"] = table_methyl_transformed_annotated.apply(lambda row : combine_visit_order_and_severity_info(row[col_meta_visit], row[col_meta_sev]), axis = 1)
    ax_methyl = plot_pca(table_methyl_transformed_annotated, pca_methyl, ax, list_pc)
    ax_methyl.set_xlabel(ax_methyl.get_xlabel(), fontsize=16)
    ax_methyl.set_ylabel(ax_methyl.get_ylabel(), fontsize=16)
    ax_methyl.set_title("DNA Methylation", fontsize=16)
    ax_methyl.set_xlim(-25, 20)
    ax_methyl.set_ylim(-5, 5)
    ax_methyl.legend(bbox_to_anchor=[1.2, 0.5], loc='center', fontsize=12, frameon=False).set_visible(True)
    ax_methyl.annotate(fig_letter,
                xy=letter_pos,  # Adjusted further left to avoid overlap
                xycoords='axes fraction',
                xytext=(0, 0),
                textcoords='offset points',
                size=plt.rcParams["font.size"]+10, 
                ha='left', 
                va='center',
                fontweight="bold",
                color="black")
    
    return ax_methyl

# %%
if __name__ == "__main__":
    path_meta = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Resources/MetaTable/COVID19_master_table_20231007.Methyl_Overlap.with_Severity.20240402.txt"
    col_meta_id = "Project_ID_Alias"
    col_meta_visit = "Visit_order"
    col_meta_sev = "Severity"
    
    path_methyl = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/MethylCpGTable/Infectomics.Copy_From_HjRyu/MethylCpGTable.Control.Mild.Case.Severe.Filtered.DMP.Hyper_Hypo.Sev_vs_Mild.Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.tsv"
    cols_feat_methyl = ["chr", "start", "end"]
    
    path_rna_table = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Methylation/Marker_Selection_Severe_Mild_DMP_Low_Cutoff/Results/RNAExpressionTable/RNAExpression.COVID19.RNA_samples_with_Methyl.filter_Genes.Mild_Sev_Visit1.LOO_common.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_6.tsv"
    col_rna_sampleid = "Sample_ID"
    col_rna_geneid = "Gene_ID"
    col_rna_exp = "TPM"

    nrow = 1
    ncol = 2
    figsize=(15, 5)

    fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=figsize)  
    # Flatten axes to handle the case when it's a 2D array
    if nrow > 1 or ncol > 1:
        axes = axes.flatten()
    else:
        axes = np.array([axes]).flatten()

    ax_rna = plot_pca_rna(path_meta, col_meta_id, col_meta_visit, col_meta_sev, axes[0], list_pc = [1,2], fig_letter="A", letter_pos=(-0.15, 1.03))
    ax_methyl = plot_pca_methyl(path_meta, col_meta_id, col_meta_visit, col_meta_sev, path_methyl, cols_feat_methyl, axes[1], list_pc = [1,2], fig_letter="B", letter_pos=(-0.15, 1.03))
    
    plt.tight_layout()
    plt.show()
    plt.close()

# %%
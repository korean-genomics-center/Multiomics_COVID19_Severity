# %%
import os
import pandas as pd
from pathlib import Path
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from statannotations.Annotator import Annotator

# %%
WORKDIR = str(Path(__file__).parents[3])
path_sev_info = f"{WORKDIR}/Results/9_clinical/Infectomics_Severity_Information_Methyl_20240102.tsv"
dir_methylcpgmin = f"{WORKDIR}/Results/10_methyl/MethylCpGMin"
infilenamehyper = f"{WORKDIR}/Results/10_methyl/Epigenetic_changes/first/hyper/dictionary_marker_freqC_all_samples_20240220.pk.gz"
infilenamehypo = f"{WORKDIR}/Results/10_methyl/Epigenetic_changes/first/hypo/dictionary_marker_freqC_all_samples_20240220.pk.gz"
file_deg_dmp_overlap = f"{WORKDIR}/Results/11_dmp/Yes_LOO/Covariate_Sex_Age/Methyl_RNA_Correlation.Filtered.DMP_Covariate_Sex_Age_Yes_LOO_common_across_9folds_or_more.cis.distance_cf.1000000.Overlap_DEG.corr_0.5.log2fc_1.3.loo_7.20240326.tsv"
path_drop_samples = f"{WORKDIR}/Resources/Scripts/Final/list_remove_samples.txt"
outdir = f"{WORKDIR}/Results/Paper"
os.makedirs(outdir, exist_ok=True)
path_methyl_target = f"{WORKDIR}/Results/11_dmp/Yes_LOO/Covariate_Sex_Age/list_new_methyl_markers.txt"


# %%
from methylation_severity import (get_dict_deg_dmp_overlap_markers,
                                  main,
                                  make_dataframe_normality_test,
                                  make_dataframe_stat_test,
                                  save_norm_test_result, 
                                  save_stat_test_result,
                                  plot_methyl_difference)

df_beta_all_hyper, df_beta_all_hypo = main(path_sev_info, dir_methylcpgmin, infilenamehyper, infilenamehypo)
list_metainfo = list(df_beta_all_hyper.iloc[:, -8:].columns)
list_marker_hyper = list(filter(lambda x: x not in list_metainfo, list(df_beta_all_hyper.columns)))
df_beta_all_hyper_hypo = pd.concat([df_beta_all_hyper[list_marker_hyper], df_beta_all_hypo], axis=1)
dict_markers_overlap = get_dict_deg_dmp_overlap_markers(file_deg_dmp_overlap)

marker_overlap_hyper = list(set(dict_markers_overlap.keys()).intersection(set(list(df_beta_all_hyper.columns))))
marker_overlap_hypo = list(set(dict_markers_overlap.keys()).intersection(set(list(df_beta_all_hypo.columns))))

with open(path_methyl_target, mode="r") as fr:
    list_target_methyl = list(map(lambda x: x.rstrip(), fr.readlines()))

with open(path_drop_samples, mode="r") as fr:
    list_drop_samples = list(map(lambda x: x.rstrip(), fr.readlines()))

dict_marker_target = dict()
for marker, markers in dict_markers_overlap.items():
    if marker in list_target_methyl:
        target = "-".join([x.split("_")[-1] for x in markers])
        dict_marker_target[marker] = target
     
df_beta_selec_target = df_beta_all_hyper_hypo[list_target_methyl + list_metainfo]

# Paths to the marker lists
path_hypo = "/BiO/Access/kyungwhan1998/Infectomics/Results/11_dmp/Yes_LOO/Covariate_Sex_Age/list_new_methyl_markers_hypo.txt"
path_hyper = "/BiO/Access/kyungwhan1998/Infectomics/Results/11_dmp/Yes_LOO/Covariate_Sex_Age/list_new_methyl_markers_hyper.txt"

# Get marker list
def get_list_marker(path_marker):
    with open(path_marker, mode="r") as fr:
        list_marker = list(map(lambda x: x.rstrip("\n"), fr.readlines()))
    return list_marker

# Load hypo and hyper markers
list_marker_hypo = get_list_marker(path_hypo)
list_marker_hyper = get_list_marker(path_hyper)

# Filter the dataframe based on severity
def get_filtered_data(severity, list_marker):
    df_filtered = df_beta_selec_target[ 
        (df_beta_selec_target["Severity"] == int(severity)) & 
        (df_beta_selec_target["Visit"] != "Convalescent")
    ]
    return df_filtered[["Subject_ID", "Visit", "Severity"] + list_marker]

# Create figure with gridspec
fig = plt.figure(figsize=(12, 6))
# Adjusting the height_ratios to make A and C larger
gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[1.3, 1], height_ratios=[1, 1], wspace=0.3, hspace=0.5)

# Conditions for subplots (Severity, Marker type)
conditions = [
    ("3", list_marker_hypo, "A", -0.1),
    ("4", list_marker_hypo, "B", -0.13),
    ("3", list_marker_hyper, "C", -0.1),
    ("4", list_marker_hyper, "D", -0.13)
]

# Loop through conditions and create plots
for i, (severity, markers, label, label_xpos) in enumerate(conditions):
    # Get grid position for the current subplot grid (2x2 grid)
    ax = fig.add_subplot(gs[i])

    # Remove x and y axis for the overall grid section (A, B, C, D)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.axis('off')

    # Adjust the label positioning (move the label a bit to the left to avoid overlap)
    ax.text(label_xpos, 1.08, label, transform=ax.transAxes, fontsize=16, weight="bold", va='top', ha='right')

    df_filtered = get_filtered_data(severity, markers)

    # Subjects with paired visits
    subjects = df_filtered["Subject_ID"].value_counts()
    paired_subjects = subjects[subjects > 1].index.tolist()

    # Number of rows and columns in the current grid for subjects
    n_subjects = len(paired_subjects)
    nrows = 1  # Number of rows (adjust if necessary)
    ncols = int(np.ceil(n_subjects / nrows))  # Number of columns

    # Create a subplot for each subject within the subgridspec
    subgs = gs[i].subgridspec(nrows, ncols, hspace=1.5, wspace=0.3)  # Increase hspace for more vertical spacing

    for j, subject_id in enumerate(paired_subjects):
        # Create a new subplot for each subject within the subgridspec
        ax_subject = fig.add_subplot(subgs[j])

        # Filter data for the current subject
        df_subject = df_filtered[df_filtered["Subject_ID"] == subject_id]

        # Prepare melt data for plotting
        df_melt = df_subject[["Visit"] + markers].set_index("Visit").T.reset_index()
        df_melt.columns.name = None
        df_melt = df_melt.rename(columns={"index": "Genomic_Pos"})
        df_melt = df_melt.melt(id_vars="Genomic_Pos", var_name="Visit", value_name="Methyl")
        palette = {"First": "#F4A7B9", "Last": "#ADD8E6"}
        # Plot boxplot with light grey color and turn off outlier dots
        sns.boxplot(data=df_melt, x="Visit", y="Methyl", ax=ax_subject, width=0.4, showfliers=False, color='lightgrey', zorder=2, palette=palette)

        # Add scatter plot with white facecolor and black edgecolor, with narrower edge
        sns.stripplot(data=df_melt, x="Visit", y="Methyl", ax=ax_subject, color='white', marker='o', size=4, facecolor=None, edgecolor='black', linewidth=0.5, jitter=False, alpha=0.6, zorder=2, palette=palette, label=None)

        # Add the missing lineplot
        sns.lineplot(data=df_melt, x="Visit", y="Methyl", estimator=None, units="Genomic_Pos", color="grey", alpha=0.6, linewidth=0.7, ax=ax_subject, zorder=2, palette=palette, label=None)

        # Add annotation
        pairs = [("First", "Last")]
        annotator = Annotator(ax_subject, pairs, data=df_melt, x="Visit", y="Methyl")
        annotator.configure(test='t-test_paired', text_format='star', loc='inside', verbose=0)
        annotator.apply_and_annotate()

        # Set the title for the subplot (subject label)
        ax_subject.set_title(subject_id, fontsize=10, weight="bold")

        # Add grid on the y-axis
        ax_subject.yaxis.grid(True)

        # Rotate xticks for better readability
        ax_subject.set_xticklabels(["Acute", "Recovery"], rotation=45, ha='right', rotation_mode="anchor")
        ax_subject.set_xlabel("", fontsize=0)
        
        # Set y-tick labels and ylabel only on the first subplot of each grid (A, C)
        if j == 0:  # First subplot in the subgrid (A or C)
            ax_subject.set_ylabel("Methylation (%)", fontsize=10)
            ax_subject.set_ylim(-1, 109)
        else:  # No y-tick labels for other subplots
            ax_subject.set_ylabel("", fontsize=0)
            ax_subject.set_yticklabels([])

# Adjust layout and remove unused axes
plt.tight_layout()  # Adjust padding between subplots
plt.subplots_adjust()  # Make room at the top and bottom if necessary
plt.savefig(os.path.join(outdir, "SupplementaryFigure8.png"), dpi=300)
plt.savefig(os.path.join(outdir, "SupplementaryFigure8.pdf"), dpi=300)
plt.show()
plt.close()
# %%

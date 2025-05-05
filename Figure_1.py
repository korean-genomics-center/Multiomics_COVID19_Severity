# %%
import math
import os
import warnings
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import clinical
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


# %% [clinicals]
def main(path_ehr, path_drop_samples, path_rename, outfigname):
    list_drop_samples = get_drop_samples(path_drop_samples)
    parse_clin = clinical.Parse(path_ehr, list_drop_samples, path_rename)
    df_crf = parse_clin.give_parsed_dataframe()
    df_crf["Subject_ID"] = df_crf["Sample_ID"].apply(lambda x: "-".join(x.split("-")[:-1]))
    df_crf["Enrolled"] = df_crf["Enrolled"].apply(lambda x: x.replace("\n", "") if str(x).startswith("\n") else x)
    df_crf_columns_target_fillna = fillna_by_first_row(df_crf)
    df_crf_columns_target_fillna = fillna_no_values_by_test_date(df_crf_columns_target_fillna)
    df_crf_columns_target_fillna = df_crf_columns_target_fillna.applymap(convert_timestamp_to_datetime)
    
    df_subj_timeline = df_crf_columns_target_fillna.copy()
    df_subj_timeline = df_subj_timeline.set_index("Sample_ID")
    list_filt_type_col = list(filter(lambda x: str(type(df_subj_timeline[x][0]))=="<class 'datetime.date'>", list(df_subj_timeline.columns)))
    lit_filt_col = ["Severity_group"] + list_filt_type_col
    df_subj_timeline = df_subj_timeline[lit_filt_col]
    df_subj_timeline = df_subj_timeline.reset_index(drop=False)
    df_subj_timeline = df_subj_timeline[df_subj_timeline["Sample_ID"].str.startswith("C19-C")]
    df_subj_timeline = df_subj_timeline[df_subj_timeline["Sample_ID"].str.endswith("V1")]
    
    list_raw_timeline = get_raw_timeline(df_subj_timeline)
    list_raw_timeline = sorted(list_raw_timeline, key=lambda x: x[1])
    list_timeline = get_nafilt_timeline(list_raw_timeline)
    list_list_ind = get_list_indices_no_na(list_raw_timeline)
    list_annot = get_annotation(df_subj_timeline, list_list_ind)
    
    plot_timeline(list_timeline, list_annot, outfigname)

def get_drop_samples(path_drop_samples):
    with open(path_drop_samples, mode="r") as frmv:
        list_drop_samples = [x.rstrip() for x in frmv.readlines()]
    
    return list_drop_samples

def fillna_by_first_row(df_crf_columns_target):
    df_crf_columns_target_fillna = df_crf_columns_target.copy()
    list_columns = list(df_crf_columns_target.columns)
    list_annot = list_columns[2:-1]
    for annot in list_annot:
        df_crf_columns_target_fillna[annot] = df_crf_columns_target[annot].fillna(df_crf_columns_target.groupby("Subject_ID")[annot].transform("first"))

    return df_crf_columns_target_fillna

def fillna_no_values_by_test_date(df_crf_columns_target):
    df_crf_columns_target_fillna = df_crf_columns_target.copy()
    df_crf_columns_target_fillna["StartOxygenTreatment"] = df_crf_columns_target_fillna["StartOxygenTreatment"].fillna(value=df_crf_columns_target_fillna["TestedPositive"])
    df_crf_columns_target_fillna["EndOxgenTreatment"] = df_crf_columns_target_fillna["EndOxgenTreatment"].fillna(value=df_crf_columns_target_fillna["TestedPositive"])
    
    return df_crf_columns_target_fillna

def convert_timestamp_to_datetime(x):
    if type(x) == pd._libs.tslibs.timestamps.Timestamp:
        timestamp = str(x)
        convert_x = datetime.strptime(str(timestamp), "%Y-%m-%d %H:%M:%S").date()

    else:
        convert_x = x
    
    return convert_x

def get_raw_timeline(df_subj_timeline):
    list_raw_timeline = list()
    for _, row in df_subj_timeline.iterrows():
        list_row = list(row)
        list_raw_timeline.append(list_row)
    
    return list_raw_timeline

def get_nafilt_timeline(list_raw_timeline):
    list_nafilt_timeline = list()
    for timeline in list_raw_timeline:
        list_elem = list()
        for idx, elem in enumerate(timeline):
            if str(elem) != "nan":
                list_elem.append(elem)
        
        list_nafilt_timeline.append(list_elem)
    
    return list_nafilt_timeline

def get_list_indices_no_na(list_raw_timeline):
    list_list_idx = list()
    for timeline in list_raw_timeline:
        list_idx = list()
        for idx, elem in enumerate(timeline):
            if str(elem) != "nan":
                list_idx.append(idx)
        
        list_list_idx.append(list_idx)
    
    return list_list_idx

def get_annotation(df_subj_timeline, list_list_idx):
    list_list_annotation = list()
    list_columns = list(df_subj_timeline.columns)
    for list_idx in list_list_idx:
        list_annot = list()
        for idx in list_idx:
            annot = list_columns[idx]
            list_annot.append(annot)
        list_list_annotation.append(list_annot)
    
    return list_list_annotation

def plot_timeline(list_timeline, list_annot, outfigname):
    from matplotlib.lines import Line2D
    legend_elements = [
                    #    Line2D([0], [0], color='orange',label='Enrolled~StartOxygen'), 
                    #    Line2D([0], [0], color='red',label='Start~EndOxygen'), 
                    #    Line2D([0], [0], color='green',label='EndOxygen~Released'), 
                       Line2D([0], [0], marker='^', color='k', markerfacecolor='purple', label='BloodDrawn', markersize=10),
                       Line2D([0], [0], marker='*', color='k', markerfacecolor='white', label='Admission/Release', markersize=10),]

    fig, axes = plt.subplots(nrows=len(list_timeline), ncols=1, figsize=(10, 20), constrained_layout = True, sharex=True, sharey=True)
    plt.subplots_adjust(wspace=0, hspace=0)
    ax = axes.flatten()

    for ind, (timeline, annot) in enumerate(zip(list_timeline, list_annot)):
        subj_id = str(timeline[0])
        severity = str(timeline[1])
        timepoints = timeline[2:]
        start = timepoints[0]
        timeperiods = list(map(lambda x: (x - start).days, timepoints))
        timeperiods_blooddrawn = timeperiods[5:]
        ax[ind].scatter(timeperiods[1], 0, s=150, c="white", edgecolors="k", marker= "*", zorder=1)
        ax[ind].scatter(timeperiods[4], 0, s=150, c="white", edgecolors="k", marker= "*", zorder=1)
        ax[ind].plot(timeperiods_blooddrawn, np.zeros_like(timeperiods_blooddrawn), "^", color="k", markerfacecolor="purple", markersize=15, linewidth=1, alpha=1, zorder=0)
        ax[ind].plot(timeperiods, np.zeros_like(timepoints), "-o", color="k", markerfacecolor="w", zorder=0)
        ax[ind].set_xticks(np.arange(0, 32, 1))
        ax[ind].set_xlabel("Infection Period (Days)", fontsize=16)
        plt.setp(ax[ind].xaxis.get_major_ticks(), visible=False)
        plt.setp(ax[ind].get_yticklabels(), visible=False)
        plt.setp(ax[ind].yaxis.get_major_ticks(), visible=False)
        if severity == "Mild":
            plt.setp(ax[ind].set_ylabel(subj_id, labelpad=40), rotation=360, ha="left", va="center", color="forestgreen", weight="bold", fontsize=12)
        else:
            plt.setp(ax[ind].set_ylabel(subj_id, labelpad=40), rotation=360, ha="left", va="center", color="firebrick", weight="bold", fontsize=12)

        list_xticks = list()
        for x in ax[ind].get_xticks():
            if x%7 == 0:
                list_xticks.append(x)
        [ax[ind].axvline(x, color='k', lw=1, linestyle="dashed", zorder=-1) for x in list_xticks]
        ax[ind].spines[["left", "top", "right", "bottom"]].set_visible(False)

    ax[ind].spines[["bottom"]].set_visible(True)
    plt.legend(handles= legend_elements, bbox_to_anchor=(1.05, 1.05))
    plt.setp(ax[ind].get_xticklabels(), visible=True)
    plt.setp(ax[ind].xaxis.get_major_ticks(), visible=True)
    fig.supylabel("Patients", fontsize=16)
    plt.tight_layout()
    plt.savefig(outfigname, dpi=600, bbox_inches='tight')
    plt.show()
    plt.close()


def read_crf_file(path_excel, list_drop_samples):
    df_crf = pd.read_excel(path_excel, engine="openpyxl", sheet_name="21-22등록 대상자_modify", skiprows=1)
    df_crf = correct_header_dataframe(df_crf)
    df_crf = select_dataframe(df_crf, num=2)
    df_crf = filter_dataframe(df_crf, list_drop_samples)

    return df_crf

def correct_header_dataframe(df_crf_raw):
    list_columns = list(df_crf_raw.columns)
    list_new_columns = list()
    for colname in list_columns:
        new_colname = colname.replace("\n", "_").replace(" ", "_").replace("\t", "_").replace("3-1. ", "").replace("3-2. ", "").replace("3-3. ", "").rstrip("_")
        list_new_columns.append(new_colname)
    
    df_crf_raw.columns = list_new_columns

    return df_crf_raw

def filter_dataframe(df_crf, list_drop_samples):
    if len(list_drop_samples) > 0:
        df_crf_filtered = df_crf[~df_crf.iloc[:, 0].isin(list_drop_samples)]
    else:
        df_crf_filtered = df_crf
        
    return df_crf_filtered

def select_dataframe(df_crf, num=2):
    df_crf_select = df_crf.iloc[:, num:]

    return df_crf_select


# %%
if __name__ == "__main__":
    WORKDIR = str(Path(__file__).parents[3])
    path_ehr = f"{WORKDIR}/Resources/Data/EHR/infectomics_CRF_UNIST_20250227.xlsx"
    path_drop_samples = f"{WORKDIR}/Resources/Scripts/Final/list_remove_samples.txt"
    path_rename = f"{WORKDIR}/Resources/Scripts/Final/list_variable_name_change.json"
    outfigname = f"_{WORKDIR}/Results/Paper/timeline_covid19.png"
    os.makedirs(os.path.dirname(outfigname), exist_ok=True)
    main(path_ehr, path_drop_samples, path_rename, outfigname)

# %%


# %%
import math
import warnings
from collections import defaultdict
from datetime import datetime

import matplotlib.dates as mdat
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")


# %% [clinicals]
def main(path_excel, list_drop_samples, dict_rename, outfigname):
    df_crf = read_crf_file(path_excel, list_drop_samples)
    list_columns_needed = [x for x in df_crf.columns if x in dict_rename.keys()]
    df_crf_columns_needed = df_crf[list_columns_needed]
    df_crf_columns_renamed = df_crf_columns_needed.rename(columns=dict_rename)
    df_crf_columns_sorted = df_crf_columns_renamed[list(dict_rename.values())]
    df_crf_columns_target = df_crf_columns_sorted[~df_crf_columns_sorted["Sample_ID"].str.contains("C19-R")]
    df_crf_columns_target = df_crf_columns_target[~df_crf_columns_target["Sample_ID"].str.contains("-L1")]
    df_crf_columns_target["Subject_ID"] = df_crf_columns_target["Sample_ID"].apply(lambda x: "-".join(x.split("-")[:-1]))
    df_crf_columns_target["Enrolled"] = df_crf_columns_target["Enrolled"].apply(lambda x: x.replace("\n", "") if str(x).startswith("\n") else x)
    df_crf_columns_target_fillna = fillna_by_first_row(df_crf_columns_target)
    df_crf_columns_target_fillna = fillna_no_values_by_test_date(df_crf_columns_target_fillna)
    df_crf_columns_target_fillna = df_crf_columns_target_fillna.applymap(convert_timestamp_to_datetime)
    df_crf_columns_target_fillna["Visit_Number"] = df_crf_columns_target_fillna["Sample_ID"].apply(lambda x: x.split("-")[-1])
    df_subj_timeline_final = collapse_by_subjid(df_crf_columns_target_fillna)
    list_raw_timeline = get_raw_timeline(df_subj_timeline_final)
    list_raw_timeline = sorted(list_raw_timeline, key=lambda x: x[1])
    list_timeline = get_nafilt_timeline(list_raw_timeline)
    list_list_ind = get_list_indices_no_na(list_raw_timeline)
    list_annot = get_annotation(df_subj_timeline_final, list_list_ind)
    
    plot_timeline(list_timeline, list_annot, outfigname)


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
    if type(x) != str and type(x) != int:
        timestamp = str(x)
        convert_x = datetime.strptime(str(timestamp), "%Y-%m-%d %H:%M:%S").date()

    else:
        convert_x = x
    
    return convert_x

def collapse_by_subjid(df_crf_columns_target_fillna):
    dict_subj_timeline = defaultdict(dict)
    prev = list()
    for ind, row in df_crf_columns_target_fillna.iterrows():
        dict_row = row.to_dict()
        id = (dict_row["Subject_ID"])
        visit = dict_row["Visit_Number"]
        blood = dict_row["BloodDrawn"]

        if id not in prev:
            dict_subj_timeline[id].update(dict_row)
            prev.append(id)

        dict_subj_timeline[id].update({visit: blood})
            
        prev = list()

    df_subj_timeline = pd.DataFrame.from_dict(dict_subj_timeline, orient="index")
    df_subj_timeline = df_subj_timeline.drop(columns=["BloodDrawn", "Visit_Number", "Sample_ID", "Subject_ID"])
    df_subj_timeline_final = df_subj_timeline.reset_index(drop=False).rename(columns={"index": "Subject_ID"})

    return df_subj_timeline_final

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
        list_annot_date_only = list_annot[2:]
        list_list_annotation.append(list_annot_date_only)
    
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
        severity = int(timeline[1])
        timepoints = timeline[2:]
        start = timepoints[0]
        timeperiods = list(map(lambda x: (x - start).days, timepoints))
        timeperiods_enrol_treat = timeperiods[1: 3]
        # timeperiods_treat_end = timeperiods[2: 4]
        timeperiods_end_release = timeperiods[3: 5]
        timeperiods_blooddrawn = timeperiods[5:]


        # levels = np.tile([-5, 5, -3, 3, -1, 1],
        #          int(np.ceil(len(timeperiods)/6)))[:len(timeperiods)]
        # ax[ind].vlines(timeperiods, 0, levels, color="tab:red")
        ax[ind].scatter(timeperiods[1], 0, s=150, c="white", edgecolors="k", marker= "*", zorder=1)
        ax[ind].scatter(timeperiods[4], 0, s=150, c="white", edgecolors="k", marker= "*", zorder=1)
        # if timeperiods_treat_end == [0, 0]:
        #     ax[ind].plot(timeperiods_enrol_treat, np.zeros_like(timeperiods_enrol_treat), "-o", color="orange", markerfacecolor="orange", markersize=10, linewidth=5, alpha=0, zorder=-1)
        # else:
        #     ax[ind].plot(timeperiods_enrol_treat, np.zeros_like(timeperiods_enrol_treat), "-o", color="orange", markerfacecolor="orange", markersize=10, linewidth=5, alpha=0.5, zorder=-1)
        # if timeperiods_treat_end == [0, 0]:
        #     ax[ind].plot(timeperiods_treat_end, np.zeros_like(timeperiods_treat_end), "-o", color="red", markerfacecolor="red", markersize=10, linewidth=5, alpha=0, zorder=-1)
        # else:
        #     ax[ind].plot(timeperiods_treat_end, np.zeros_like(timeperiods_treat_end), "-o", color="red", markerfacecolor="red", markersize=10, linewidth=5, alpha=0.5, zorder=-1)
        # if timeperiods_treat_end == [0, 0]:
        #     ax[ind].plot(timeperiods_end_release, np.zeros_like(timeperiods_end_release), "-o", color="green", markerfacecolor="green", markersize=10, linewidth=5, alpha=0, zorder=-1)
        # else:
        #     ax[ind].plot(timeperiods_end_release, np.zeros_like(timeperiods_end_release), "-o", color="green", markerfacecolor="green", markersize=10, linewidth=5, alpha=0.5, zorder=-1)
        ax[ind].plot(timeperiods_blooddrawn, np.zeros_like(timeperiods_blooddrawn), "^", color="k", markerfacecolor="purple", markersize=15, linewidth=1, alpha=1, zorder=0)
        ax[ind].plot(timeperiods, np.zeros_like(timepoints), "-o", color="k", markerfacecolor="w", zorder=0)
        # for d, l, r in zip(timeperiods, levels, annot):
        #     ax[ind].annotate(r, xy=(d, l),
        #                 xytext=(0, np.sign(l)*3), textcoords= "offset points",
        #                 horizontalalignment="left",
        #                 verticalalignment="bottom" if l > 0 else "top")
        ax[ind].set_xticks(np.arange(0, 32, 1))
        ax[ind].set_xlabel("Infection Period (Days)", fontsize=16)
        plt.setp(ax[ind].xaxis.get_major_ticks(), visible=False)
        plt.setp(ax[ind].get_yticklabels(), visible=False)
        plt.setp(ax[ind].yaxis.get_major_ticks(), visible=False)
        if severity == 1 or severity == 2:
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
    mode = "Methyl"
    path_excel = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/infectomics_CRF_20230410_edit.xlsx"
    list_drop_samples = ["C19-C045-V2",
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
    dict_rename = {"Subject_NO.(고유번호)":"Sample_ID", "중증도분류":"Severity", "진단일자":"TestedPositive", "본원_입원일자":"Enrolled", "Visit_date":"BloodDrawn", "산소시작일자":"StartOxygenTreatment", "산소중단일자":"EndOxgenTreatment","퇴원일자":"Released"}
    outfigname = f"/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/InfectomicsPaper1/timeline_covid19_20240621.png"
    main(path_excel, list_drop_samples, dict_rename, outfigname)

# %%

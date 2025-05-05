# %%
import json
import math
import os
import re
from collections import Counter, defaultdict
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ranksums, t
from statannot import add_stat_annotation
from statsmodels.stats.proportion import proportions_ztest


# %%
class Parse():

    def __init__(self, path_crf, list_drop_samples, path_rename):
        self.path_crf = path_crf
        self.list_drop_samples = list_drop_samples
        self.path_rename = path_rename
        self.df_crf = self.read_excel(engine="openpyxl", 
                                      sheet_name="21-22등록 대상자_modify", 
                                      skiprows=1)
        self.df_crf = self.correct_header_dataframe()
        self.df_crf = self.rename_column_names_dataframe()
        self.df_crf = self.select_dataframe()
        self.df_crf = self.filter_dataframe()
        self.df_crf = self.add_info_severity_dataframe()
        self.df_crf = self.extract_important_columns_dataframe()
        
    def read_excel(self, *args, **kwargs):
        df_excel = pd.read_excel(self.path_crf, *args, **kwargs)

        return df_excel

    def correct_header_dataframe(self):
        list_columns = list(self.df_crf.columns)
        list_new_columns = list(map(lambda x: x.replace("\n", "_"), list_columns))
        list_new_columns = list(map(lambda x: x.replace("\t", "_"), list_new_columns))
        list_new_columns = list(map(lambda x: x.replace("\xa0", "_"), list_new_columns))
        list_new_columns = list(map(lambda x: x.replace(".1", "_Quantitative"), list_new_columns))
        list_new_columns = list(map(lambda x: x.rstrip(), list_new_columns))
        list_new_columns = list(map(lambda x: x.replace(" ", "_"), list_new_columns))
        self.df_crf.columns = list_new_columns

        return self.df_crf

    def rename_column_names(self):
        with open(self.path_rename, mode="rb") as fj:
            dict_rename = json.load(fj)
        
        list_colnames = list(self.df_crf.columns)
        
        list_new_colnames = list()
        for colname in list_colnames:
            if colname not in dict_rename.keys():
                new_colname = colname
                list_new_colnames.append(new_colname)
            else:
                new_colname = dict_rename[colname]
                list_new_colnames.append(new_colname)
        
        return list_new_colnames

    def rename_column_names_dataframe(self):
        list_new_colnames = self.rename_column_names()    
        self.df_crf.columns = list_new_colnames
        
        return self.df_crf
    
    def filter_dataframe(self, col_sample="Sample_ID"):
        if len(self.list_drop_samples) > 0:
            df_crf_filtered = self.df_crf[~self.df_crf.loc[:, col_sample].isin(self.list_drop_samples)]
        else:
            df_crf_filtered = self.df_crf
            
        return df_crf_filtered

    def select_dataframe(self, col_start="Sample_ID", col_excl="NA"):
        list_cols_init = list(self.df_crf.columns)
        ind_start = list_cols_init.index(col_start)
        list_cols_start = list_cols_init[ind_start:]
        list_cols_select = list(filter(lambda x: x != col_excl, list_cols_start))
        df_crf_select = self.df_crf.loc[:, list_cols_select]

        return df_crf_select

    def change_severity_class_to_group(self, sev_class):
        sev_class = int(sev_class)
        if sev_class == 1 or sev_class == 2:
            sev_group = "Mild"

        else:
            sev_group = "Severe"
        
        return sev_group

    def add_info_severity_dataframe(self, col_sample="Sample_ID", col_target="Severity_group"):
        self.df_crf = self.df_crf.dropna(subset=[col_target])
        self.df_crf[col_target] = self.df_crf[col_target].apply(self.change_severity_class_to_group)
        for sampleid in self.df_crf.loc[:, col_sample]:
            if str(sampleid).split("-")[1][0] == "R":
                df_crf_indexed = self.df_crf.set_index(col_sample)
                df_crf_indexed.loc[sampleid, col_target] = "Convalescent"
                df_crf_add_sev = df_crf_indexed.reset_index(drop=False)

        return df_crf_add_sev
    
    def extract_important_columns_dataframe(self):
        list_new_colnames = self.rename_column_names()
        list_english_cols = list(filter(lambda x: re.search("[a-zA-Z\s]+", x[0]) is not None, list_new_colnames)) 
        list_english_cols = list(filter(lambda x: re.search("[a-zA-Z\s]+", x[-1]) is not None, list_english_cols))  
        df_crf_ext = self.df_crf[list_english_cols]
        
        return df_crf_ext
    
    def give_parsed_dataframe(self):
        
        return self.df_crf
    
    def lookup_dict_sample_visit(self, df_crf, col_sample="Sample_ID"):
        dict_sample_phase = dict()
        list_sampleid = df_crf.loc[:, col_sample].to_list()
        for sampleid in list_sampleid:
            subjectid = "-".join(sampleid.split("-")[:2])
            phase = sampleid.split("-")[-1]
            if dict_sample_phase.get(subjectid) == None:
                dict_sample_phase[subjectid] = list()
            dict_sample_phase[subjectid].append(phase)

        return dict_sample_phase
    
    def get_clinical_assoc_severity(self, df_crf, col_target="Severity_group", col_sample="Sample_ID"):
        list_df_compare = list()
        list_columns = list(df_crf.columns)
        df_crf_selec = df_crf.copy()
        for colname in list_columns:
            list_crf_values = df_crf_selec[colname].to_list()
            list_crf_values_nona = list(filter(lambda x : not isinstance(x, float) or not math.isnan(x), list_crf_values))
            list_crf_values_nona = list(map(lambda x: x.replace("PY", "") if "PY" in str(x) else x, list_crf_values_nona))
            if list_crf_values_nona == list():
                continue

            is_numeric = isinstance(list_crf_values_nona[0], int) or isinstance(list_crf_values_nona[0], float)
            if is_numeric:
                try:
                    df_crf_selec[colname] = df_crf_selec[colname].apply(float)
                    df_compare = df_crf_selec[[col_target, col_sample, colname]].dropna()
                    list_df_compare.append(df_compare)
                
                except Exception as e:
                    print(e)
        
        return list_df_compare

class Timeline():
    def __init__(self, df_crf):
        self.df_crf = df_crf

    def fillna_by_first_row(self):
        df_crf_fillna = self.df_crf.copy()
        list_columns = list(self.df_crf.columns)
        list_annot = list_columns[2:-1]
        for annot in list_annot:
            df_crf_fillna[annot] = self.df_crf[annot].fillna(self.df_crf.groupby("Subject_ID")[annot].transform("first"))

        return df_crf_fillna

    def fillna_no_values_by_test_date(self, df_crf_fillna):
        df_crf_fillna["StartOxygenTreatment"] = df_crf_fillna["StartOxygenTreatment"].fillna(value=df_crf_fillna["TestedPositive"])
        df_crf_fillna["EndOxgenTreatment"] = df_crf_fillna["EndOxgenTreatment"].fillna(value=df_crf_fillna["TestedPositive"])
        
        return df_crf_fillna

    def convert_timestamp_to_datetime(x):
        if type(x) != str and type(x) != int:
            timestamp = str(x)
            convert_x = datetime.strptime(str(timestamp), "%Y-%m-%d %H:%M:%S").date()

        else:
            convert_x = x
        
        return convert_x

    def collapse_by_subjid(self, df_crf_fillna):
        dict_subj_timeline = defaultdict(dict)
        prev = list()
        for ind, row in df_crf_fillna.iterrows():
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

    def get_annotation(self, df_subj_timeline, list_list_idx):
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

    def tmp(self, df_crf):
        self.df_crf = df_crf
        df_crf["SubjectID"] = df_crf["SampleID"].apply(lambda x: "-".join(x.split("-")[:-1]))
        df_crf["Enrolled"] = df_crf["Visit date"].apply(lambda x: x.replace("\n", "") if str(x).startswith("\n") else x)
        df_crf_fillna = self.fillna_by_first_row(df_crf)
        df_crf_fillna = self.fillna_no_values_by_test_date(df_crf_fillna)
        df_crf_fillna = df_crf_fillna.applymap(self.convert_timestamp_to_datetime)
        df_crf_fillna["Visit_Number"] = df_crf_fillna["Sample_ID"].apply(lambda x: x.split("-")[-1])
        df_subj_timeline_final = self.collapse_by_subjid(df_crf_fillna)
        list_raw_timeline = self.get_raw_timeline(df_subj_timeline_final)
        list_raw_timeline = sorted(list_raw_timeline, key=lambda x: x[1])
        list_timeline = self.get_nafilt_timeline(list_raw_timeline)
        list_list_ind = self.get_list_indices_no_na(list_raw_timeline)
        list_annot = self.get_annotation(df_subj_timeline_final, list_list_ind)

class Stats():
    
    def __init__(self, df_compare, outdir):
        self.df_compare = df_compare 
        self.target = "Severity_group"
        self.outdir = outdir
    
    def get_series_values_per_severity(self):
        value = list(self.df_compare.columns)[-1]
        series_sev_val = self.df_compare.groupby(self.target)[value].apply(np.array)

        return series_sev_val
    
    def is_prop(self, list_control, list_case):
        set_case_ctrl = set(list_control).union(set(list_case))
        num_cnt = len(set_case_ctrl)
        if (int(num_cnt) == 2) and (set_case_ctrl=={1,2}):
            return True
        
        return False
    
    def get_counts_by_severity(self, series_sev_val):
        mild_yes = dict(Counter(series_sev_val.get("mild", []))).get(2.0, 0)
        mild_no = dict(Counter(series_sev_val.get("mild", []))).get(1.0, 0)
        sev_yes = dict(Counter(series_sev_val.get("severe", []))).get(2.0, 0)
        sev_no = dict(Counter(series_sev_val.get("severe", []))).get(1.0, 0)
        
        return mild_yes, mild_no, sev_yes, sev_no
    
    def calculate_diff_prop(self, series_sev_val): 
        mild_yes, mild_no, sev_yes, sev_no = self.get_counts_by_severity(series_sev_val)
        if (mild_yes + mild_no) > 0:
            mild_prop = mild_yes / (mild_yes + mild_no)
        else:
            mild_prop = 0
            
        if (sev_yes + sev_no) > 0:
            sev_prop = sev_yes / (sev_yes + sev_no)
        else:
            sev_prop = 0

        diff_prop = sev_prop - mild_prop
        
        return diff_prop
    
    def calculate_diff_median(self, list_control, list_case):
        median_control = np.median(list_control)
        median_case = np.median(list_case)
        diff_median = (median_case - median_control)
        
        return diff_median
    
    def calculate_diff_stats(self):
        series_sev_val = self.get_series_values_per_severity()
        mild_yes, mild_no, sev_yes, sev_no = self.get_counts_by_severity(series_sev_val)
        control = series_sev_val.index[0]
        list_control = series_sev_val[control]
        case = series_sev_val.index[1]
        list_case = series_sev_val[case]
        if self.is_prop(list_control, list_case):
            diff_prop = self.calculate_diff_prop(series_sev_val)
            count = np.array([mild_yes, mild_no])
            nobs = np.array([(mild_yes + mild_no), (sev_yes + sev_no)])
            stat, pval = proportions_ztest(count, nobs)
            dict_stat_res = {
                        "compar": list(series_sev_val.index),
                        "test": "2_samp_prop",
                        "diff": diff_prop, 
                        "stat": stat, 
                        "pval": pval}
        else:
            diff_median = self.calculate_diff_median(list_control, list_case)
            stat, pval = ranksums(list_case, list_control)
            dict_stat_res = {
                        "compar": list(series_sev_val.index), 
                        "test": "ranksum",
                        "diff": diff_median, 
                        "stat": stat, 
                        "pval": pval}
    
        return dict_stat_res

class Plot(): 
    
    def __init__(self, dict_stat_res, dict_palette, outdir, ax):
        self.dict_stat_res = dict_stat_res
        self.dict_palette = dict_palette
        self.outdir = outdir
        self.ax = ax
        self.target = "Severity_group"

    def get_sample_size_group(self, df_compare):
        list_target = df_compare[self.target].to_list()
        count_target = dict(sorted(dict(Counter(list_target)).items()))
        
        return count_target     

    def edit_xticklabels(self, count_target, dict_rename={"Mild": "Mild-Moderate", "Severe": "Severe-Critical"}):
        new_xticklabels = list()
        for _, count_pair in enumerate(count_target.items(), start = 1):
            key, cnt = count_pair
            shortkey = key[0].upper()+key[1:]
            finalkey = dict_rename.get(shortkey, shortkey)
            xticklabel = f"{finalkey}\n(N={cnt})"
            new_xticklabels.append(xticklabel)

        return new_xticklabels

    def draw_boxplot(self, 
                     df_compare, 
                     colname, 
                     order=["Mild", "Severe"], 
                     dict_unit={
                        'CT(PCR)_E': "a.u.", 
                        'CT(PCR)_R': "a.u.", 
                        'CT(PCR)_N': "a.u.",
                        'RBC':"×10⁶/mm³",
                        'WBC':"×10³/mm³",
                        'Neutrophil': "%",
                        'Lymphocytes': "%",
                        'Monocytes': "%",
                        'Basophils': "%",
                        'Eosinophils': "%",
                        'Glucose': "mg/dL",
                        'CRP': "mg/dL"}
                    ):
        flierprops = dict(marker='o', markerfacecolor='None', markersize=5,  markeredgecolor='black')
        palette = list(map(lambda x: self.dict_palette[x], order))

        g = sns.boxplot(data=df_compare, 
                         x=self.target, 
                         y=colname, 
                         palette=palette, 
                         order=order, 
                         flierprops=flierprops,
                         ax=self.ax,
                         width=0.3,
                         zorder=3)
        
        g.set_ylabel(f"{colname}\n({dict_unit.get(colname, None)})", fontsize=16)
        g.tick_params(axis='y',labelsize=14)

        box_pairs = [tuple(sorted(self.dict_stat_res["compar"], reverse=True))]
        pvalues = [self.dict_stat_res["pval"]]

        add_stat_annotation(g, 
                            data=df_compare, 
                            x=self.target, 
                            y=colname, 
                            order=order,
                            box_pairs=box_pairs,
                            perform_stat_test=False, 
                            pvalues=pvalues,
                            test=None, 
                            text_format='star',
                            fontsize=12,  
                            loc='inside', 
                            verbose=0)

        count_target = self.get_sample_size_group(df_compare)
        new_xticklabels = self.edit_xticklabels(count_target)
        g.set_xticks(ticks=list(range(0, len(new_xticklabels))))
        g.set_xticklabels(labels=new_xticklabels, rotation=45, ha="right", rotation_mode='anchor', fontsize=14)
        g.set_xlabel("", fontsize=0)
        g.spines[['right', 'top']].set_visible(False)
        return g, pvalues

    def draw_timeline(self, list_timeline, list_annot):
        from matplotlib.lines import Line2D
        legend_elements = [Line2D([0], [0], marker='^', color='k', markerfacecolor='purple', label='BloodDrawn', markersize=10),
                           Line2D([0], [0], marker='*', color='k', markerfacecolor='white', label='Admission/Release', markersize=10),]

        _, axes = plt.subplots(nrows=len(list_timeline), ncols=1, figsize=(10, 20), constrained_layout = True, sharex=True, sharey=True)
        plt.subplots_adjust(wspace=0, hspace=0)
        ax = axes.flatten()

        for ind, (timeline, annot) in enumerate(zip(list_timeline, list_annot)):
            subj_id = str(timeline[0])
            severity = int(timeline[1])
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
            ax[ind].xaxis.get_major_ticks()
            ax[ind].get_yticklabels()
            ax[ind].yaxis.get_major_ticks()
            if severity == 1 or severity == 2:
                ax[ind].set_ylabel(subj_id, labelpad=40, rotation=360, ha="left", va="center", color="forestgreen", weight="bold", fontsize=12)
            else:
                ax[ind].set_ylabel(subj_id, labelpad=40, rotation=360, ha="left", va="center", color="firebrick", weight="bold", fontsize=12)

            list_xticks = list()
            for x in ax[ind].get_xticks():
                if x%7 == 0:
                    list_xticks.append(x)
            [ax[ind].axvline(x, color='k', lw=1, linestyle="dashed", zorder=-1) for x in list_xticks]
            ax[ind].spines[["left", "top", "right", "bottom"]].set_visible(False)

        ax[ind].spines[["bottom"]].set_visible(True)
        ax.legend(handles= legend_elements, bbox_to_anchor=(1.05, 1.05))
        ax[ind].get_xticklabels()
        ax[ind].xaxis.get_major_ticks()

        return ax
    

class Correlate():
    
    def __init__(self, omics, etc):
        self.omics = omics
        self.etc = etc
# %%


# %%
import glob
import os
import warnings

import pandas as pd

warnings.filterwarnings("ignore")

# %%
def read_excel(*args, **kwargs):
    df_excel = pd.read_excel(*args, **kwargs)

    return df_excel

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

def read_crf_file(path_excel, list_drop_samples):
    df_crf = pd.read_excel(path_excel, engine="openpyxl", sheet_name="21-22등록 대상자_modify", skiprows=1)
    df_crf = correct_header_dataframe(df_crf)
    df_crf = select_dataframe(df_crf, num=2)
    df_crf = filter_dataframe(df_crf, list_drop_samples)

    return df_crf

def get_list_files_methylcpgmin(dir_methylcpgmin):
    list_files_methylcpgmin = glob.glob(f"{dir_methylcpgmin}/**/*pair_merged.methyl_cpg_min.tsv", recursive=True)

    return list_files_methylcpgmin

def get_list_name_sample(list_files_methylcpgmin):
    list_name_sample = list()
    for file_methylcpgmin in list_files_methylcpgmin:
        dir_methylcpgmin = os.path.basename(os.path.dirname(file_methylcpgmin))
        if dir_methylcpgmin == "HealthyControl":
            name_sample = os.path.basename(file_methylcpgmin).split(".")[0]
        else:
            name_sample = os.path.basename(dir_methylcpgmin)
        list_name_sample.append(name_sample)

    return list_name_sample

# %%
mode = "Methyl"
path_excel = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/Vaccination/Resources/Data/infectomics_CRF_20230410_edit.xlsx" 
dir_methylcpgmin = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Resources/Infectomics_COVID-19_RNA/Backup/Copy_from_Shrimp/COVID19Infected/Results/10_methyl/MethylCpGMin"
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
                    'C19-C061-V3',
                    'C19-C002-L1',
                    'C19-C010-L1',
                    'C19-C017-L1',
                    'C19-C022-L1',
                    'C19-R001-L1',
                    'C19-R096-L1',
                    'C19-R092-L1',
                    'C19-R081-L1',
                    'C19-R078-L1',
                    'C19-R028-L1',
                    'C19-R024-L1',
                    'C19-R010-L1',
                    'C19-R091-L1',
                    'C19-R030-L1']
outfilename = f"/BiO/Access/kyungwhan1998/Infectomics/Results/9_clinical/Infectomics_Severity_Information_20240926.tsv"

# %%
df_crf = read_crf_file(path_excel, list_drop_samples)
dict_rename = {"Subject_NO.(고유번호)":"Sample_ID", "성별":"Sample_Sex", "만나이":"Sample_Age", "중증도분류":"Severity"}
list_columns_needed = [x for x in df_crf.columns if x in dict_rename.keys()]
df_crf_columns_needed = df_crf[list_columns_needed]
df_sev_info = df_crf_columns_needed.rename(columns=dict_rename)

def get_visit_order(sample_id):
    visit_order = sample_id.split("-")[-1]
    visit_order = visit_order.replace("V", "Visit")
    
    return visit_order

def get_sex_binary(sex_text):
    if str(sex_text).startswith("남"):
        sex_binary = 1
    else:
        sex_binary = 2
    
    return sex_binary

def is_conval_sample(sample_id):
    flag_conval = str(sample_id).split("-")[1][0]
    if flag_conval == "R":
        return True

    return False

# %%
df_sev_info["Visit_order"] = df_sev_info["Sample_ID"].apply(get_visit_order)
df_sev_info["Sample_Sex"] = df_sev_info["Sample_Sex"].apply(get_sex_binary)
df_sev_info["is_conval"] = df_sev_info["Sample_ID"].apply(is_conval_sample)
df_sev_info["Severity_visit"] = df_sev_info["Severity"].astype(str) + "_" + df_sev_info["Visit_order"].astype(str)
filt_conval = (df_sev_info["is_conval"])
df_sev_info.loc[filt_conval, "Visit_order"] = "Visit5"
df_sev_info.loc[filt_conval, "Severity_visit"] = "Convalescent"
# %%
list_files_methylcpgmin = get_list_files_methylcpgmin(dir_methylcpgmin)
list_name_sample = get_list_name_sample(list_files_methylcpgmin)
list_name_subject = list(map(lambda x: "-".join(x.split("-")[:2]), list_name_sample))
list_healthy = list(filter(lambda x: not str(x).startswith("C19"), list_name_subject))

# %%
df_healthy = pd.DataFrame(columns=df_sev_info.columns)
df_healthy["Sample_ID"] = list_healthy
df_healthy["Severity"] = "0"
df_healthy["Visit_order"] = "Visit0"
df_healthy["Severity_visit"] = "Healthy"
df_sev_info_added_healthy = pd.concat([df_sev_info, df_healthy], axis=0)

df_sev_info_final = df_sev_info_added_healthy.sort_values(by=["Sample_ID"], ascending=True)
df_sev_info_final.to_csv(outfilename, sep="\t", index=False)

# %%

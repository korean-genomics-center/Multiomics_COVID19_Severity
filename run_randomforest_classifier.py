#%%
import os
import pickle
import random
from copy import deepcopy

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold

PARAM_TYPE = {
    "max_depth" : int,
    "min_samples_split" : int,
    "min_samples_leaf" : int,
    "n_estimators" : int
}

# %%    
def make_save_dirs(path_save_model, path_save_test_result):
    os.makedirs(os.path.dirname(path_save_model), exist_ok=True)
    os.makedirs(os.path.dirname(path_save_test_result), exist_ok=True)
    
def train_final_randomforest_model(train_x, train_y, n_jobs, random_state):
    trained_model = _train_randomforest_model(train_x, train_y, n_jobs, random_state)
    
    return trained_model

def _train_randomforest_model(train_x, train_y, n_jobs, random_state):
    randomforest_params = set_default_randomforest_param(n_jobs, random_state)
    rfmodel = get_randomforest_model_with_param(**randomforest_params)
    rfmodel.fit(train_x, train_y)
    return rfmodel  
    
def set_default_randomforest_param(n_jobs, random_state, **params):
    cv_param = {
        "n_jobs" : n_jobs, 
        "random_state" : random_state
    }
    cv_param.update(params)
    return cv_param

def set_grid_tuning_param(**parameter_tuning):
    dict_param = {
        "max_depth" : np.linspace(5, 100, 6), # step=13
        "min_samples_split" : np.linspace(2, 10, 3), # step=4
        "min_samples_leaf" : np.linspace(1, 5, 3), # step=2
        "n_estimators" : np.linspace(100, 1000, 4) # step=300
    }
    for param, value in parameter_tuning.items():
        if len(value) == 1:
            dict_param[param] = [value]
        else:
            start, end, num = value
            dict_param[param] = np.linspace(start, end, int(num))
            
    for param, datrange in dict_param.items():
        list_datrange = make_np_array_to_list(datrange, PARAM_TYPE[param])
        dict_param[param] = list_datrange
    return dict_param

def make_np_array_to_list(nparray, dattype):
    list_value = list(map(dattype, nparray))
    return list_value

def get_randomforest_model_with_param(**params):
    return RandomForestClassifier(**params)

def test_randomforest_model(rfmodel, test_sample_id, test_x, test_y):
    pred_y = rfmodel.predict(test_x)
    pred_y_proba = rfmodel.predict_proba(test_x)
    
    accuracy = accuracy_score(test_y, pred_y)
    print(f"Final accuracy : {accuracy}", flush = True)
    table_proba = pd.DataFrame(pred_y_proba)
    table_proba["Sample_ID"] = test_sample_id
    table_proba["Answer"] = test_y
    table_proba["Prediction"] = pred_y
    table_proba["Correct"] = table_proba["Answer"] == table_proba["Prediction"]    
    table_proba = table_proba.rename(columns= {0 : "Proba_Control", 1 : "Proba_Case"})
    table_proba["Answer"] = table_proba["Answer"].apply(lambda x : "Case" if x == 1 else "Control")
    table_proba["Prediction"] = table_proba["Prediction"].apply(lambda x : "Case" if x == 1 else "Control")
    return table_proba

def save_obj_pickle(obj, path_save):
    with open(path_save, 'wb') as fw:
        pickle.dump(obj, fw)

def load_obj_pickle(path_save):
    with open(path_save, 'rb') as fr:
        return pickle.load(fr)

#%% 
WORKDIR = "/BiO/Access/kyungwhan1998/Infectomics/Results/InfectomicsPaper1/20240906"
path_table = f"{WORKDIR}/Methylation_Beta_32_Markers_COVID19.txt"
path_meta = f"{WORKDIR}/Metadata_Train_Test_Split_COVID19.txt"
col_meta_pred = "Severity_group"
col_meta_sampleid = "Sample_ID"
col_meta_dataset = "Dataset"
n_jobs = 20

cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)

table_values = pd.read_csv(path_table, sep = '\t', index_col=[0])
table_meta = pd.read_csv(path_meta, sep = '\t', index_col=[0])

table_X = table_values.T.values
table_y = np.array(list(map(lambda x: 1 if x == "Severe" else 0, table_meta[col_meta_pred].to_list())))

auroc_scores = []

for i, (train_index, test_index) in enumerate(cv.split(table_X, table_y), start=1):
    train_X, test_X = table_X[train_index], table_X[test_index]
    train_y, test_y = table_y[train_index], table_y[test_index]
    test_sample_id = list(np.array(table_values.columns)[test_index])
    rf_model = train_final_randomforest_model(train_X, train_y, n_jobs, random_state=42)
    path_save_model = f"{WORKDIR}/RF_Classifier_Model_Severity_COVID19_3Fold_CV_Fold{i}.pk"
    save_obj_pickle(rf_model, path_save_model)
    table_proba = test_randomforest_model(rf_model, test_sample_id, test_X, test_y)
    path_save_test_result = f"{WORKDIR}/RF_Classification_Result_COVID19_3Fold_CV_Fold{i}.txt"
    table_proba.to_csv(path_save_test_result, sep="\t")
    

# %%

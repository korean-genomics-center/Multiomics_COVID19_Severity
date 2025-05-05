# %%
import os
import pickle
import random
from copy import deepcopy

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold

# %%
    
def make_save_dirs(path_save_model, path_save_test_result):
    os.makedirs(os.path.dirname(path_save_model), exist_ok=True)
    os.makedirs(os.path.dirname(path_save_test_result), exist_ok=True)
    
def get_logisticclassifier_model_with_param(**params):
    
    return LogisticRegression(solver="liblinear", **params)

def _train_logisitcclassifer_model(train_X, train_y, n_jobs, random_state):
    parameters = dict()
    parameters.update({"n_jobs": n_jobs, "random_state": random_state})
    logimodel = get_logisticclassifier_model_with_param(**parameters)
    logimodel.fit(train_X, train_y)
    
    return logimodel

def train_final_logisticClassifier_model(train_X, train_y, n_jobs, random_state):
    list_trained_model = _train_logisitcclassifer_model(train_X, train_y, n_jobs, random_state)
    
    return list_trained_model

def test_logi_model(logimodel, test_sample_id, test_x, test_y):
    pred_y = logimodel.predict(test_x)
    pred_y_proba = logimodel.predict_proba(test_x)

    accuracy = accuracy_score(test_y, pred_y)
    print(accuracy)
    table_proba = pd.DataFrame(pred_y_proba)
    table_proba["Sample_ID"] = test_sample_id
    table_proba["Answer"] = test_y
    table_proba["Prediction"] = pred_y
    table_proba["Correct"] = (table_proba["Answer"] == table_proba["Prediction"])
    table_proba = table_proba.rename(columns = {0: "Proba_Control", 1: "Proba_Case"})
    
    return table_proba

def save_obj_pickle(obj, path_save):
    with open(path_save, "wb") as fw:
        pickle.dump(obj, fw)
        
# %%
WORKDIR = "/BiO/Access/kyungwhan1998/Infectomics/Results/InfectomicsPaper1/20240906"
path_table = f"{WORKDIR}/Methylation_Beta_32_Markers_COVID19.txt"
path_meta = f"{WORKDIR}/Metadata_Train_Test_Split_COVID19.txt"
path_save_model = f"{WORKDIR}/Logisitic_Classifier_Severity_COVID19.pk"
path_save_test_result = f"{WORKDIR}/Logisitic_Classification_Result_COVID19.txt"
col_meta_pred = "Severity_group"
col_meta_sampleid = "Sample_ID"
col_meta_dataset = "Dataset"
n_jobs = 3

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
    logi_model = train_final_logisticClassifier_model(train_X, train_y, n_jobs, random_state=42)
    path_save_model = f"{WORKDIR}/Logistic_Classifier_Model_Severity_COVID19_3Fold_CV_Fold{i}.pk"
    save_obj_pickle(logi_model, path_save_model)
    table_proba = test_logi_model(logi_model, test_sample_id, test_X, test_y)
    path_save_test_result = f"{WORKDIR}/Logistic_Classification_Result_COVID19_3Fold_CV_Fold{i}.txt"
    table_proba.to_csv(path_save_test_result, sep="\t")
    

    
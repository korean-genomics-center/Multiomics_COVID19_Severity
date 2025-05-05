# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (auc, classification_report, confusion_matrix,
                             roc_curve)


def get_metrics(dir_res, classifier, fold, threshold=0.5):
    path_res = f"{dir_res}/{classifier}_Classification_Result_COVID19_3Fold_CV_Fold{fold}.txt"
    df_res = pd.read_csv(path_res, sep="\t")
    y_truth = df_res["Answer"].to_list()
    y_test = y_truth
    # y_test = list(map(lambda x: 1 if x == "Case" or x == "Severe" else 0, y_truth))
    y_scores = df_res["Proba_Case"].to_list()
    
    # Compute ROC curve
    fpr, tpr, thresholds = roc_curve(y_test, y_scores)
    roc_auc = auc(fpr, tpr)
    
    # Compute predictions based on threshold
    y_pred = [1 if prob >= threshold else 0 for prob in y_scores]
    
    # Confusion matrix
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    
    # Sensitivity (Recall/TPR) and Specificity (TNR)
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    
    return fpr, tpr, roc_auc, thresholds, sensitivity, specificity


# %%
# Path and classifier setup
dir_res = "/BiO/Access/kyungwhan1998/Infectomics/Results/InfectomicsPaper1/20240906"
list_classifier = ["XGBoost", "RF", "Logistic"]
classifier = list_classifier[2]

# Prepare to store results
mean_fpr = np.linspace(0, 1, 100)
tpr_interpolated = []
aucs = []
sensitivities = []
specificities = []

plt.figure(figsize=(7, 7))

# Process each fold
for fold in range(1, 4):
    fpr, tpr, roc_auc, thresholds, sensitivity, specificity = get_metrics(dir_res, classifier, fold)
    plt.plot(fpr, tpr, lw=2, label=f'{classifier} Fold {fold} (AUC = {roc_auc:.2f})', alpha=0.7)
    
    aucs.append(roc_auc)
    sensitivities.append(sensitivity)
    specificities.append(specificity)
    
    # Interpolate TPR for mean ROC
    interp_tpr = np.interp(mean_fpr, fpr, tpr)
    interp_tpr[0] = 0.0
    tpr_interpolated.append(interp_tpr)

# Compute mean TPR, sensitivity, and specificity
mean_tpr = np.mean(tpr_interpolated, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_tpr = np.std(tpr_interpolated, axis=0)

mean_sensitivity = np.mean(sensitivities)
mean_specificity = np.mean(specificities)

# Plot mean ROC curve
plt.plot(
    mean_fpr,
    mean_tpr,
    color='tab:red',
    label=f'Mean ROC (AUC = {mean_auc:.2f})',
    lw=3,
    alpha=0.8
)

# Plot standard deviation as shaded area
plt.fill_between(
    mean_fpr,
    np.maximum(mean_tpr - std_tpr, 0),
    np.minimum(mean_tpr + std_tpr, 1),
    color='tab:red',
    alpha=0.2,
    label='± 1 std. dev.'
)

# Random guessing line
plt.plot([0, 1], [0, 1], color='grey', linestyle='--', lw=1, label='Random')

# Determine Youden index (J)
youden_index = mean_tpr - mean_fpr
ind_optimal = np.argmax(youden_index)
optimal_threshold = mean_fpr[ind_optimal]
tpr_optimal = mean_tpr[ind_optimal]
fpr_optimal = mean_fpr[ind_optimal]
sensitivity_optimal = tpr_optimal
specificity_optimal = 1 - fpr_optimal

# Plot Youden index point
plt.scatter([fpr_optimal], [tpr_optimal], marker="x", color='k', facecolor=None, linewidth=3, s=100, zorder=5, label='Youden Index Optimal Point')

# Add arrows for sensitivity and specificity
plt.annotate(
    f'Sensitivity: {sensitivity_optimal:.2f}\nSpecificty: {specificity_optimal:.2f}',
    xy=(1-specificity_optimal+0.01, sensitivity_optimal-0.01),
    xytext=(1-specificity_optimal+0.1, sensitivity_optimal-0.1),
    arrowprops=dict(facecolor='k', arrowstyle="->"),
    fontsize=10,
    ha="left",
    va="top",
    color='k'
)

# Final adjustments
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate (FPR)', fontsize=12)
plt.ylabel('True Positive Rate (TPR)', fontsize=12)
plt.title('Mean Receiver Operating Characteristic (ROC) Curve', fontsize=14)
plt.legend(loc="lower right")
plt.grid(alpha=0.3)
plt.show()

print(f"Optimal Sensitivity (Youden Index): {sensitivity_optimal:.2f}")
print(f"Optimal Specificity (Youden Index): {specificity_optimal:.2f}")
print(f"Optimal Threshold: {optimal_threshold}")



# %%
path_res = f"{dir_res}/{classifier}_Classification_Result_COVID19_3Fold_CV_Fold{fold}.txt"
df_res = pd.read_csv(path_res, sep="\t")
y_truth = df_res["Answer"].to_list()
y_test = y_truth
# y_test = list(map(lambda x: 1 if x == "Case" or x == "Severe" else 0, y_truth))
y_scores = df_res["Proba_Case"].to_list()
# Use the optimal threshold for classification
optimal_predictions = (y_scores >= optimal_threshold).astype(int)

y_real = list(map(lambda x: "Severe" if x == 1 else "Mild", y_test))
y_pred = list(map(lambda x: "Severe" if x == 1 else "Mild", optimal_predictions))
pred_res = list(y_test == optimal_predictions)
df_pred = pd.DataFrame()
df_pred["Real"] = y_real
df_pred["Pred"] = y_pred
df_pred["Outcome"] = pred_res 

# Evaluate the model with the optimal threshold
print("Confusion Matrix:")
print(confusion_matrix(y_test, optimal_predictions))
print("\nClassification Report:")
print(classification_report(y_test, optimal_predictions))
# %%

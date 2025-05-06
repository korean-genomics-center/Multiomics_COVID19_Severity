import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ranksums

file_path = "/BiO/Research/Project2/Infectomics_COVID-19_Host/Analysis/Infectomics_COVID-19_Methyl_Severity/Analysis/Figures/Resource/SupplemantaryFigure_5/infectomics_clinical_info_for_celltype_analy.xlsx"
df = pd.read_excel(file_path)

cell_types = ['Neutrophil', 'Lymphocytes', 'Monocytes', 'Eosinophils', 'Basophils']
id_col = "Subject NO.(고유번호)"

cols_needed = [id_col, '중증도분류', 'VISIT_no.'] + cell_types
df_filtered = df[cols_needed].copy()

def classify_severity(severity):
    if severity in [1, 2]:
        return 'Mild-Moderate'
    elif severity in [3, 4]:
        return 'Severe-Critical'
    else:
        return None

df_filtered['Severity Group'] = df_filtered['중증도분류'].apply(classify_severity)
df_filtered = df_filtered[
    df_filtered['Severity Group'].notna() &
    df_filtered['VISIT_no.'].isin([1, 2])
]

melted_df = df_filtered.melt(
    id_vars=[id_col, 'VISIT_no.', 'Severity Group'],
    value_vars=cell_types,
    var_name='Cell Type', value_name='Value'
)

palette = {
    'Mild-Moderate': 'forestgreen',
    'Severe-Critical': 'firebrick'
}
order = ['Mild-Moderate', 'Severe-Critical']

n_plots = len(cell_types) * 2
ncols = 4
nrows = (n_plots + ncols - 1) // ncols

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 10), sharex=False, sharey=False)
axes = axes.flatten()


for idx, (i, cell_type) in enumerate(enumerate(cell_types)):
    for j, visit in enumerate([1, 2]):
        plot_idx = idx * 2 + j
        ax = axes[plot_idx]
        data = melted_df[
            (melted_df['Cell Type'] == cell_type) &
            (melted_df['VISIT_no.'] == visit)
        ]

        # sns.boxplot(x='Severity Group', y='Value', data=data,
        #             palette=palette, order=order, ax=ax, width=0.35,
        sns.boxplot(x='Severity Group', y='Value', data=data,
                    palette=palette, order=order, ax=ax, zorder=1, width=0.35,
                    fliersize=3, linewidth=1.5,
                    flierprops=dict(marker='o', markerfacecolor='none', markeredgecolor='black', markersize=5, linestyle='none'))

        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        ax.set_xlabel('')
        if j == 0:
            ax.set_ylabel(f"{cell_type} (%)", fontsize=14, fontweight='bold')
        else:
            ax.set_ylabel('')
        ax.grid(True, axis='y', linestyle='--', linewidth=0.5, color='gray', alpha=0.7, zorder=0)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        for spine in ax.spines.values():
            spine.set_linewidth(1.0)
            spine.set_color("black")

        group1 = data[data['Severity Group'] == 'Mild-Moderate']['Value'].dropna()
        group2 = data[data['Severity Group'] == 'Severe-Critical']['Value'].dropna()
        if len(group1) > 0 and len(group2) > 0:
            stat, pval = ranksums(group1, group2)
            y_max = max(data['Value'].dropna())
            y_min = min(data['Value'].dropna())
            y, h = y_max + (y_max - y_min) * 0.1, (y_max - y_min) * 0.03

            ax.plot([0, 0, 1, 1], [y, y + h, y + h, y], lw=1.2, color='black')
            p_text = 'ns' if pval >= 0.05 else ("*" if pval < 0.05 else "")
            ax.text(0.5, y + h * 1.1, p_text, ha='center', va='bottom', fontsize=12)

# Hide any unused subplots
for k in range(n_plots, len(axes)):
    fig.delaxes(axes[k])

for col in range(ncols):
    if col % 2 == 0:
        fig.text(x=0.13 + (col * 0.26), y=0.94, s='Baseline', ha='center', fontsize=14, fontweight='bold')
    else:
        fig.text(x=0.13 + (col * 0.26), y=0.94, s='After 1 week', ha='center', fontsize=14, fontweight='bold')

plt.tight_layout(rect=[0.03, 0.03, 1, 0.93], w_pad=1)

#plt.tight_layout(rect=[0.03, 0.03, 1, 0.97], w_pad=1)
#plt.savefig("boxplot_MM_vs_SC_celltype_updated_1.png", dpi=300, bbox_inches='tight')
plt.savefig("SupplementaryFigure5.png", dpi=300, bbox_inches='tight')
plt.show()

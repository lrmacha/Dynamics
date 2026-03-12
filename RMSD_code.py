import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re
from scipy.stats import shapiro, kruskal, f_oneway
import scikit_posthocs as sp
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# ---------------------------
# Superscript function
# ---------------------------
def superscript_mutant(name):
    superscript_map = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    return re.sub(r'(\d+)', lambda x: x.group().translate(superscript_map), name)

# ---------------------------
# Test selection
# ---------------------------
def choose_test(df_data):
    normality = []
    for col in df_data.columns:
        stat, p = shapiro(df_data[col].dropna())
        normality.append(p > 0.05)
    all_normal = all(normality)
    
    if all_normal:
        groups = [df_data[col].dropna() for col in df_data.columns]
        stat, p = f_oneway(*groups)
        df_long = df_data.melt(var_name='Receptor', value_name='Value')
        tukey = pairwise_tukeyhsd(df_long['Value'], df_long['Receptor'])
        posthoc = pd.DataFrame(
            data=tukey._results_table.data[1:],
            columns=tukey._results_table.data[0]
        )
        test_type = 'ANOVA'
    else:
        groups = [df_data[col].dropna() for col in df_data.columns]
        stat, p = kruskal(*groups)
        df_long = df_data.melt(var_name='Receptor', value_name='Value')
        posthoc = sp.posthoc_dunn(df_long, val_col='Value', group_col='Receptor', p_adjust='bonferroni')
        test_type = 'Kruskal-Wallis'
    
    return stat, p, posthoc, test_type, all_normal

# ---------------------------
# Process RMSD CSV
# ---------------------------
def process_rmsd_csv(file_path, val_name='RMSD'):
    df = pd.read_csv(file_path)
    df_data = df.iloc[:, 1:]
    df_data.columns = df_data.columns.map(superscript_mutant)
    
    df_long = df_data.melt(var_name='Receptor', value_name=val_name)
    means = df_data.mean()
    stds = df_data.std()
    
    stat, p, posthoc, test_type, normal = choose_test(df_data)
    WT_label = superscript_mutant('WT')
    print(f"{file_path} -> {test_type}: statistic={stat:.3f}, p={p:.4f}, normal={normal}")
    
    return df_long, means, stds, posthoc, WT_label, test_type, normal, df_data

# ---------------------------
# p-value stars
# ---------------------------
def p_to_stars(pval):
    if pval < 0.001:
        return '***'
    elif pval < 0.01:
        return '**'
    elif pval < 0.05:
        return '*'
    else:
        return 'ns'

# ---------------------------
# Plot RMSD bar chart
# ---------------------------
def plot_rmsd(ax, df_long, means, stds, posthoc, WT_label, val_name='RMSD', panel_label='A', show_ylabel=True, normal=True):
    original_order = list(means.index)
    palette = ['#0072B2' if col == WT_label else '#E69F00' for col in original_order]
    
    sns.barplot(
        data=df_long,
        x='Receptor',
        y=val_name,
        order=original_order,
        palette=palette,
        ci='sd',
        edgecolor='black',
        linewidth=1.5,
        capsize=0.15,
        ax=ax
    )
    
    if show_ylabel:
        ax.set_ylabel(f'Average {val_name} (Å)')
    ax.set_xlabel('Substitution')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    # Remove gridlines
    ax.grid(False)
    ax.yaxis.grid(False)
    ax.xaxis.grid(False)
    
    for i, receptor in enumerate(original_order):
        if receptor != WT_label:
            if normal:
                row = posthoc[
                    ((posthoc['group1'] == WT_label) & (posthoc['group2'] == receptor)) |
                    ((posthoc['group1'] == receptor) & (posthoc['group2'] == WT_label))
                ]
                if not row.empty:
                    pval = float(row['p-adj'].iloc[0])
                else:
                    pval = 1.0
            else:
                pval = posthoc.loc[WT_label, receptor]
            
            bar_height = means[receptor]
            error = stds[receptor]
            y = bar_height + error + 0.005 * np.max(means)
            ax.text(i, y, p_to_stars(pval), ha='center', va='bottom', fontsize=12, color='black')
    
    ax.text(-0.1, 1.05, panel_label, transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='right')

# ---------------------------
# Process both RMSD files
# ---------------------------
df_long_L, means_L, stds_L, posthoc_L, WT_label_L, test_L, normal_L, df_data_L = process_rmsd_csv(
    "C:/Users/lrmacha/Downloads/RMSD_L.csv", val_name='RMSD_L'
)
df_long_R, means_R, stds_R, posthoc_R, WT_label_R, test_R, normal_R, df_data_R = process_rmsd_csv(
    "C:/Users/lrmacha/Downloads/RMSD_R.csv", val_name='RMSD_R'
)

# ---------------------------
# Create figure with 2 bar plots stacked vertically
# ---------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharey=False)

plot_rmsd(ax1, df_long_L, means_L, stds_L, posthoc_L, WT_label_L,
          val_name='RMSD_L', panel_label='A', show_ylabel=True, normal=normal_L)
plot_rmsd(ax2, df_long_R, means_R, stds_R, posthoc_R, WT_label_R,
          val_name='RMSD_R', panel_label='B', show_ylabel=True, normal=normal_R)

plt.tight_layout()

# ---------------------------
# Save figure only
# ---------------------------
output_pdf = "C:/Users/lrmacha/Downloads/RMSD_LR_panels.pdf"
output_svg = "C:/Users/lrmacha/Downloads/RMSD_LR_panels.svg"

fig.savefig(output_pdf, format='pdf', bbox_inches='tight')
fig.savefig(output_svg, format='svg', bbox_inches='tight')

print(f"Figure saved as PDF: {output_pdf}")
print(f"Figure saved as SVG: {output_svg}")

plt.show()
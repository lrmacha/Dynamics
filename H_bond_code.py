import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kruskal
import scikit_posthocs as sp
import re
import numpy as np

# ---------------------------
# Load CSV
# ---------------------------
file_path = "C:/Users/lrmacha/Downloads/H_Bonds.csv"
df = pd.read_csv(file_path)

# Drop time column
df_data = df.iloc[:, 1:]

# Convert to long format
df_long = df_data.melt(var_name='Receptor', value_name='H-Bonds')

# ---------------------------
# Superscript numbers in mutant names
# ---------------------------
def superscript_mutant(name):
    superscript_map = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    return re.sub(r'(\d+)', lambda x: x.group().translate(superscript_map), name)

df_long['Receptor'] = df_long['Receptor'].apply(superscript_mutant)

# Compute mean and SD
means = df_data.mean()
stds = df_data.std()
means.index = df_data.columns.map(superscript_mutant)
stds.index = df_data.columns.map(superscript_mutant)

# ---------------------------
# Color-blind-friendly palette
# WT = blue, mutants = orange
# ---------------------------
WT_label = superscript_mutant('WT')
palette = ['#0072B2' if col == 'WT' else '#E69F00' for col in df_data.columns]

# ---------------------------
# Plot figure
# ---------------------------
fig, ax = plt.subplots(figsize=(10, 5))

sns.barplot(
    data=df_long,
    x='Receptor',
    y='H-Bonds',
    palette=palette,
    ci='sd',
    edgecolor='black',
    linewidth=1.5,
    capsize=0.15,
    ax=ax
)

ax.set_ylabel('Average Number of H-Bonds')
ax.set_xlabel('Substitution')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
ax.grid(False)
ax.yaxis.grid(False)
ax.xaxis.grid(False)

# ---------------------------
# Kruskal-Wallis test
# ---------------------------
groups = [df_data[col].dropna() for col in df_data.columns]
stat, p = kruskal(*groups)
print(f"Kruskal-Wallis H-statistic = {stat:.3f}, p = {p:.4f}")

# Post-hoc Dunn test
df_posthoc = sp.posthoc_dunn(df_long, val_col='H-Bonds', group_col='Receptor', p_adjust='bonferroni')

# ---------------------------
# Add significance stars above bars
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

for i, receptor in enumerate(df_long['Receptor'].unique()):
    if receptor != WT_label:
        pval = df_posthoc.loc[WT_label, receptor]
        bar_height = means[receptor]
        error = stds[receptor]
        y = bar_height + error + 0.005 * np.max(means)
        ax.text(i, y, p_to_stars(pval), ha='center', va='bottom', fontsize=12, color='black')

plt.tight_layout()

# ---------------------------
# Save figure
# ---------------------------
output_pdf = "C:/Users/lrmacha/Downloads/H_Bonds_plot.pdf"
output_svg = "C:/Users/lrmacha/Downloads/H_Bonds_plot.svg"

fig.savefig(output_pdf, format='pdf', bbox_inches='tight')
fig.savefig(output_svg, format='svg', bbox_inches='tight')

print(f"Figure saved as PDF: {output_pdf}")
print(f"Figure saved as SVG: {output_svg}")

plt.show()
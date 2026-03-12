import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from scipy.stats import f_oneway

# ---------------------------
# Helper functions
# ---------------------------
superscript_map = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

def superscript_codons(label):
    """Format codons with superscript numbers."""
    if label == 'WT':
        return label
    parts = re.findall(r'[A-Za-z]+|\d+|/', label)
    return ''.join([
        p.upper() if p.isalpha() else 
        p.translate(superscript_map) if p.isdigit() else 
        p for p in parts
    ])

def format_p(p):
    """Format p-values in journal style."""
    if p < 0.0001:
        return "p < 0.0001"
    return f"p = {p:.4f}"

# ---------------------------
# Data
# ---------------------------
data = {
    'WT':          [-90.603, -76.675, -62.747],
    'H148P':       [-86.168, -85.2021, -84.236],
    'H148A':       [-84.0843, -66.3143, -48.5443],
    'W149R':       [-77.516, -81.2264, -84.9364],
    'W149A':       [-69.552, -69.8835, -70.216],
    'H148P/W149R': [-77.093, -88.4217, -99.749],
    'H148A/W149A': [-59.643, -61.7823, -63.921],
}

# ---------------------------
# Data processing
# ---------------------------
original_labels = list(data.keys())
sorted_labels = ['WT'] + sorted([x for x in original_labels if x != 'WT'])
formatted_sorted_labels = [superscript_codons(label) for label in sorted_labels]

df = pd.DataFrame(data)
df_melted = df.melt(var_name='Mutant', value_name='ΔG')

# ANOVA
groups = [df[col] for col in sorted_labels]
anova_stat, anova_p = f_oneway(*groups)
p_text = format_p(anova_p)

# ---------------------------
# Colour-blind-friendly palette
# ---------------------------
palette = ['#0072B2' if label == 'WT' else '#E69F00' for label in sorted_labels]

# ---------------------------
# Plotting
# ---------------------------
fig, ax = plt.subplots(figsize=(9, 6))

sns.barplot(
    data=df_melted,
    x='Mutant',
    y='ΔG',
    order=sorted_labels,
    ax=ax,
    errorbar='sd',
    palette=palette,
    edgecolor='black',
    linewidth=1.2,
    capsize=0.15
)

# Formatting
ax.set_xticks(range(len(formatted_sorted_labels)))
ax.set_xticklabels(formatted_sorted_labels, rotation=0, fontsize=10)
ax.set_ylabel('ΔG (kcal/mol)', fontsize=12)
ax.set_xlabel('')
ax.grid(False)



plt.tight_layout()

# ---------------------------
# Save figure
# ---------------------------
output_pdf = "C:/Users/lrmacha/Downloads/deltaG_plot.pdf"
output_svg = "C:/Users/lrmacha/Downloads/deltaG_plot.svg"

fig.savefig(output_pdf, format='pdf', bbox_inches='tight')
fig.savefig(output_svg, format='svg', bbox_inches='tight')

print(f"Figure saved as PDF: {output_pdf}")
print(f"Figure saved as SVG: {output_svg}")

plt.show()
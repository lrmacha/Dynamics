import matplotlib.pyplot as plt
import numpy as np
import re

# ---------------------------
# Helper function
# ---------------------------
def superscript_mutant(name):
    superscript_map = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    return re.sub(r'(\d+)', lambda x: x.group().translate(superscript_map), name)

# ---------------------------
# Data reconstructed from figure
# ---------------------------
labels = ['W149R', 'WT', 'H148P/W149R', 'H148A/W149A', 'H148P', 'H148A', 'W149A']
values = [2.80, 3.80, 4.10, 4.50, 5.30, 8.50, 10.60]

formatted_labels = [superscript_mutant(x) for x in labels]

# ---------------------------
# Colour-blind-friendly palette
# WT = blue, mutants = orange
# ---------------------------
colors = ['#0072B2' if x == 'WT' else '#E69F00' for x in labels]

# ---------------------------
# Plot
# ---------------------------
fig, ax = plt.subplots(figsize=(10, 5))

bars = ax.bar(
    range(len(labels)),
    values,
    color=colors,
    edgecolor='black',
    linewidth=1.2
)

# Value labels above bars
for bar, val in zip(bars, values):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 0.15,
        f"{val:.2f}",
        ha='center',
        va='bottom',
        fontsize=9
    )

# Axis formatting
ax.set_ylabel('Destabilization Time (ns)')
ax.set_xlabel('')
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(formatted_labels, rotation=45, ha='right')
ax.set_ylim(0, 11.2)

# Clean style to match prior figures
ax.grid(False)
ax.yaxis.grid(False)
ax.xaxis.grid(False)

# Remove top/right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

# ---------------------------
# Save outputs
# ---------------------------
output_pdf = "C:/Users/lrmacha/Downloads/destabilization_time_plot.pdf"
output_svg = "C:/Users/lrmacha/Downloads/destabilization_time_plot.svg"

fig.savefig(output_pdf, format='pdf', bbox_inches='tight')
fig.savefig(output_svg, format='svg', bbox_inches='tight')

print(f"Figure saved as PDF: {output_pdf}")
print(f"Figure saved as SVG: {output_svg}")

plt.show()
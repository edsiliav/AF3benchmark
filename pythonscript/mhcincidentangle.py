# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 15:36:11 2024

@author: edsil
"""

import pandas as pd
import matplotlib.pyplot as plt
import colorsys

def generate_distinct_colors(n):
    colors = []
    golden_ratio_conjugate = (1 + 5 ** 0.5) / 2
    for i in range(n):
        hue = (i * golden_ratio_conjugate) % 1
        r, g, b = colorsys.hsv_to_rgb(hue, 0.7, 0.9)  # Adjust saturation and value as needed
        colors.append((r, g, b))
    return colors

df = pd.read_csv('../data/mergeddata.csv')
grouped_df = df.groupby('MHC Name')

for name, group in grouped_df:
    print(f"Group: {name}")
    print(group)
    print()
    

plt.figure(figsize=(8, 5))

unique_genes = df['MHC Name'].unique()

# Create a color map
colors = generate_distinct_colors(len(unique_genes))


    
# Plot each TRAV gene with a different color
for i, gene in enumerate(unique_genes):
    subset = df[df['MHC Name'] == gene]
    plt.scatter([gene] * len(subset), subset['incident_angle'], color=colors[i], label=gene)
    
    
plt.title('Incident Angles by MHC Name')
plt.xlabel('MHC Name')
plt.ylabel('Incident Angle (degrees)')
plt.xticks(rotation=90)  # Rotate x-axis labels by 90 degrees
#plt.legend(title='TRAV Gene')
plt.tight_layout()
#plt.show()
plt.savefig('../plots/mhcincidentangle.png', format='png', dpi=300)
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:51:10 2024

@author: edsil
"""

import pandas as pd
import matplotlib.pyplot as plt
from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.io import output_notebook
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
grouped_df = df.groupby('TRAV gene')

for name, group in grouped_df:
    print(f"Group: {name}")
    print(group)
    print()
    

plt.figure(figsize=(8, 5))

unique_genes = df['TRAV gene'].unique()

# Create a color map
colors = generate_distinct_colors(len(unique_genes))

# Plot each TRAV gene with a different color
for i, gene in enumerate(unique_genes):
    subset = df[df['TRAV gene'] == gene]
    plt.scatter([gene] * len(subset), subset['incident_angle'], color=colors[i], label=gene)
    
    
plt.title('Incident Angles by TRAV Gene')
plt.xlabel('TRAV Gene')
plt.ylabel('Incident Angle (degrees)')
plt.xticks(rotation=90)  # Rotate x-axis labels by 90 degrees
#plt.legend(title='TRAV Gene')
plt.tight_layout()
#plt.show()
plt.savefig('../plots/incidentangletravgene.png', format='png', dpi=300)

source = ColumnDataSource(df)
output_file("../plots/incidentanglebokeh.html")
# Create a scatter plot
p = figure(title="Scatter Plot of Incident Angle vs Incident Angle",
           x_axis_label='Calculated Incident Angle (degrees)',
           y_axis_label='Database Incident Angle (degrees)')

# Add points using the scatter method
p.scatter(x='incident_angle', y='Incident angle', size=10, source=source, color='blue', alpha=0.6)

# Configure the HoverTool
hover = HoverTool()
hover.tooltips = [("Model Name", "@modelname")]

# Add the HoverTool to the plot
p.add_tools(hover)
#show(p)
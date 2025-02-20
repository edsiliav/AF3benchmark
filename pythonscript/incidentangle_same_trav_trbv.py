# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 13:37:21 2024

@author: edsil
"""

import pandas as pd
import matplotlib.pyplot as plt
import colorsys
from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, HoverTool

from math import pi


def generate_distinct_colors(n):
    colors = []
    golden_ratio_conjugate = (1 + 5 ** 0.5) / 2
    for i in range(n):
        hue = (i * golden_ratio_conjugate) % 1
        r, g, b = colorsys.hsv_to_rgb(hue, 0.7, 0.9)  # Adjust saturation and value as needed
        colors.append((r, g, b))
    return colors

# Load your dataset
# Assuming the file has columns ['TRAV', 'TRBV', 'crossing_angle']
data = pd.read_csv('../data/classI_complexes.csv')

# Create a new column that combines TRAV and TRBV for easier grouping
data['TRAV_TRBV'] = data['TRAV gene'] + '_' + data['TRBV gene']

# Group data by TRAV_TRBV and filter it so they only show if there is more than one data point and calculate crossing angles
data = data.dropna(subset=['Incident angle'])
filtered_data = data.groupby('TRAV_TRBV').filter(lambda x: len(x) > 1)

# Get unique pairs of TRAV_TRBV after filtering
unique_pairs = filtered_data['TRAV_TRBV'].unique()

# Scatter plot with TRAV_TRBV on the x-axis and crossing_angle on the y-axis
# Use `unique TRAV_TRBV` pairs as x-axis labels, assigning each a unique position

colors = generate_distinct_colors(len(unique_pairs))

pair_to_color = {pair: f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}' 
                 for pair, (r, g, b) in zip(unique_pairs, colors)}

x_positions = {pair: i for i, pair in enumerate(unique_pairs)}
data['x_pos'] = data['TRAV_TRBV'].map(x_positions)

plt.figure(figsize=(15, 8))


for i, pair in enumerate(unique_pairs):
    subset = data[data['TRAV_TRBV'] == pair]
    plt.scatter([x_positions[pair]] * len(subset), subset['Incident angle'], color=colors[i], label=pair)
    

# Plot each crossing angle for each TRAV_TRBV pair


# Set custom x-axis ticks and labels
plt.xticks(ticks=range(len(unique_pairs)), labels=unique_pairs, rotation=90)

# Set plot title and labels
plt.title('Incident Angles for Each TRAV-TRBV Pair')
plt.xlabel('TRAV_TRBV Pair')
plt.ylabel('Incident Angle')

# Show plot
plt.tight_layout()
#plt.show()
plt.savefig('../plots/incidentangle_sametravtrbv_filtered.png', format='png', dpi=300)


filtered_data['color'] = filtered_data['TRAV_TRBV'].map(pair_to_color)
source = ColumnDataSource(data=dict(
    TRAV_TRBV=filtered_data['TRAV_TRBV'],
    Incident_angle=filtered_data['Incident angle'],
    PDB_file=filtered_data['PDB ID'],
    MHC_name=filtered_data['MHC Name'],
    color=filtered_data['color'],
    peptide=filtered_data['Epitope']
))

# Set output to file
output_file("../plots/incidentanglesametravtrbvbokeh_filtered.html")


# Create a scatter plot
p = figure(title="Scatter Plot of TRAV_TRBV vs Incident Angle",
           x_axis_label='TRAV_TRBV pair',
           y_axis_label='Incident Angle (degrees)',
           x_range=filtered_data['TRAV_TRBV'].unique(),  # Categorical x-axis
           width=800, height=400)


# Add points using the scatter method
p.scatter(x='TRAV_TRBV', y='Incident_angle', size=10, source=source, color='color', alpha=0.6)

p.xaxis.major_label_orientation = pi/2

# Configure the HoverTool
hover = HoverTool()
hover.tooltips = [("TRAV_TRBV Pair", "@TRAV_TRBV"), ("Incident Angle", "@Incident_angle"),("PDB file", "@PDB_file"),("MHC name", "@MHC_name"),("peptide", "@peptide")]

# Add the HoverTool to the plot
p.add_tools(hover)
show(p)
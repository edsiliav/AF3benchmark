# -*- coding: utf-8 -*-
"""
Name: crossangle_sametrav_trbv_colouredbypep_length.py
Function: Plots the crossing angle for TRAV-TRBV combinations with more than one data point in the dataset coloured by peptide length
Date: 28-01-2025
Author: Edsilia Vermeulen
"""

import pandas as pd
import matplotlib.pyplot as plt
import colorsys
from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, HoverTool
from math import pi


def main():
    # Input and output paths
    input_file = '../data/classI_complexes.csv'
    matplotlib_output = '../plots/pep_length_crossangle_sametravtrbv.png'
    bokeh_output = '../plots/pep_length_crossanglesametravtrbv.html'

    # Load dataset
    data = pd.read_csv(input_file)

    # Combine TRAV and TRBV columns
    data['TRAV_TRBV'] = data['TRAV gene'] + '_' + data['TRBV gene']

    # Compute peptide length
    data['Peptide_Length'] = data['Epitope'].str.len()

    # Filter dataset to show TRAV-TRBV combinations for which there is more than one data point available
    data = data.dropna(subset=['Docking angle'])
    filtered_data = data.groupby('TRAV_TRBV').filter(lambda x: len(x) > 1)

    # Get unique peptide lengths
    unique_lengths = sorted(filtered_data['Peptide_Length'].unique())

    # Generate colors
    colors = generate_distinct_colors(len(unique_lengths))

    # Create a mapping for peptide length to color
    length_to_color = {length: f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}' 
                       for length, (r, g, b) in zip(unique_lengths, colors)}

    # Add a color column to the data
    filtered_data['color'] = filtered_data['Peptide_Length'].map(length_to_color)

    # Create matplotlib scatter plot
    scatter_plot_matplotlib(filtered_data, length_to_color, matplotlib_output)

    # Create bokeh scatter plot
    scatter_plot_bokeh(filtered_data, length_to_color, bokeh_output)


def generate_distinct_colors(n):
    colors = []
    golden_ratio_conjugate = (1 + 5 ** 0.5) / 2
    for i in range(n):
        hue = (i * golden_ratio_conjugate) % 1
        r, g, b = colorsys.hsv_to_rgb(hue, 0.7, 0.9)  # Adjust saturation and value as needed
        colors.append((r, g, b))
    return colors


def scatter_plot_matplotlib(data, length_to_color, output_path):
    plt.figure(figsize=(15, 8))

    for length, color in length_to_color.items():
        subset = data[data['Peptide_Length'] == length]
        plt.scatter(subset['TRAV_TRBV'], subset['Docking angle'], color=color, label=f"Peptide Length {length}", alpha=0.7)

    plt.xticks(rotation=90)
    plt.title('Crossing Angles Colored by Peptide Length')
    plt.xlabel('TRAV_TRBV Pair')
    plt.ylabel('Crossing Angle')
    plt.legend(title="Peptide Length", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_path, format='png', dpi=300)
    plt.show()


def scatter_plot_bokeh(data, length_to_color, output_path):
    source = ColumnDataSource(data=dict(
        TRAV_TRBV=data['TRAV_TRBV'],
        Docking_angle=data['Docking angle'],
        PDB_file=data['PDB ID'],
        Peptide_Length=data['Peptide_Length'].astype(str),  # Convert to string for better tooltip display
        color=data['color'],
        peptide=data['Epitope']
    ))

    output_file(output_path)

    p = figure(title="Scatter Plot of TRAV_TRBV vs Crossing Angle",
               x_axis_label='TRAV_TRBV pair',
               y_axis_label='Crossing Angle (degrees)',
               x_range=data['TRAV_TRBV'].unique(),  # Categorical x-axis
               width=800, height=400)

    p.scatter(x='TRAV_TRBV', y='Docking_angle', size=10, source=source, color='color', alpha=0.6)
    p.xaxis.major_label_orientation = pi / 2

    hover = HoverTool()
    hover.tooltips = [("TRAV_TRBV Pair", "@TRAV_TRBV"), ("Crossing Angle", "@Docking_angle"),
                      ("PDB file", "@PDB_file"), ("Peptide Length", "@Peptide_Length"), ("Peptide", "@peptide")]

    p.add_tools(hover)
    show(p)


if __name__ == "__main__":
    main()

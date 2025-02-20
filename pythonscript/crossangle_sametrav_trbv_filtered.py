# -*- coding: utf-8 -*-
"""
Name: crossangle_sametrav_trbv_filtered.py
Function: Plots the crossing angle for TRAV-TRBV combinations with more than one data point in the dataset
Date: 6-11-2024
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
    matplotlib_output = '../plots/crossangle_sametravtrbv_filtered.png'
    bokeh_output = '../plots/crossanglesametravtrbvbokeh_filtered.html'

    # Load dataset
    data = pd.read_csv(input_file)

    # Combine TRAV and TRBV columns
    data['TRAV_TRBV'] = data['TRAV gene'] + '_' + data['TRBV gene']

    # Filter dataset to show TRAV-TRBV combinations for which there is more than one data point available
    data = data.dropna(subset=['Docking angle'])
    filtered_data = data.groupby('TRAV_TRBV').filter(lambda x: len(x) > 1)

    # Get unique TRAV_TRBV pairs
    unique_pairs = filtered_data['TRAV_TRBV'].unique()

    # Generate colors
    colors = generate_distinct_colors(len(unique_pairs))

    # Create matplotlib scatter plot
    scatter_plot_matplotlib(filtered_data, unique_pairs, colors, matplotlib_output)

    # Create bokeh scatter plot
    scatter_plot_bokeh(filtered_data, unique_pairs, colors, bokeh_output)


def generate_distinct_colors(n):
    colors = []
    golden_ratio_conjugate = (1 + 5 ** 0.5) / 2
    for i in range(n):
        hue = (i * golden_ratio_conjugate) % 1
        r, g, b = colorsys.hsv_to_rgb(hue, 0.7, 0.9)  # Adjust saturation and value as needed
        colors.append((r, g, b))
    return colors


def scatter_plot_matplotlib(data, unique_pairs, colors, output_path):
    pair_to_color = {pair: f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}' 
                     for pair, (r, g, b) in zip(unique_pairs, colors)}

    x_positions = {pair: i for i, pair in enumerate(unique_pairs)}
    data['x_pos'] = data['TRAV_TRBV'].map(x_positions)

    plt.figure(figsize=(15, 8))
    for i, pair in enumerate(unique_pairs):
        subset = data[data['TRAV_TRBV'] == pair]
        plt.scatter([x_positions[pair]] * len(subset), subset['Docking angle'], color=colors[i], label=pair)

    plt.xticks(ticks=range(len(unique_pairs)), labels=unique_pairs, rotation=90)
    plt.title('Crossing Angles for Each TRAV-TRBV Pair')
    plt.xlabel('TRAV_TRBV Pair')
    plt.ylabel('Crossing Angle')
    plt.tight_layout()
    plt.savefig(output_path, format='png', dpi=300)


def scatter_plot_bokeh(data, unique_pairs, colors, output_path):
    pair_to_color = {pair: f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}' 
                     for pair, (r, g, b) in zip(unique_pairs, colors)}

    data['color'] = data['TRAV_TRBV'].map(pair_to_color)
    source = ColumnDataSource(data=dict(
        TRAV_TRBV=data['TRAV_TRBV'],
        Docking_angle=data['Docking angle'],
        PDB_file=data['PDB ID'],
        MHC_name=data['MHC Name'],
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
                      ("PDB file", "@PDB_file"), ("MHC name", "@MHC_name"), ("peptide", "@peptide")]

    p.add_tools(hover)
    show(p)


if __name__ == "__main__":
    main()


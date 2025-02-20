# -*- coding: utf-8 -*-
"""
Name: crossanglebybgene.py
Function: Plots the crossing angle by TRBV gene excluding cases with reverse docking
Date: 29-10-2024
Author: Edsilia Vermeulen
"""

import pandas as pd
import matplotlib.pyplot as plt
import colorsys

def main():
    # Input file and output file paths
    input_file = '../data/classI_complexes.csv'
    output_file = '../plots/crossinganglebybgene2.png'
    
    # Exclude PDB entries which have reverse docking to improve visibility of the figure
    excluded_pdbs = ['6UZ1', '5SWS', '5SWZ', '7JWI']

    # Load the data
    df = pd.read_csv(input_file)

    # Filter out rows with excluded PDB entries
    df = df[~df['PDB ID'].isin(excluded_pdbs)]

    # Plot the data
    plot_crossing_angles(df, output_file)

def generate_distinct_colors(n):
    colors = []
    golden_ratio_conjugate = (1 + 5 ** 0.5) / 2
    for i in range(n):
        hue = (i * golden_ratio_conjugate) % 1
        r, g, b = colorsys.hsv_to_rgb(hue, 0.7, 0.9)  # Adjust saturation and value as needed
        colors.append((r, g, b))
    return colors

def plot_crossing_angles(df, output_file):
    #Create a scatter plot
    plt.figure(figsize=(8, 5))

    unique_genes = df['TRBV gene'].unique()

    # Create a color map
    colors = generate_distinct_colors(len(unique_genes))

    # Plot each TRBV gene with a different color
    for i, gene in enumerate(unique_genes):
        subset = df[df['TRBV gene'] == gene]
        plt.scatter([gene] * len(subset), subset['Docking angle'], color=colors[i], label=gene)
        
    plt.title('Crossing Angles by TRBV Gene')
    plt.xlabel('TRBV Gene')
    plt.ylabel('Crossing Angle (degrees)')
    plt.xticks(rotation=90)  # Rotate x-axis labels by 90 degrees
    plt.tight_layout()
    plt.savefig(output_file, format='png', dpi=300)

if __name__ == "__main__":
    main()
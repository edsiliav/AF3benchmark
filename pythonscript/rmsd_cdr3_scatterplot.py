# -*- coding: utf-8 -*-
"""
Name: rmsd_cdr3_scatterplot.py
Function: Creates a scatter plot of the RMSD data of the CDR3 loop coloured based on the model rank
Date: 29-01-2025
Name: Edsilia Vermeulen
"""

import pandas as pd
import matplotlib.pyplot as plt
import re


def main():
    # Input: File path
    file_path = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\cdr3_rmsd_results\cdr3_rmsd.csv"
    output_path = "../plots/rmsd_pdb_scatter_colored_CDR3B.png"

    # Load dataset
    df = pd.read_csv(file_path)

    # Process data
    df['Base_PDB'] = df['PDB_ID'].apply(lambda x: re.sub(r'_\d+$', '', x).upper())
    df['Model'] = df['PDB_ID'].apply(lambda x: x.split('_')[-1])

    # Assign a unique color for each model rank (_0, _1, etc.)
    model_colors = {'0': 'red', '1': 'blue', '2': 'green', '3': 'purple', '4': 'orange'}

    # Create scatter plot
    plt.figure(figsize=(12, 6))
    
    # Plot each point, grouped by PDB and colored by model
    for model, color in model_colors.items():
        subset = df[df['Model'] == model]
        plt.scatter(subset['Base_PDB'], subset['RMSD_CDR3_B'], alpha=0.7, color=color, label=f"Model _{model}")

    plt.xticks(rotation=90)
    plt.xlabel("PDB ID")
    plt.ylabel("RMSD CDR3 B")
    plt.title("RMSD of CDR3 B loops Distribution by PDB")
    plt.legend(title="Model Type")
    plt.grid(True, linestyle="--", alpha=0.6)

    # Save and show
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.show()

if __name__ == "__main__":
    main()


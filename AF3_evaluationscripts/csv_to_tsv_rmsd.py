# -*- coding: utf-8 -*-
"""
Name: cdr3loop_csv_to_tsv.py
Function: Convert CSV files of RMSD to TSV files and arrange files to make success plots
Date: 10-12-2024
Author: Edsilia Vermeulen
"""

import pandas as pd
import os

def main():
    # Input & Output Settings
    data_file = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\rmsd_good\dockq_results.tsv"
    output_folder = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\rmsd_good"

    # Ensure output directory exists
    os.makedirs(output_folder, exist_ok=True)

    # Load the dataset
    try:
        df = pd.read_csv(data_file, sep="\t", header=None, names=["PDB ID", "Fnat", "LRMSD", "IRMSD"])
    except Exception as e:
        print(f"Error reading the file: {e}")
        return

    # Convert numeric columns to float
    try:
        df[["Fnat", "LRMSD", "IRMSD"]] = df[["Fnat", "LRMSD", "IRMSD"]].astype(float)
    except Exception as e:
        print(f"Error converting columns to float: {e}")
        return

    # Group by PDB ID
    grouped = df.groupby("PDB ID")

    # Output files
    output_files = {
        "Fnat": f"{output_folder}/fnat_scores.tsv",
        "LRMSD": f"{output_folder}/lrmsd_scores.tsv",
        "IRMSD": f"{output_folder}/irmsd_scores.tsv"
    }

    # Write grouped data to files
    for column, output_file in output_files.items():
        write_to_file(grouped, column, output_file)

    print("Files created:", ", ".join(os.path.basename(f) for f in output_files.values()))

def write_to_file(grouped_data, column, output_file):
    """Writes grouped data to a file."""
    try:
        with open(output_file, "w") as f:
            for pdb_id, group in grouped_data:
                values = group[column].tolist()
                f.write(f"{pdb_id}\t" + "\t".join(map(str, values)) + "\n")
    except Exception as e:
        print(f"Error writing to file {output_file}: {e}")

if __name__ == "__main__":
    main()


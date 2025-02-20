# -*- coding: utf-8 -*-
"""
Name: cdr3loop_csv_to_tsv.py
Function: Convert CSV files of CDR3 loops of the TCR alpha and beta chains to TSV files
Date: 6-1-2025
Author: Edsilia Vermeulen
"""
import pandas as pd

def main():
    # Load the CSV file
    input_file = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\cdr3_rmsd_results\cdr3_rmsd.csv"  
    
    # Specify output file paths
    output_file_a = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\cdr3_rmsd_results\rmsd_cdr3_A.tsv"
    output_file_b = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\cdr3_rmsd_results\rmsd_cdr_B.tsv"

    # Read the CSV data
    data = pd.read_csv(input_file)

    # Create the first TSV file with PDB_ID and RMSD_CDR3_A
    data[['PDB_ID', 'RMSD_CDR3_A']].to_csv(output_file_a, sep='\t', index=False)

    # Create the second TSV file with PDB_ID and RMSD_CDR3_B
    data[['PDB_ID', 'RMSD_CDR3_B']].to_csv(output_file_b, sep='\t', index=False)

    print(f"TSV files created successfully:\n{output_file_a}\n{output_file_b}")

if __name__ == "__main__":
    main()


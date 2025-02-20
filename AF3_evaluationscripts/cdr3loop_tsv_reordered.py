# -*- coding: utf-8 -*-
"""
Name: cdr3loop_tsv_reordered.py
Function: Organising the data to create a success plot
Date: 6-1-2025
Author: Edsilia Vermeulen
"""


def main():
    # Specify the file paths
    input_file = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\cdr3_rmsd_results\rmsd_cdr_B.tsv"
    output_file = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\cdr3_rmsd_results\rmsd_cdr3_B_good.tsv"  

    # Read the input TSV file
    with open(input_file, mode='r') as infile:
        # Read the data (skip the header)
        lines = infile.readlines()[1:]  # Skip the first header line
        
        # Create a dictionary to store the RMSD values grouped by PDB
        pdb_dict = {}

        # Iterate through the lines and group by PDB_ID
        for line in lines:
            pdb_id, rmsd_value = line.strip().split('\t')
            pdb_base = pdb_id.split('_')[0]  # Extract the base PDB ID (e.g., '7l1d')
            
            if pdb_base not in pdb_dict:
                pdb_dict[pdb_base] = []
            
            pdb_dict[pdb_base].append(rmsd_value)  # Append the RMSD value

    # Write the transformed data to the output TSV file
    with open(output_file, mode='w') as outfile:
        # Write the data (no header)
        for pdb_base, rmsd_values in pdb_dict.items():
            row = [pdb_base] + rmsd_values
            outfile.write('\t'.join(row) + '\n')

    print("Output saved to", output_file)

# Run the main function
if __name__ == "__main__":
    main()

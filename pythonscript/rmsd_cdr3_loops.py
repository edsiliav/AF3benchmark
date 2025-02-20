# -*- coding: utf-8 -*-
"""
Name: rmsd_cdr3_loops.py
Function: Calculates RMSD of CDR3 loops which are IMGT numbered
Date: 11-12-2024
Author: Edsilia Vermeulen
"""

from pdb2sql import pdb2sql
import numpy as np
import os
import csv



def main():
    model_file = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\merged_renumbered_models\7rk7_0_merged.pdb"
    reference_file = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\merged_renumbered_references\7rk7_merged.pdb"
    
    pdb_id = os.path.basename(model_file).split('_merged')[0]
    
    #opens the files in pdb2sql
    db_model = pdb2sql(model_file)
    db_reference = pdb2sql(reference_file)
    

    #extract coordinates of the CDR3 loops for the IMGT numbered models and references, structures should already be aligned
    model_cdr3_a = db_model.get('x,y,z', chainID= ['A'], resSeq=['105', '117'], name=['CA'])
    model_cdr3_b = db_model.get('x,y,z', chainID= ['B'], resSeq=['105', '117'], name=['CA'])
    reference_cdr3_a = db_reference.get('x,y,z', chainID= ['D'], resSeq=['105', '117'], name=['CA'])
    reference_cdr3_b = db_reference.get('x,y,z', chainID= ['E'], resSeq=['105', '117'], name=['CA'])

    # Convert the coordinates to numpy arrays to calculate RMSD
    np_model_cdr3_a = np.array(model_cdr3_a)
    np_reference_cdr3_a = np.array(reference_cdr3_a)
    np_model_cdr3_b = np.array(model_cdr3_b)
    np_reference_cdr3_b = np.array(reference_cdr3_b)
    

    
    # Calculate RMSD
    rmsd_cdr3_a = round(calc_rmsd(np_model_cdr3_a, np_reference_cdr3_a), 3)
    rmsd_cdr3_b = round(calc_rmsd(np_model_cdr3_b, np_reference_cdr3_b), 3)


    # Save results to CSV
    save_results_to_csv(pdb_id, rmsd_cdr3_a, rmsd_cdr3_b)
    

def calc_rmsd(coords1, coords2):
    """Calc rmsd without superposition

    coord1, coord2: numpy array with coordinates [(x1,y1,z1), (xi,yi,zi)]
    Return: rmsd(float): the rmsd.
    """
    rmsd = np.sqrt(((((coords1 - coords2)** 2))*3).mean())
    return rmsd

def save_results_to_csv(pdb_id, rmsd_cdr3_a, rmsd_cdr3_b):
    """Save the PDB ID and RMSD values to a CSV file."""
    csv_file = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\cdr3_rmsd_results\cdr3_rmsd.csv"
    print(f"Saving results to {csv_file}...")

    # Check if the file exists
    file_exists = os.path.isfile(csv_file)

    # Open the file in append mode
    with open(csv_file, mode='a', newline='') as file:
        writer = csv.writer(file)

        # Write the header if the file is new
        if not file_exists:
            writer.writerow(["PDB_ID", "RMSD_CDR3_A", "RMSD_CDR3_B"])

        # Write the data row
        writer.writerow([pdb_id, rmsd_cdr3_a, rmsd_cdr3_b])


if __name__ == "__main__":
    main()

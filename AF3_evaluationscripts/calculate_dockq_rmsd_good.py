# -*- coding: utf-8 -*-
"""
Name: calculate_dockq_rmsd.py
Function: Uses DockQ to calculate fnat, LRMSD, and iRMSD of predicted models based on a reference structure. 
Date: 4-12-2024
Author: Edsilia Vermeulen
"""


import os
import csv
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

def main():
    # Paths to predicted structures and reference files
    reres_files_path = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\reres_files"
    reres_references_path = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\reres_references"

    # New directory for saving results
    output_dir = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\rmsd_good"
    os.makedirs(output_dir, exist_ok=True)  # Create the directory if it doesn't exist

    # Output CSV file in the new directory
    output_csv = os.path.join(output_dir, "dockq_results_good.csv")
    # Run the processing function
    process_files(reres_files_path, reres_references_path, output_csv)


# Function to process each PDB file and run DockQ
def process_files(reres_files_path, reres_references_path, output_csv):
    results = []
    for root, dirs, files in os.walk(reres_files_path):
        for file in files:
            if file.endswith(".pdb"):
                pdb_id = os.path.basename(root)  # Extract PDB ID from folder name
                ref_file = os.path.join(reres_references_path, f"{pdb_id}_ref_reres.pdb")
                model_file = os.path.join(root, file)

                if os.path.isfile(ref_file):
                    print(f"Processing model: {model_file} with reference: {ref_file}")
                    model = load_PDB(model_file)
                    native = load_PDB(ref_file)

                    chain_map = {"A": "A", "D": "D"} #compares chain A of the model to chain A of native, and chain D of model to chain D of native

                    result_tuple = run_on_all_native_interfaces(model, native, chain_map=chain_map)
                    print(f"Result for {file}: {result_tuple}")  # Debug print

                    if isinstance(result_tuple, tuple) and isinstance(result_tuple[0], dict):
                        interface_results = result_tuple[0]
                        dockq_score = round(result_tuple[1], 3)  # Round overall DockQ score

                        # Extract and round DockQ metrics from the first interface found (e.g., 'AD')
                        first_interface = next(iter(interface_results.values()))
                        fnat = round(first_interface.get('fnat', 0.0), 3)
                        lrmsd = round(first_interface.get('LRMSD', 0.0), 3)
                        irmsd = round(first_interface.get('iRMSD', 0.0), 3)

                        results.append([pdb_id, file, dockq_score, fnat, lrmsd, irmsd])
                    else:
                        print(f"Unexpected result format: {result_tuple}")
                        continue

    # Write results to a CSV file
    with open(output_csv, mode="w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["PDB ID", "Model File", "DockQ Score", "Fnat", "LRMSD", "iRMSD"])
        writer.writerows(results)
    
    print(f"Results saved to {output_csv}")

if __name__ == "__main__":
    main()


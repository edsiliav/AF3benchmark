# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 16:01:57 2024

@author: edsil
"""

from Bio.PDB import MMCIFParser, PDBIO

# Path to the CIF file
cif_file = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\alphafold\8wul\fold_8wul_2024_11_20_14_11_model_4.cif"
# Desired output PDB file
pdb_file = r"C:\Users\edsil\Documents\coevolutionproject\difficult_cases\alphafold\pdb_files\8wul_4.pdb"

# Initialize the MMCIFParser and PDBIO
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("structure", cif_file)

# Create a PDBIO object to write the structure to a PDB file
io = PDBIO()
io.set_structure(structure)
io.save(pdb_file)

print(f"Conversion complete. The PDB file is saved as {pdb_file}")

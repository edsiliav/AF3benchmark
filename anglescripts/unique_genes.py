# -*- coding: utf-8 -*-
"""
Name: unique_genes.py
Function: Determine unique TRAV, TRBV, MHC or TRAV_TRBV genes in a dataset for PDB files released after AF3 cut-off date
Date: 29-11-2024
Author: Edsilia Vermeulen
"""


import pandas as pd

df = pd.read_csv('../data/classI_complexes.csv')
df['TRAV_TRBV'] = df['TRAV gene'] + '_' + df['TRBV gene']
df['date'] = pd.to_datetime(df['Release<BR>date'], format='%Y-%m-%d')

filtered_df = df.loc[(df['date'] > '2021-09-30')].copy()
training_set = df.loc[(df['date'] < '2021-09-30')]

# Unique TRAV genes

unique_trav = df['TRAV gene'].unique()
unique_trav_df = pd.DataFrame(unique_trav, columns=['Unique TRAV'])

merged_trav = pd.merge(unique_trav_df, df, left_on='Unique TRAV', right_on='TRAV gene', suffixes=('_left', '_right'))
filtered_trav = merged_trav[~merged_trav["Unique TRAV"].isin(training_set["TRAV gene"])]


# Unique TRBV genes

unique_trbv = df['TRBV gene'].unique()
unique_trbv_df = pd.DataFrame(unique_trbv, columns=['Unique TRBV'])

merged_trbv = pd.merge(unique_trbv_df, df, left_on='Unique TRBV', right_on='TRBV gene', suffixes=('_left', '_right'))
filtered_trbv = merged_trbv[~merged_trbv["Unique TRBV"].isin(training_set["TRBV gene"])]

#Unique MHCs

unique_mhc = df['MHC Name'].unique()
unique_mhc_df = pd.DataFrame(unique_mhc, columns=['Unique MHC'])

merged_mhc = pd.merge(unique_mhc_df, df, left_on='Unique MHC', right_on='MHC Name', suffixes=('_left', '_right'))
filtered_mhc = merged_mhc[~merged_mhc["Unique MHC"].isin(training_set["MHC Name"])]

#filtered_df.to_csv('../data/uniqe_travtrbvmhc_data.csv')


# Unique TRAV-TRBV gene pairs

unique_vgene = df['TRAV_TRBV'].unique()
unique_vgene_df = pd.DataFrame(unique_vgene, columns=['Unique TRAV_TRBV'])

merged_vgene = pd.merge(unique_vgene_df, df, left_on='Unique TRAV_TRBV', right_on='TRAV_TRBV', suffixes=('_left', '_right'))
filtered_vgene = merged_vgene[~merged_vgene["Unique TRAV_TRBV"].isin(training_set["TRAV_TRBV"])]



# Create a new column for each unique type (TRAV, TRBV, MHC, and TRAV_TRBV)
filtered_df['Unique_TRAV'] = filtered_df['TRAV gene'].apply(lambda x: 'X' if x in unique_trav and x not in training_set['TRAV gene'].values else '')
filtered_df['Unique_TRBV'] = filtered_df['TRBV gene'].apply(lambda x: 'X' if x in unique_trbv and x not in training_set['TRBV gene'].values else '')
filtered_df['Unique_MHC'] = filtered_df['MHC Name'].apply(lambda x: 'X' if x in unique_mhc and x not in training_set['MHC Name'].values else '')
filtered_df['Unique_TRAV_TRBV'] = filtered_df['TRAV_TRBV'].apply(lambda x: 'X' if x in unique_vgene and x not in training_set['TRAV_TRBV'].values else '')

# Select the relevant columns and arrange them in the desired order
output_df = filtered_df[['PDB ID', 'Release<BR>date', 'MHC Name', 'TRAV gene', 'TRBV gene', 
                         'Docking angle', 'Unique_TRAV', 'Unique_TRBV', 'Unique_MHC', 'Unique_TRAV_TRBV']]



# Save the final DataFrame to a CSV file
output_df.to_csv('../data/unique_travtrbvmhc_data.csv', index=False)








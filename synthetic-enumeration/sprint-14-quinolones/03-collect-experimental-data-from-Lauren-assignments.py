#!/bin/env python

"""
Collect experimental data from Lauren's reassignments via CSV file

"""

import numpy as np
import json
import math
import itertools
import datetime
from rich.progress import track
from openeye import oechem

xchem_project = 'Mpro'
creator = 'John Chodera <john.chodera@choderalab.org>'
creation_date = datetime.datetime.now()
prefix = 'sprint-12'
description = 'COVID Moonshot Sprint 12 for optimizing 5-spiro compounds'

csv_filename = 'experimental-data/Fl_agg_data_all_data_11_01_2022_11_13_20-cleaned-reassigned_isomers.csv'


#
# Now pull in all submitted designs
#

def smiles_is_racemic(suspected_smiles):
    """
    Return True if compound is racemic.

    Examples:
    "CNC(=O)CN1Cc2ccc(Cl)cc2[C@@]2(CCN(c3cncc4c3CCCC4)C2=O)C1 |o1:14|" : compound is enantiopure, but stereochemistry is uncertain
    "CNC(=O)CN1Cc2ccc(Cl)cc2[C@@]2(CCN(c3cncc4c3CCCC4)C2=O)C1" : compound is enantiopure, stereochemistry is certain
    "CNC(=O)CN1Cc2ccc(Cl)cc2[C]2(CCN(c3cncc4c3CCCC4)C2=O)C1" : compound is racemic
    """
    smiles = suspected_smiles.split()[0] # truncate suffix
    return stereochemistry_is_uncertain(smiles)

def stereochemistry_is_uncertain(suspected_smiles):
    """
    Return True if there is uncertainty in the enantiopure compound or mixture is racemic.

    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

    rdmol = Chem.MolFromSmiles(suspected_smiles)
    smi_list = []
    opts = StereoEnumerationOptions(unique=True)
    isomers = tuple(EnumerateStereoisomers(rdmol, options=opts))
    for smi in sorted(Chem.MolToSmiles(isomer, isomericSmiles=True) for isomer in isomers):
        smi_list.append(smi)

    if len(smi_list) > 1:
        return True
    else:
        return False

# Read all submitted designs
print('Reading CSV export...')
compounds_with_experimental_data = list()
# Drop columns that cause trouble for OpenEye
import pandas as pd
df = pd.read_csv(csv_filename, dtype=str)
# Drop columns
#drop_columns = []
#df.drop(columns=drop_columns, inplace=True)
# Replace suspected_SMILES with SMILES
#df['suspected_SMILES'].fillna(df['SMILES'], inplace=True)
# Exchange columns so suspected_SMILES is first
#title_column_index = df.columns.get_loc("Canonical PostEra ID")
#smiles_column_index = df.columns.get_loc("suspected_SMILES")
#cols = df.columns.tolist()
#cols = cols[smiles_column_index:(smiles_column_index+1)] + cols[title_column_index:(title_column_index+1)] + cols[:]
#df = df[cols]
# Replace < and > with limits
#df.applymap(lambda x: str(x))
#df.applymap(lambda x: 0.050 if "<" in str(x) else x)
#df.applymap(lambda x: 99.0 if ">" in str(x) else x)
# Eliminate stuff after spaces
#df = df.applymap(lambda x: str(x).split()[0])

ncompounds_dropped_due_to_uncertain_stereochemistry = 0
ncompounds_racemic = 0

# Iterate over molecules
# Fields: compound_name,compound_structure,measurement,qualifier,reassigned_structure
# Format: PostEra ID,SMILES,pIC50,comparator,reassigned_structure
delta_pIC50 = 0.2 # 95% CI is this many units in either direction
from fah_xchem.schema import ExperimentalCompoundData, ExperimentalCompoundDataUpdate
for index, row in df.iterrows():
    row = row.to_dict()
    suspected_smiles = row['compound_structure']
    compound_id = row['compound_name']
    is_racemic = smiles_is_racemic(suspected_smiles)
    # Skip inequalities
    if row['qualifier'] != '=':
        continue

    pIC50 = float(row['measurement'])
    pIC50_lower = pIC50 - delta_pIC50
    pIC50_upper = pIC50 + delta_pIC50

    # Canonicalize with OpenEye SMILES
    suspected_smiles = suspected_smiles.split()[0] # truncate stuff after whitespace
    oemol = oechem.OEGraphMol()
    oechem.OESmilesToMol(oemol, suspected_smiles)
    suspected_smiles = oechem.OEMolToSmiles(oemol)

    experimental_data = dict()
    experimental_data['pIC50'] = pIC50
    experimental_data['pIC50_lower'] = pIC50_lower
    experimental_data['pIC50_upper'] = pIC50_upper

    if is_racemic:
        ncompounds_racemic += 1

    # Store compound experimental data
    experimental_compound_data = ExperimentalCompoundData(
        compound_id=compound_id,
        smiles=suspected_smiles,
        is_racemic=is_racemic,
        experimental_data=experimental_data,
    )
    compounds_with_experimental_data.append(experimental_compound_data)

print(f'{len(compounds_with_experimental_data)} measurements read and retained')
print(f'{ncompounds_dropped_due_to_uncertain_stereochemistry} enantiopure compounds with uncertain stereochemistry dropped.')
print(f'{ncompounds_racemic} compounds assayed as racemates')

dataset = ExperimentalCompoundDataUpdate(compounds=compounds_with_experimental_data)
print(f'There are {len(compounds_with_experimental_data)} compounds in this sprint with in-range IC50 measurements')

# Write JSON
def write_json(compound_series, json_filename):
    print(f'Writing JSON to {json_filename}')
    if '.bz2' in json_filename:
        import bz2
        with bz2.open(json_filename, "wt") as f:
            f.write(compound_series.json())
    elif '.gz' in json_filename:
        import gzip
        with gzip.open(json_filename, "wt") as f:
            f.write(compound_series.json())
    else:
        with open(json_filename, "wt") as f:
            f.write(compound_series.json())

import os
os.makedirs('json', exist_ok=True)
print(f'Generating experimental data JSON for {prefix}...')
json_filename = f'json/{prefix}-experimental-data.json' # output filename
write_json(dataset, json_filename)

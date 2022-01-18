#!/bin/env python

"""
Collect experimental data from GitHub repository and compile into a separate experimental-data.json file

https://github.com/postera-ai/COVID_moonshot_submissions

Experimental data available from: https://github.com/postera-ai/COVID_moonshot_submissions/blob/master/covid_submissions_all_info.csv

Presumed SMILES is written instead of registered SMILES

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
prefix = 'sprint-11A'
description = 'COVID Moonshot Sprint 11A for optimizing 5-spiro compounds'

# TODO: Auto-download and timestamp
#cdd_csv_filename = 'experimental-data/CDD CSV Export - 2022-01-01 22h24m13s.csv'
cdd_csv_filename = 'experimental-data/CDD CSV Export - 2022-01-15 23h16m57s.csv'

def replace_non_numerical_sddata(mol, tagname, value):
    """Replace any non-numerical SDData values with the specified value.

    Parameters
    ----------
    mol : openeye.oechem.OEMol
        The molecule
    tagname : str
        The SDData tag name
    value : float
        The value to replace the SDData by
    """
    from openeye import oechem
    try:
        numerical_value = float(oechem.OEGetSDData(mol, tagname))
    except ValueError as e:
        oechem.OESetSDData(mol, tagname, str(value))

#
# Define data model for export
# TODO: Move this into fah-xchem
#

from fah_xchem.schema import Model, Field
from typing import List, Dict

class ExperimentalCompoundData(Model):
    compound_id: str = Field(
        None, description="The unique compound identifier (PostEra or enumerated ID)"
    )

    smiles: str = Field(
        None,
        description="OpenEye canonical isomeric SMILES string defining suspected SMILES of racemic mixture (with unspecified stereochemistry) or specific enantiopure compound (if is_racemic=False); may differ from what is registered under compound_id.",
    )

    is_racemic: bool = Field(
        False,
        description="If True, this experiment was performed on a racemate; if False, the compound was enantiopure.",
    )

    experimental_data: Dict[str, float] = Field(
        dict(), description='Experimental data fields, including "pIC50" and uncertainty (either "pIC50_stderr" or  "pIC50_{lower|upper}"',
    )

class ExperimentalCompoundDataUpdate(Model):
    """A bundle of experimental data for compounds (racemic or enantiopure)."""
    compounds: List[ExperimentalCompoundData]

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
print('Reading CDD export...')
compounds_with_experimental_data = list()
# Drop columns that cause trouble for OpenEye
import pandas as pd
df = pd.read_csv(cdd_csv_filename, dtype=str)
# Drop columns
drop_columns = []
df.drop(columns=drop_columns, inplace=True)
# Replace suspected_SMILES with SMILES
df['suspected_SMILES'].fillna(df['SMILES'], inplace=True)
# Exchange columns so suspected_SMILES is first
title_column_index = df.columns.get_loc("Canonical PostEra ID")
smiles_column_index = df.columns.get_loc("suspected_SMILES")
cols = df.columns.tolist()
cols = cols[smiles_column_index:(smiles_column_index+1)] + cols[title_column_index:(title_column_index+1)] + cols[:]
df = df[cols]
# Replace < and > with limits
#df.applymap(lambda x: str(x))
#df.applymap(lambda x: 0.050 if "<" in str(x) else x)
#df.applymap(lambda x: 99.0 if ">" in str(x) else x)
# Eliminate stuff after spaces
#df = df.applymap(lambda x: str(x).split()[0])

ncompounds_dropped_due_to_uncertain_stereochemistry = 0
ncompounds_racemic = 0

# Iterate over molecules
for index, row in df.iterrows():
    row = row.to_dict()
    suspected_smiles = row['suspected_SMILES']
    compound_id = row['Canonical PostEra ID']
    is_racemic = smiles_is_racemic(suspected_smiles)
    IC50 = row['ProteaseAssay_Fluorescence_Dose-Response_Weizmann: IC50 (µM)']
    IC50_lower = row["ProteaseAssay_Fluorescence_Dose-Response_Weizmann: IC50 CI (Lower) (µM)"]
    IC50_upper = row["ProteaseAssay_Fluorescence_Dose-Response_Weizmann: IC50 CI (Upper) (µM)"]

    # DEBUG: Overrides for important compounds
    if (compound_id == 'VLA-UCB-50c39ae8-2'):
        suspected_smiles = 'c1ccc2c(c1)cncc2N3C(=O)C[C@@]4(C3=O)CCOc5c4cc(cc5)Cl'

    # Drop compounds that are enantiopure but stereochemistry is uncertain
    if (not is_racemic) and stereochemistry_is_uncertain(suspected_smiles):
        print(f'Dropping {compound_id} because enantiomeric stereochemistry is uncertain')
        ncompounds_dropped_due_to_uncertain_stereochemistry += 1
        continue

    try:
        # TODO: Handle case where IC50 is given as '< 0.0500' or '> 99' to expand our dynamic range

        IC50 = float(IC50) * 1.0e-6
        IC50_lower = float(IC50_lower) * 1.0e-6
        IC50_upper = float(IC50_upper) * 1.0e-6
    except ValueError as e:
        print(e)
        continue

    # Canonicalize with OpenEye SMILES
    suspected_smiles = suspected_smiles.split()[0] # truncate stuff after whitespace
    oemol = oechem.OEGraphMol()
    oechem.OESmilesToMol(oemol, suspected_smiles)
    suspected_smiles = oechem.OEMolToSmiles(oemol)

    experimental_data = dict()
    experimental_data['pIC50'] = - np.log10(IC50)
    experimental_data['pIC50_lower'] = - np.log10(IC50_upper)
    experimental_data['pIC50_upper'] = - np.log10(IC50_lower)

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

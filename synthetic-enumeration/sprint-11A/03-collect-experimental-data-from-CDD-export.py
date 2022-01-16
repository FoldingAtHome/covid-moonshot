#!/bin/env python

"""
Collect experimental data from GitHub repository and compile into a separate experimental-data.json file

https://github.com/postera-ai/COVID_moonshot_submissions

Experimental data available from: https://github.com/postera-ai/COVID_moonshot_submissions/blob/master/covid_submissions_all_info.csv

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
cdd_csv_filename = 'experimental-data/CDD CSV Export - 2022-01-01 22h24m13s.csv'

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

# Read all submitted designs: Compounds with the key substructure will be retained
print('Reading submitted designs...')
compounds = dict()
# Drop columns that cause trouble for OpenEye
import pandas as pd
df = pd.read_csv(cdd_csv_filename, dtype=str)
# Drop columns
drop_columns = []
df.drop(columns=drop_columns, inplace=True)
# Exchange columns so suspected_SMILES is first
title_column_index = df.columns.get_loc("Canonical PostEra ID")
smiles_column_index = df.columns.get_loc("suspected_SMILES")
cols = df.columns.tolist()
cols = cols[smiles_column_index:(smiles_column_index+1)] + cols[title_column_index:(title_column_index+1)] + cols[:]
df = df[cols]
# Replace suspected_SMILES with SMILES
df['suspected_SMILES'].fillna(df['SMILES'], inplace=True)
# Replace < and > with limits
#df.applymap(lambda x: str(x))
#df.applymap(lambda x: 0.050 if "<" in str(x) else x)
#df.applymap(lambda x: 99.0 if ">" in str(x) else x)
# Write CSV file
import tempfile
import csv
with tempfile.NamedTemporaryFile(suffix='.csv') as csv_file:
    df.to_csv(csv_file.name, header=True, index=False, quoting=csv.QUOTE_NONNUMERIC)
    # Read file
    with oechem.oemolistream(csv_file.name) as ifs:
        mol = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, mol):            
            try:
                float(oechem.OEGetSDData(mol, "ProteaseAssay_Fluorescence_Dose-Response_Weizmann: IC50 (µM)"))
            except ValueError as e:
                print(e)
                continue

            replace_non_numerical_sddata(mol, "ProteaseAssay_Fluorescence_Dose-Response_Weizmann: IC50 CI (Lower) (µM)", 0.050)
            replace_non_numerical_sddata(mol, "ProteaseAssay_Fluorescence_Dose-Response_Weizmann: IC50 CI (Upper) (µM)", 100.0)
            title = oechem.OEGetSDData(mol, "Canonical PostEra ID")
            mol.SetTitle(title)
            # Store the molecule
            compounds[title] = mol.CreateCopy()
print(f'{len(compounds)} molecules read')

# Compile experimental data on compounds in this sprint
compounds_smi_filename = f'docked/{prefix}-compounds.smi' # compounds in this sprint

from fah_xchem.schema import Model, CompoundMetadata
from typing import List
class ExperimentalCompoundData(Model):
    """Experimental data for compounds."""
    compounds: List[CompoundMetadata]

from openeye import oechem
print('Processing compounds...')
from openeye import oechem
import numpy as np
compounds_with_experimental_data = list()
# TODO: Use CSV filename instead of SMI, since we need it to retain experimental data
with oechem.oemolistream(compounds_smi_filename) as ifs:
    for oemol in ifs.GetOEGraphMols():
        compound_id = oemol.GetTitle()
        smiles = oechem.OEMolToSmiles(oemol)
        # Check if this molecule has experimental data
        if compound_id in compounds:
            experimental_data = dict()
            compound = compounds[compound_id]
            
            IC50 = float(oechem.OEGetSDData(compound, "ProteaseAssay_Fluorescence_Dose-Response_Weizmann: IC50 (µM)")) * 1e-6 # Molar
            IC50_lower = float(oechem.OEGetSDData(compound, "ProteaseAssay_Fluorescence_Dose-Response_Weizmann: IC50 CI (Lower) (µM)")) * 1e-6 # Molar
            IC50_upper = float(oechem.OEGetSDData(compound, "ProteaseAssay_Fluorescence_Dose-Response_Weizmann: IC50 CI (Upper) (µM)")) * 1e-6 # Molar
            
            experimental_data['pIC50'] = - np.log10(IC50)
            experimental_data['pIC50_lower'] = - np.log10(IC50_upper)
            experimental_data['pIC50_upper'] = - np.log10(IC50_lower)

            # Store compound experimental data
            compound_metadata = CompoundMetadata(
                compound_id=compound_id,
                smiles=smiles,
                experimental_data=experimental_data,
            )
            compounds_with_experimental_data.append(compound_metadata)
        
dataset = ExperimentalCompoundData(compounds=compounds_with_experimental_data)
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


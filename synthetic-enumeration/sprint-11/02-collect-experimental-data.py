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
prefix = 'sprint-11-2021-12-26'
description = 'COVID Moonshot Sprint 11 for optimizing spiro compounds'

# TODO: Auto-download and timestamp
submissions_csv_filename = 'submissions/submissions-2021-12-26.csv.gz'

# Read all submitted designs: Compounds with the key substructure will be retained
print('Reading submitted designs...')
compounds = dict()
# Drop columns that cause trouble for OpenEye
import pandas as pd
drop_columns = ['Submission Rationale', 'Submission Notes']
df = pd.read_csv(submissions_csv_filename, dtype=str)
df.drop(columns=drop_columns, inplace=True)
import tempfile
with tempfile.NamedTemporaryFile(suffix='.csv') as csv_file:
    df.to_csv(csv_file.name, header=True, index=False)
    # Read file
    with oechem.oemolistream(csv_file.name) as ifs:
        mol = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, mol):
            # Clear SD tags
            #oechem.OEClearSDData(mol)
            # Store the molecule
            compounds[mol.GetTitle()] = mol.CreateCopy()
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
            if oechem.OEHasSDData(compound,'f_avg_IC50'):
                IC50 = oechem.OEGetSDData(compound, 'f_avg_IC50')
                if IC50 != '':
                    IC50 = float(IC50)
                    if IC50 < 99: # dynamic range of assay
                        IC50 *= 1.0e-6 # convert to uM
                        pIC50 = - np.log10(IC50)
                        experimental_data['pIC50'] = pIC50
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


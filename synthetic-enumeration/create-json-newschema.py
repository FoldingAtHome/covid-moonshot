#!/bin/env python

"""
Create fah-xhcem compound set JSON file using new schema

https://github.com/choderalab/fah-xchem/pull/77

"""

import numpy as np
import json
import math
import itertools
import datetime
from rich.progress import track

xchem_project = 'Mpro'
creator = 'John Chodera <john.chodera@choderalab.org>'
ff = 'openff-1.2.0'
creation_date = datetime.datetime.now() # TODO : Use Sep  7 02:07
series_name = '2020-09-06-ugi-tBu-x3110-3v3m-2020-04-Jacobs'
description = 'COVID Moonshot Sprint 4 to prioritize Ugi compounds based on x3110 (LON-WEI-adv59df6-2) to optimize substituents in P2 pocket'
sdf_filename = f'sprint-4/{series_name}.sdf'
json_filename = f'sprint-4/{series_name}-newschema.json'
smarts = 'N(C(=O)c2ccco2)C(C(=O)NC)c2cccnc2' # SMARTS for common core scaffold
xchem_fragment_id = 'x3110'
receptors = f'../receptors/monomer/Mpro-{xchem_fragment_id}_0_bound-protein-thiolate.pdb'
receptor_variant = {'catalytic-dyad' : 'His41(+) Cys145(-)'}
temperature = 300.0 # kelvin
pH = 7.3 # pH (fluorescence assay)

def get_compound_id(microstate_id):
    import re
    match = re.match('^(?P<compound_id>\S+)_(?P<microstate_suffix>\d+)$', microstate_id)
    if match is None:
        # No warts; compound and microstate are identical
        compound_id = microstate_id
    else:
        # Remove the wart
        compound_id = match.group('compound_id')
    return compound_id

def get_moldb(series):
    lig_file_name = ligand_dict[series]
    from openeye import oechem
    moldb = oechem.OEMolDatabase()
    moldb.Open(lig_file_name)
    return moldb

# Project pair
from fah_xchem.schema import ProjectPair
fah_projects = ProjectPair(
    complex_phase=13426,
    solvent_phase=13427
)

# Compound series metadata
from fah_xchem.schema import CompoundSeriesMetadata
series_metadata = CompoundSeriesMetadata(
    name=series_name,
    description=description,
    creator=creator,
    created_at=creation_date,
    xchem_project=xchem_project,
    receptor_variant=receptor_variant,
    temperature_kelvin=temperature,
    ionic_strength_millimolar=70.0,
    pH=pH,
    fah_projects=fah_projects
)

# Compounds and microstates
microstates = list()
compounds = dict()
from fah_xchem.schema import Microstate, CompoundMetadata, Compound, CompoundMicrostate
from openeye import oechem
moldb = oechem.OEMolDatabase()
moldb.Open(sdf_filename)
nmols = moldb.NumMols()
oemol = oechem.OEGraphMol()
smiles_flag = oechem.OESMILESFlag_Canonical | oechem.OESMILESFlag_ISOMERIC
for molecule_index in track(range(nmols), description='Processing molecules...'):
    moldb.GetMolecule(oemol, molecule_index)
    microstate_id = oemol.GetTitle()
    smiles = oechem.OECreateSmiString(oemol, smiles_flag)
    # Determine if our molecule has warts
    compound_id = get_compound_id(microstate_id)
    # Compile information about the microstate
    microstate = Microstate(microstate_id=microstate_id, smiles=smiles)
    microstates.append(microstate)
    # Extract experimental data, if present
    experimental_data = dict()
    if oechem.OEHasSDData(oemol,'f_avg_pIC50'):
        experimental_data['pIC50'] = oechem.OEGetSDData(oemol, 'f_avg_pIC50')
    # Extract information about the compound
    compound_metadata = CompoundMetadata(
        compound_id=oemol.GetTitle(),
        smiles=oechem.OECreateSmiString(oemol, smiles_flag),
        experimental_data=experimental_data,
    )
    # Add microstate to compound if it already exists
    compound_microstates = list() # previous compound microstates
    if compound_id in compounds:
        compound_microstates = compounds[compound_id].microstates
    # Create new (version of) compound
    compound = Compound(
        metadata=compound_metadata,
        microstates=compound_microstates + [microstate]
    )
    # Store compound
    compounds[compound_id] = compound


# Add transformations
from fah_xchem.schema import Transformation, CompoundMicrostate
reference_molecule_index = 0
run_id = 0
transformations = list()
for molecule_index in track(range(1, nmols), description='Creating transformations...'):
    transformation = Transformation(
        run_id=run_id,
        xchem_fragment_id=xchem_fragment_id,
        initial_microstate=CompoundMicrostate(
            compound_id=get_compound_id(microstates[molecule_index].microstate_id),
            microstate_id=microstates[molecule_index].microstate_id
        ),
        final_microstate=CompoundMicrostate(
            compound_id=get_compound_id(microstates[reference_molecule_index].microstate_id),
            microstate_id=microstates[reference_molecule_index].microstate_id
        )
    )
    transformations.append(transformation)
    run_id += 1

# Compile compound series
from fah_xchem.schema import CompoundSeries
compound_series = CompoundSeries(
    metadata=series_metadata,
    compounds=list(compounds.values()),
    transformations=transformations
)

print(f'Writing JSON to {json_filename}')
with open(json_filename, "w") as f:
    f.write(compound_series.json())

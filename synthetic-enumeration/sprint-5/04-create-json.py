#!/bin/env python

"""
Create fah-xchem compound set JSON for Sprint 5 file using new schema
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
ff = 'openff-1.2.0'
creation_date = datetime.datetime.now() # TODO : Use Sep  7 02:07

#xchem_fragment_id = 'x11498'
#reference_compound_id = 'MAT-POS-b3e365b9-1'

xchem_fragment_id = 'x12073'
reference_compound_id  = 'MAT-POS-8a69d52e-7'

series_name = f'sprint-5-{xchem_fragment_id}-monomer-neutral'
description = f"COVID Moonshot Sprint 5 to prioritize benzopyran-isoquinoline series based on {xchem_fragment_id} ({reference_compound_id}) to optimize substituents in the P1' pocket with Mpro monomer and neutral Cys145:His41"
microstates_sdf_filename = f'docked/sprint-5-microstates-{xchem_fragment_id}-sorted.sdf' # microstates with docked geometries
compounds_sdf_filename = f'docked/sprint-5-compounds.sdf' # compounds with annotation
json_filename = f'json/{series_name}.json' # output filename
smarts = 'C(=O)Nc1cncc2ccccc12' # SMARTS for common core scaffold : linker:isoquinoline
receptors = f'../receptors/monomer/Mpro-{xchem_fragment_id}_0_bound-protein-thiolate.pdb'
receptor_variant = {'catalytic-dyad' : 'His41(0) Cys145(0)'}
temperature = 300.0 # kelvin
pH = 7.3 # pH (fluorescence assay)
ionic_strength_millimolar = 70.0 # millimolar
reference_microstate_id = f'{reference_compound_id}_1' # microstate id for reference for transformations

def get_compound_id(microstate_id):
    """
    Extract the compound ID from a microstate ID (which includes a wart suffix like '_1', '_2')

    Parameters
    ----------
    microstate_id : str
        The microstate ID, which includes a wart suffix (e.g. 'MAT-POS-8a69d52e-7_1')

    Returns
    -------
    compound_id : str
        The compound ID (e.g. 'MAT-POS-8a69d52e-7_1')
    """
    import re
    match = re.match('^(?P<compound_id>\S+)_(?P<microstate_suffix>\d+)$', microstate_id)
    if match is None:
        # No warts; compound and microstate are identical
        compound_id = microstate_id
    else:
        # Remove the wart
        compound_id = match.group('compound_id')
    return compound_id

# Project pair
from fah_xchem.schema import ProjectPair
fah_projects = ProjectPair(
    complex_phase=13428,
    solvent_phase=13429
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
    ionic_strength_millimolar=ionic_strength_millimolar,
    pH=pH,
    fah_projects=fah_projects
)

# Compounds
from fah_xchem.schema import Compound, CompoundMetadata
smiles_flag = oechem.OESMILESFlag_Canonical | oechem.OESMILESFlag_ISOMERIC

from openeye import oechem
compounds_moldb = oechem.OEMolDatabase()
compounds_moldb.Open(compounds_sdf_filename)
compounds = dict()
ncompounds = compounds_moldb.NumMols()
oemol = oechem.OEGraphMol()
for compound_index in track(range(ncompounds), description='Processing compounds...'):
    # Get compound
    compounds_moldb.GetMolecule(oemol, compound_index)
    # Set ID and SMILES
    compound_id = oemol.GetTitle()
    smiles = oechem.OECreateSmiString(oemol, smiles_flag)
    # Extract experimental data, if present
    experimental_data = dict()
    if oechem.OEHasSDData(oemol,'f_avg_pIC50'):
        pIC50 = oechem.OEGetSDData(oemol, 'f_avg_pIC50')
        if pIC50 != '':
            pIC50 = float(pIC50)
            experimental_data['pIC50'] = pIC50
    # Extract information about the compound
    compound_metadata = CompoundMetadata(
        compound_id=compound_id,
        smiles=oechem.OECreateSmiString(oemol, smiles_flag),
        experimental_data=experimental_data,
    )
    # Create new compound
    compound = Compound(
        metadata=compound_metadata,
        microstates=list()
    )
    # Store compound
    compounds[compound_id] = compound

# Microstates
from fah_xchem.schema import CompoundMicrostate, Microstate
microstates_moldb = oechem.OEMolDatabase()
microstates_moldb.Open(microstates_sdf_filename)
microstates = list()

nmicrostates = microstates_moldb.NumMols()
oemol = oechem.OEGraphMol()
for microstate_index in track(range(nmicrostates), description='Processing microstates...'):
    microstates_moldb.GetMolecule(oemol, microstate_index)
    microstate_id = oemol.GetTitle()
    smiles = oechem.OECreateSmiString(oemol, smiles_flag)
    # Determine if our molecule has warts
    compound_id = oechem.OEGetSDData(oemol, 'compound')
    # Compile information about the microstate
    microstate = Microstate(microstate_id=microstate_id, smiles=smiles)
    microstates.append(microstate)
    # Add microstate to compound if it already exists
    compound_microstates = list() # previous compound microstates
    if compound_id in compounds:
        compound_microstates = compounds[compound_id].microstates
    else:
        raise Exception(f'Microstate {microstate_id} supposedly belongs to compound {compound_id}, but compound not found')
    # Create new (version of) compound
    compound = Compound(
        metadata=compound_metadata,
        microstates=compound_microstates + [microstate]
    )
    # Store compound
    compounds[compound_id] = compound

# Find reference molecule index for transformations
print(f'Identifying reference microstate {reference_microstate_id} for transformations...')
reference_microstate_found = False
for reference_microstate_index, microstate in enumerate(microstates):
    if microstate.microstate_id == reference_microstate_id:
        reference_microstate_found = True
        break
if not reference_microstate_found:
    raise Exception(f'Could not find reference microstate id {reference_microstate_id} among microstates')

# Create transformations
from fah_xchem.schema import Transformation, CompoundMicrostate
run_id = 0
transformations = list()
for microstate_index in track(range(nmicrostates), description='Creating transformations to reference microstate...'):
    # Skip the self-transformation
    if microstate_index == reference_microstate_index:
        continue

    # Create the transformation
    transformation = Transformation(
        run_id=run_id,
        xchem_fragment_id=xchem_fragment_id,
        initial_microstate=CompoundMicrostate(
            compound_id=get_compound_id(microstates[microstate_index].microstate_id),
            microstate_id=microstates[microstate_index].microstate_id
        ),
        final_microstate=CompoundMicrostate(
            compound_id=get_compound_id(microstates[reference_microstate_index].microstate_id),
            microstate_id=microstates[reference_microstate_index].microstate_id
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

# Write JSON
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

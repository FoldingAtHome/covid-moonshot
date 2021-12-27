#!/bin/env python

"""
Create fah-xchem compound set JSON for Sprint 5 file using new schema

DO NOT USE - it uses the compounds SDF 2D file which may screw up stereochemistryx
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
ff = 'openff-2.0.0'
creation_date = datetime.datetime.now()

reference_compound_ids = {
    'P1800_0A' : 'VLA-UCB-50c39ae8-2', # chromane-5spiro-isoquinoline
    'P2113_0B' : 'JOH-MSK-1f2dff76-4', # tetrahydroisoquinoline-5spiro-isoquinoline
    'P2222_0A' : 'JOH-MSK-1f2dff76-4', # tetrahydroisoquinoline-5spiro-isoquinoline
}
prefix = 'sprint-11-2021-12-26'
run_id = 0
description = 'COVID Moonshot Sprint 11 for optimizing spiro compounds'

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

import os
os.makedirs('json', exist_ok=True)

#assembly_states = ['monomer', 'dimer']
#charge_states = ['neutral', 'charged']

rmsd_restraints = ['unrestrained', 'restrained']
assembly_states = ['dimer']
charge_states = ['neutral', 'charged']

for charge_state in charge_states:
    for assembly_state in assembly_states:
        for rmsd_restraint in rmsd_restraints:
            for xchem_fragment_id in reference_compound_ids.keys():
                reference_compound_id = reference_compound_ids[xchem_fragment_id]
                series_name = f'{prefix}-{xchem_fragment_id}-{assembly_state}-{charge_state}-{rmsd_restraint}'
                print(f'Generating JSON for {series_name}...')
                description = f"{description} based on {xchem_fragment_id} using reference compound {reference_compound_id} with Mpro {assembly_state} and {charge_state} Cys145:His41 catalytic dyad"

                json_filename = f'json/{series_name}.json' # output filename
                #microstates_sdf_filename = f'docked/{prefix}-{xchem_fragment_id}-{assembly_state}-{charge_state}.sdf' # microstates with docked geometries
                microstates_sdf_filename = f'docked/{prefix}-{xchem_fragment_id}-{assembly_state}-neutral.sdf' # microstates with docked geometries
                compounds_smi_filename = f'docked/{prefix}-compounds.smi' # compounds with annotation
                smarts = 'a1aaa2a(a1)AAAC23*~*~*(c4cncc5ccccc45)~*3' # SMARTS for common core scaffold : tetralinlike-5spiro-isoquinoline
                suffix = 'charged' if charge_state=='charged' else ''
                receptors = f'../receptors/{assembly_state}/Mpro-{xchem_fragment_id}_bound-protein{suffix}.pdb'
                if charge_state=='neutral':
                    receptor_variant = {'catalytic-dyad' : 'His41(0) Cys145(0)'}
                else:
                    receptor_variant = {'catalytic-dyad' : 'His41(+) Cys145(-)'}
                receptor_variant['rmsd_restraint'] = rmsd_restraint
                temperature = 300.0 # kelvin
                pH = 7.3 # pH (fluorescence assay)
                ionic_strength_millimolar = 70.0 # millimolar
                reference_microstate_id = f'{reference_compound_id}_1' # microstate id for reference for transformations

                # Project pair
                from fah_xchem.schema import ProjectPair
                fah_projects = ProjectPair(
                    complex_phase=13458, # complex
                    solvent_phase=13459 # solvent
                )

                # Compound series metadata
                # Encode dict of assembly and charge states?
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

                from openeye import oechem
                print('Processing compounds...')
                compounds = dict()
                from openeye import oechem
                import numpy as np
                # TODO: Use CSV filename instead of SMI, since we need it to retain experimental data
                with oechem.oemolistream(compounds_smi_filename) as ifs:
                    for oemol in ifs.GetOEGraphMols():
                        # Set ID and SMILES
                        compound_id = oemol.GetTitle()
                        #smiles = Molecule.from_openeye(oemol, allow_undefined_stereo=True).to_smiles()
                        smiles = oechem.OEMolToSmiles(oemol)
                        # We are migrating experimental data to another file
                        experimental_data = dict()
                        # Extract information about the compound
                        compound_metadata = CompoundMetadata(
                            compound_id=compound_id,
                            smiles=smiles,
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
                print('Processing microstates...')
                from fah_xchem.schema import CompoundMicrostate, Microstate
                microstates = list()
                with oechem.oemolistream(microstates_sdf_filename) as ifs:
                    for oemol in ifs.GetOEGraphMols():
                        microstate_id = oemol.GetTitle()
                        smiles = oechem.OEMolToSmiles(oemol)
                        # Identify parent compound
                        compound_id = oechem.OEGetSDData(oemol, 'compound')
                        if not compound_id in compounds:
                            raise Exception(f'Microstate {microstate_id} supposedly belongs to compound {compound_id}, but compound not found')
                        compound = compounds[compound_id]
                        # No experimental data will be present in compound metadata
                        experimental_data = dict()
                        # Rebuild compound_metadata with experimental data
                        compound_metadata = CompoundMetadata(
                            compound_id=compound_id,
                            smiles=compound.metadata.smiles,
                            experimental_data=experimental_data,
                        )
                        # Compile information about the microstate
                        microstate = Microstate(microstate_id=microstate_id, smiles=smiles)
                        microstates.append(microstate)
                        # Add microstate to compound
                        compound = Compound(
                            metadata=compound_metadata,
                            microstates=compound.microstates + [microstate]
                        )
                        # Store compound
                        compounds[compound_id] = compound

                # Find reference molecule index for transformations
                print(f'Identifying reference microstate {reference_microstate_id} for transformations...')
                reference_microstate = None
                for reference_microstate_index, microstate in enumerate(microstates):
                    if microstate.microstate_id == reference_microstate_id:
                        reference_microstate = microstate
                if reference_microstate is None:
                    print(compounds)
                    print(microstates)
                    raise Exception(f'Could not find reference microstate id {reference_microstate_id} among microstates')
                print(f'Reference microstate is ligand index {reference_microstate_index}')
                print(reference_microstate)

                # Create transformations
                from fah_xchem.schema import Transformation, CompoundMicrostate
                transformations = list()
                print('Creating transformations to reference microstate...')
                for microstate_index, microstate in enumerate(microstates):
                    # Skip the self-transformation
                    if microstate_index == reference_microstate_index:
                        continue

                    # Create the transformation
                    transformation = Transformation(
                        run_id=run_id,
                        xchem_fragment_id=xchem_fragment_id,
                        initial_microstate=CompoundMicrostate(
                            compound_id=get_compound_id(reference_microstate.microstate_id),
                            microstate_id=reference_microstate.microstate_id
                        ),
                        final_microstate=CompoundMicrostate(
                            compound_id=get_compound_id(microstate.microstate_id),
                            microstate_id=microstate.microstate_id
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

print(f'There are {run_id} total RUNs')

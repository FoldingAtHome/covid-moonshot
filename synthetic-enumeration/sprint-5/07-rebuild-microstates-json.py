"""
Regenerate JSON from microstates apearing in actual RUNs in case John screwed up numbering

06-extract-transformation-molecules.py must be run first

"""

def annotate_with_assay_data(oemols, assay_data_filename):
    """
    Annotate the set of molecules with activity data using SD tags

    Parameters
    ----------
    oemols : dict of OEMol
        Dict of molecules to annotate
    assay_data_filename
        Filename of CSV file containing activity data
    """
    # Load assay data
    assayed_molecules = dict()
    with oechem.oemolistream(assay_data_filename) as ifs:
        for mol in ifs.GetOEGraphMols():
            smiles = oechem.OEMolToSmiles(mol)
            assayed_molecules[smiles] = oechem.OEGraphMol(mol)
    print(f'Loaded data for {len(assayed_molecules)} assayed molecules')

    # Copy all SDData from assayed molecules that match title
    nmols_with_assay_data = 0
    for smiles in oemols:
        oemol = oemols[smiles]
        if smiles in assayed_molecules:
            assayed_molecule = assayed_molecules[smiles]
            oemol.SetTitle(assayed_molecule.GetTitle())
            oechem.OECopySDData(oemol, assayed_molecule)
            nmols_with_assay_data += 1

    print(f'Found assay data for {nmols_with_assay_data} / {len(oemols)} molecules')

if __name__ == '__main__':
    new_json_filename = 'json/sprint-5-x12073-monomer-neutral-simulated.json'
    from fah_xchem.schema import CompoundSeries, Transformation, CompoundMicrostate, CompoundMetadata, Microstate, Compound

    # Read source JSON
    source_json_filename = 'json/sprint-5-x12073-monomer-neutral.json'
    import json
    with open(source_json_filename, "r") as infile:
        compound_series = CompoundSeries.parse_obj(json.loads(infile.read()))
        
    # Read all microstates
    import os
    project = '13429'
    from glob import glob
    runs = glob(f'{project}/RUNS/RUN*')
    nruns = len(runs)
    from rich.progress import track
    runs = list()
    oemols = dict()
    for run in track(range(nruns), description='Reading initial and final molecules from all RUNs...'):
        try:
            rundata = dict()
            for prefix in ['initial', 'final']:
                filename = os.path.join(f'{project}/RUNS/RUN{run}', f'molecule-{prefix}.csv')
                from openeye import oechem
                with oechem.oemolistream(filename) as ifs:
                    for oemol in ifs.GetOEGraphMols():
                        smiles = oechem.OEMolToSmiles(oemol)
                        #oechem.OESetSDData(oemol, 'RUN', f'RUN{run}')
                        oemols[smiles] = oechem.OEGraphMol(oemol)
                        rundata[f'{prefix}_smiles'] = smiles
            if ('initial_smiles' in rundata) and ('final_smiles' in rundata):
                rundata['run'] = run
                runs.append(rundata)
        except Exception as e:
            print(e)

    # Annotate with activity data if we have it
    annotate_with_assay_data(oemols, 'activity-data/activity-data-2020-11-30.csv')

    # Generate compounds and microstates
    compounds = dict()
    for smiles in track(oemols, description='Processing compounds'):
        oemol = oemols[smiles]
        # Set ID and SMILES
        compound_id = oemol.GetTitle()
        smiles = oechem.OEMolToSmiles(oemol)
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
            smiles=oechem.OEMolToSmiles(oemol),
            experimental_data=experimental_data,
        )
        # Create a single microstate
        microstate = Microstate(microstate_id=compound_id, smiles=smiles)
        # Create new compound
        compound = Compound(
            metadata=compound_metadata,
            microstates=[microstate]
        )
        # Store compound
        compounds[compound_id] = compound

    # Index microstates
    compound_microstates = dict()
    for compound_id in track(compounds,description= 'Indexing microstates...'):
        compound = compounds[compound_id]
        for microstate in compound.microstates:
            compound_microstate = CompoundMicrostate(compound_id=compound.metadata.compound_id, microstate_id=microstate.microstate_id)
            compound_microstates[microstate.smiles] = compound_microstate

    # Regenerate transformations
    xchem_fragment_id = compound_series.transformations[0].xchem_fragment_id
    transformations = list()
    for rundata in track(runs, description=f'Generating transformations from {len(runs)} runs'):
        try:
            transformation = Transformation(
                run_id=rundata['run'],
                xchem_fragment_id=xchem_fragment_id,
                initial_microstate=compound_microstates[rundata['initial_smiles']],
                final_microstate=compound_microstates[rundata['final_smiles']]
            )
            transformations.append(transformation)
        except Exception as e:
            print(rundata)
            print('--------------')

    # Repackage
    compound_series = CompoundSeries(
        metadata=compound_series.metadata,
        compounds=[compound for compound in compounds.values()],
        transformations=transformations,
    )

    # Write
    print(f'Writing new JSON files to {new_json_filename}...')
    with open(new_json_filename, 'wt') as outfile:
        outfile.write(compound_series.json())

"""
Regenerate corrected JSON from actual RUNs in case John screwedup numbering
"""

def extract_transformation(run, compound_microstates, project):
    import perses
    import openmoltools
    import numpy as np
    from openeye import oechem
    from fah_xchem.schema import Transformation

    try:
        npz = np.load(f'{project}/RUNS/RUN{run}/htf.npz', allow_pickle=True)
        x = npz['arr_0']
        htf = x.item()
        old_smiles = oechem.OEMolToSmiles(htf._topology_proposal.old_topology.residue_oemol)
        new_smiles = oechem.OEMolToSmiles(htf._topology_proposal.new_topology.residue_oemol)
        
        if (old_smiles not in compound_microstates):
            print(f'{old_smiles} not found')
            return None
        elif (new_smiles not in compound_microstates):
            print(f'{new_smiles} not found')
            return None

        #print(run, old_smiles, new_smiles, compound_microstates[old_smiles].microstate_id, compound_microstates[new_smiles].microstate_id)
        
        transformation = Transformation(
            run_id=run,
            xchem_fragment_id=xchem_fragment_id,
            initial_microstate=compound_microstates[old_smiles],
            final_microstate=compound_microstates[new_smiles]
        )
    except Exception as e:
        print(e)
        return None

    return transformation

def canonical_smiles(smiles):
    from openeye import oechem
    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, smiles)
    return oechem.OEMolToSmiles(mol)

if __name__ == '__main__':
    from fah_xchem.schema import CompoundSeries, Transformation, CompoundMicrostate

    new_json_filename = 'json/sprint-5-x12073-monomer-neutral-corrected.json'
    from fah_xchem.schema import CompoundSeries, Transformation, CompoundMicrostate

    # Read source JSON
    source_json_filename = 'json/sprint-5-x12073-monomer-neutral.json'
    import json
    with open(source_json_filename, "r") as infile:
        compound_series = CompoundSeries.parse_obj(json.loads(infile.read()))
        
    # Index microstates
    compound_microstates = dict()
    for compound in compound_series.compounds:
        for microstate in compound.microstates:
            compound_microstate = CompoundMicrostate(compound_id=compound.metadata.compound_id, microstate_id=microstate.microstate_id)
            compound_microstates[canonical_smiles(microstate.smiles)] = compound_microstate

    xchem_fragment_id = compound_series.transformations[0].xchem_fragment_id

    # Scan RUNs
    print('Scanning RUNs...')
    project = '13429'
    from glob import glob
    runs = glob(f'{project}/RUNS/RUN*')
    nruns = len(runs)

    def extract_transformation_partial(run):
        return extract_transformation(run, compound_microstates, project)

    print("Reading RUNs...")
    from tqdm import tqdm
    from multiprocessing import Pool
    pool = Pool(64)
    transformations = list()
    for result in tqdm(pool.imap_unordered(func=extract_transformation_partial, iterable=range(nruns)), total=nruns):
        if result is not None:
            transformations.append(result)

    # Sort transformations
    transformations.sort(key=lambda transformation : transformation.run_id)

    # Repackage
    compound_series = CompoundSeries(
        metadata=compound_series.metadata,
        compounds=compound_series.compounds,
        transformations=transformations
    )

    # Write
    print(f'Writing new JSON files to {new_json_filename}...')
    with open(new_json_filename, 'wt') as outfile:
        outfile.write(compound_series.json())

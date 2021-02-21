"""
Extract initial and final molecules for each transformation.
"""

def extract_transformation(run, project):
    """
    Extract initial and final molecules from a transformation and store them within the RUN directory in CSV, SDF, and mol2 formats.

    Parameters
    ----------
    run : str
        The RUN to extract initial and final molecules from (e.g. 'RUN3')
    project : str
        Project path (e.g. '13429')   
    """
    import perses
    import numpy as np
    from openeye import oechem
    import os

    try:
        # Extract initial and final OEMols
        npz = np.load(f'{project}/RUNS/RUN{run}/htf.npz', allow_pickle=True)
        x = npz['arr_0']
        htf = x.item()
        oemols = dict()
        oemols['initial'] = htf._topology_proposal.old_topology.residue_oemol
        oemols['final'] = htf._topology_proposal.new_topology.residue_oemol

        # Write oemols
        for prefix in ['initial', 'final']:
            for suffix in ['csv', 'smi', 'sdf', 'mol2']:
                filename = os.path.join(f'{project}/RUNS/RUN{run}', f'molecule-{prefix}.{suffix}')
                with oechem.oemolostream(filename) as ofs:
                    oechem.OEWriteMolecule(ofs, oemols[prefix])
    except Exception as e:
        print(e)
        
    return None

if __name__ == '__main__':
    # Processing RUNs
    print('Scanning RUNs...')
    project = '13429'
    from glob import glob
    runs = glob(f'{project}/RUNS/RUN*')
    nruns = len(runs)

    def extract_transformation_partial(run):
        return extract_transformation(run, project)

    print("Processing RUNs...")
    from tqdm import tqdm
    from multiprocessing import Pool
    pool = Pool(64)
    transformations = list()
    for result in tqdm(pool.imap_unordered(func=extract_transformation_partial, iterable=range(nruns)), total=nruns):
        pass

"""
Generate JSON caches of parameters for all ligands
"""

small_molecule_forcefield = 'openff-1.1.0'
basepath = 'json-files'

def generate_parameters(molecule, basepath='json-files', small_molecule_forcefield='openff-1.1.0'):
    """
    Generate JSON parameter cache for a molecule in f'{basepath}/{molecule.name}.json'

    Parameters
    ----------
    molecule : openforcefield.topology.Molecule
        The molecule to parameterize

    """
    # Create generator
    import os
    cache_filename = f'parallel/{molecule.name}.json'
    if os.path.exists:
        return

    # Generate and cache parameters
    from openmmforcefields.generators import SystemGenerator
    system_generator = SystemGenerator(small_molecule_forcefield=small_molecule_forcefield, molecules=[molecule], cache=cache_filename)
    try:
        system_generator.create_system(molecule.to_topology().to_openmm())
    except Exception as e:
        print(f'FAILED: {molecule.smiles}')
        print(e)

    del system_generator

if __name__ == '__main__':

    # Read molecules
    from openforcefield.topology import Molecule
    molecules = Molecule.from_file('../docking/covid_submissions_03_26_2020 - docked.sdf', allow_undefined_stereo=True)

    # Extract unique molecules
    molecules = list(set(molecules))
    print(f'There are {len(molecules)} unique molecules...')

    import os
    os.mkdir(basepath)

    # Parameterize in parallel
    from multiprocessing import Pool
    from tqdm import tqdm
    with Pool() as pool:
        max_ = len(molecules)
        with tqdm(total=max_) as pbar:
            for i, _ in enumerate(pool.imap_unordered(generate_parameters, molecules)):
                pbar.update()

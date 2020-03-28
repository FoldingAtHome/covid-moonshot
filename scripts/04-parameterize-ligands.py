"""
Generate JSON cache of parameters for all ligands
"""

cache_filename = 'covid_submissions_03_26_2020-openff-1.1.0.json'
small_molecule_forcefield = 'openff-1.1.0'

def generate_parameters(molecule):
    # Create generator
    from openmmforcefields.generators import SystemGenerator
    system_generator = SystemGenerator(molecules=[molecule], cache=cache_filename)
    try:
        system_generator.create_system(molecule.to_topology().to_openmm())
    except Exception as e:
        with open('failures.smi','a') as outfile:
            outfile.write(molecule.to_smiles() + '\n')

if __name__ == '__main__':

    # Read molecules
    from openforcefield.topology import Molecule
    molecules = Molecule.from_file('../docking/covid_submissions_03_26_2020 - docked.sdf', allow_undefined_stereo=True)

    # Extract unique molecules
    molecules = list(set(molecules))

    # Parameterize
    from multiprocessing import Pool
    from tqdm import tqdm
    with Pool() as pool:
        max_ = len(molecules)
        with tqdm(total=max_) as pbar:
            for i, _ in enumerate(pool.imap_unordered(generate_parameters, molecules)):
                pbar.update()


    # Check JSON file
    print('Checking JSON file integrity...')
    from openmmforcefields.generators import SystemGenerator
    system_generator = SystemGenerator(molecules=[molecule], cache=cache_filename)
    for molecule in molecules:
        try:
            system_generator.create_system(molecule.to_topology().to_openmm())
        except Exception as e:
            print(e)

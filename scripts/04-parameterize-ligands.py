"""
Generate JSON cache of parameters for all ligands
"""

small_molecule_forcefield = 'openff-1.1.0'

def generate_parameters(molecule):
    # Create generator
    from openmmforcefields.generators import SystemGenerator
    import os
    if not os.path.exists('parallel'):
        os.mkdir('parallel')
    cache_filename = f'{molecule.name}.json'
    system_generator = SystemGenerator(small_molecule_forcefield=small_molecule_forcefield, molecules=[molecule], cache=cache_filename)
    try:
        system_generator.create_system(molecule.to_topology().to_openmm())
    except Exception as e:
        with open('failures.smi','a') as outfile:
            outfile.write(molecule.to_smiles() + '\n')
    del system_generator

if __name__ == '__main__':

    # Read molecules
    from openforcefield.topology import Molecule
    molecules = Molecule.from_file('../docking/covid_submissions_03_26_2020 - docked.sdf', allow_undefined_stereo=True)

    # Extract unique molecules
    unique_smiles = set()
    unique_molecules = list()
    for molecule in molecules:
        if molecule.to_smiles() not in unique_smiles:
            unique_smiles.add(molecule.to_smiles())
            unique_molecules.append(molecule)
    molecules = unique_molecules

    # Parameterize
    from multiprocessing import Pool
    from tqdm import tqdm
    with Pool() as pool:
        max_ = len(molecules)
        with tqdm(total=max_) as pbar:
            for i, _ in enumerate(pool.imap_unordered(generate_parameters, molecules)):
                pbar.update()

    # Merge JSON files
    print('Merging JSON files...')
    from tinydb import TinyDB
    cache_filename = 'covid_submissions_03_26_2020-openff-1.1.0.json'
    tinydb_kwargs = { 'sort_keys' : True, 'indent' : 4, 'separators' : (',', ': ') } # for pretty-printing
    with TinyDB(cache_filename, **tinydb_kwargs) as db:
        table = db.table(small_molecule_forcefield)
        for molecule in molecules:
            with TinyDB(f'parallel/{molecule.name}.json', **tinydb_kwargs) as moldb:
                moltable = moldb.table(small_molecule_forcefield)
                for record in moltable:
                    table.insert(record)

    # Check JSON file
    print('Checking JSON file integrity...')
    from openmmforcefields.generators import SystemGenerator
    system_generator = SystemGenerator(small_molecule_forcefield=small_molecule_forcefield, molecules=molecules, cache=cache_filename)
    for molecule in molecules:
        try:
            system_generator.create_system(molecule.to_topology().to_openmm())
        except Exception as e:
            print(e)

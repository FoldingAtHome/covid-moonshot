"""
Generate JSON cache of parameters for all ligands
"""

def generate_parameters(molecule):
    # Create generator
    small_molecule_forcefield = 'openff-1.1.0'
    from openmmforcefields.generators import SystemGenerator
    system_generator = SystemGenerator(molecules=[molecule], cache='covid_submissions_03_26_2020-openff-1.1.0.json')
    system_generator.create_system(molecule.to_topology().to_openmm())

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

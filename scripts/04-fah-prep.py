"""
Prepare for FAH simulation
"""

# Read ligands
def read_csv_molecules(filename):
    """Read molecules from the specified path

    Parameters
    ----------
    filename : str
        File from which molecules are to be read

    Returns
    -------
    molecules : list of openeye.oechem.OEMol
        The read molecules
    """

    from openeye import oechem
    mol = oechem.OEMol()
    molecules = list()
    with oechem.oemolistream(filename) as ifs:
        while oechem.OEReadCSVFile(ifs, mol):
            molecules.append(oechem.OEMol(mol))
    return molecules

def score(molecule, field='Hybrid2'):
    """Return the docking score"""
    from openeye import oechem
    value = oechem.OEGetSDData(molecule, field)
    return float(value)

if __name__ == '__main__':
    # Parse arguments
    import argparse

    parser = argparse.ArgumentParser(description='Aggregate results and coalesce Folding@Home RUNs.')
    parser.add_argument('--docked', dest='docked_basedir', type=str, default='covid_submissions_03_31_2020',
                        help='base directory for docked molecules (default: covid_submissions_03_31_2020/)')
    parser.add_argument('--output', dest='output_filename', type=str, default='covid_submissions_03_31_2020-fah.csv',
                        help='output aggregated CSV file (default: covid_submissions_03_31_2020-fah.csv)')
    parser.add_argument('--clean', dest='clean', action='store_true', default=False,
                        help='if specified, will only store minimal information for each molecule (default: False)')

    args = parser.parse_args()

    import os
    from openeye import oechem
    from tqdm import tqdm
    from glob import glob

    # Read the docked molecules as CSV
    docked_filenames = glob(f'{args.docked_basedir}/docking/*ligand.sdf')
    docked_molecules = list()
    molecule = oechem.OEGraphMol()
    for filename in tqdm(docked_filenames):
        with oechem.oemolistream(filename) as ifs:
            molecules = list()
            while oechem.OEReadMolecule(ifs, molecule):
                molecules.append( molecule.CreateCopy() )
        docked_molecules += molecules
    print(f'{len(docked_molecules)} read')

    # Sort molecules
    docked_molecules.sort(key=score)

    gromacs_basedir = os.path.join(args.docked_basedir, 'gromacs')

    # Write molecules and set up directories in sorted order to CSV
    run_index = 0
    print(f'Writing sorted molecules to {args.output_filename}')
    with oechem.oemolostream(args.output_filename) as ofs:
        for molecule in tqdm(docked_molecules):
            # Check if gro/top exists
            files_missing = False
            for phase in ['complex', 'ligand']:
                for ext in ['gro', 'top']:
                    filename = os.path.join(gromacs_basedir, f'{molecule.GetTitle()} - {phase}.{ext}')
                    if not os.path.exists(filename):
                        files_missing = True
            if files_missing:
                continue
            
            # Add RUN number
            oechem.OESetSDData(molecule, 'run', f'RUN{run_index}')

            if args.clean:
                for sdpair in oechem.OEGetSDDataPairs(molecule):
                    if sdpair.GetTag() not in ['Hybrid2', 'docked_fragment', 'fragments', 'site', 'run']:
                        oechem.OEDeleteSDData(molecule, sdpair.GetTag())

            # Copy files
            run_dir = os.path.join(args.docked_basedir, 'fah', f'RUN{run_index}')
            os.makedirs(run_dir, exist_ok=True)
            import shutil
            for phase in ['complex', 'ligand']:
                for ext in ['gro', 'top']:
                    src = os.path.join(gromacs_basedir, f'{molecule.GetTitle()} - {phase}.{ext}')
                    dst = os.path.join(run_dir, f'{phase}.{ext}')
                    shutil.copyfile(src, dst)

            oechem.OEWriteMolecule(ofs, molecule)
            
            run_index += 1

"""
Combine all docking results
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

    parser = argparse.ArgumentParser(description='Dock a molecule and (optionally) prepare it for alchemical free energy calculations.')
    parser.add_argument('--docked', dest='docked_basedir', type=str, default='docked',
                        help='base directory for docked molecules (default: docked/)')
    parser.add_argument('--output', dest='output_filename', type=str, default='docked-aggregated.csv',
                        help='output aggregated CSV file (default: docked-aggregated.csv)')
    parser.add_argument('--clean', dest='clean', action='store_true', default=False,
                        help='if specified, will only store minimal information for each molecule (default: False)')

    args = parser.parse_args()

    import os
    from openeye import oechem
    from tqdm import tqdm
    from glob import glob

    # Read the docked molecules as CSV
    docked_filenames = glob(f'{args.docked_basedir}/*ligand.sdf')
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

    # Write molecules in sorted order to CSV
    print(f'Writing sorted molecules to {args.output_filename}')
    with oechem.oemolostream(args.output_filename) as ofs:
        for molecule in tqdm(docked_molecules):
            if args.clean:
                for sdpair in oechem.OEGetSDDataPairs(molecule):
                    if sdpair.GetTag() not in ['Hybrid2', 'fragments', 'site']:
                        oechem.OEDeleteSDData(molecule, sdpair.GetTag())
            oechem.OEWriteMolecule(ofs, molecule)
    

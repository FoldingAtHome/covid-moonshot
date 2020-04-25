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

def read_original_smiles(filename):
    """Read mapping of original SMILES for each compound name

    Parameters
    ----------
    filename : str
        File from which molecules are to be read

    Returns
    -------
    original_smiles : dict
        original_smiles[compound_id] is the original SMILES
    """
    original_smiles = dict()
    with open(filename, 'r') as infile:
        for line in infile:
            elements = line.split(',')
            smiles = elements[0]
            cid = elements[1]
            original_smiles[cid] = smiles

    return original_smiles

def score(molecule, field='Hybrid2'):
    """Return the docking score"""
    from openeye import oechem
    value = oechem.OEGetSDData(molecule, field)
    return float(value)

if __name__ == '__main__':
    # Parse arguments
    import argparse

    parser = argparse.ArgumentParser(description='Aggregate docking results')
    parser.add_argument('--molecules', dest='molecules_filename', type=str, default='covid_submissions_all_info-2020-04-23.csv',
                        help='molecules CSV file to pull from (default: covid_submissions_all_info-2020-04-23.csv)')
    parser.add_argument('--docked', dest='docked_basedir', type=str, default='covid_submissions_all_info-2020-04-23-docked',
                        help='base directory for docked molecules (default: covid_submissions_all_info-2020-04-23-docked/)')
    parser.add_argument('--output', dest='output_filename', type=str, default='covid_submissions_all_info-2020-04-23-docked.csv',
                        help='output aggregated CSV file (default: covid_submissions_all_info-2020-04-23.csv)')
    parser.add_argument('--clean', dest='clean', action='store_true', default=False,
                        help='if specified, will only store minimal information for each molecule (default: False)')
    parser.add_argument('--fragalysis', dest='fragalysis', action='store_true', default=False,
                        help='if specified, format data for fragalysis (default: False)')

    args = parser.parse_args()

    import os
    from openeye import oechem
    from tqdm import tqdm
    from glob import glob

    # Read the docked molecules from SDF
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

    # Read original SMILES to annotate molecules
    if args.molecules_filename:
        original_smiles = read_original_smiles(args.molecules_filename)
    
    # Write molecules in sorted order to CSV
    # TODO: Translate into final format
    print(f'Writing sorted molecules to {args.output_filename}')
    with oechem.oemolostream(args.output_filename) as ofs:
        for molecule in tqdm(docked_molecules):
            if args.clean or args.fragalysis:
                for sdpair in oechem.OEGetSDDataPairs(molecule):
                    tag = sdpair.GetTag()
                    value = sdpair.GetValue()
                    if tag not in ['Hybrid2', 'docked_fragment', 'fragments', 'covalent_distance_min', 'covalent_distance_mean', 'covalent_distance_stderr']:
                        oechem.OEDeleteSDData(molecule, sdpair.GetTag())
                    # Translate into fragalysis
                    if args.fragalysis:
                        if tag == 'Hybrid2':
                            oechem.OESetSDData(molecule, 'Chemgauss4', value)
                            oechem.OEDeleteSDData(molecule, tag)
                        elif tag == 'fragments':
                            oechem.OESetSDData(molecule, 'ref_mols', value + '_0')
                            oechem.OEDeleteSDData(molecule, tag)
                        elif tag == 'docked_fragment':
                            oechem.OESetSDData(molecule, 'ref_pdb', value + '_0')
                            oechem.OEDeleteSDData(molecule, tag)
                        elif tag == 'covalent_distance_min':
                            oechem.OESetSDData(molecule, 'CYS145-warhead min dist (A)', value)
                            oechem.OEDeleteSDData(molecule, tag)
                        elif tag == 'covalent_distance_mean':
                            oechem.OESetSDData(molecule, 'CYS145-warhead avg dist (A)', value)
                            oechem.OEDeleteSDData(molecule, tag)
                        elif tag == 'covalent_distance_stderr':
                            oechem.OESetSDData(molecule, 'CYS145-warhead avg dist stderr (A)', value)
                            oechem.OEDeleteSDData(molecule, tag)
                        # Add original SMILES
                        oechem.OESetSDData(molecule, 'original SMILES', original_smiles[molecule.GetTitle()])

            # Write the molecule
            oechem.OEWriteMolecule(ofs, molecule)



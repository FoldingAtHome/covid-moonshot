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

if __name__ == '__main__':
    # Dock the ligands
    import os
    from openeye import oechem

    prefix = 'covid_submissions_03_26_2020'

    # Generate list of all X-ray fragments
    fragment_molecules = read_csv_molecules(os.path.join('../molecules', 'mpro_fragments_03_25_2020.csv'))
    all_fragments = [ oechem.OEGetSDData(molecule, "fragments") for molecule in fragment_molecules ]

    # Filter fragments to only retain those with prepped receptor structures
    all_fragments = [ fragment for fragment in all_fragments if os.path.exists(os.path.join(f'../receptors/Mpro-{fragment}-receptor.oeb.gz')) ]

    def score(molecule):
        """Return the docking score"""
        return float(oechem.OEGetSDData(molecule, 'Chemgauss4 Score'))

    # Read molecules, keeping only best scoring ones
    from tqdm import tqdm
    print('Reading records...')
    docked_molecules = dict()
    receptors = dict()
    import drconvert
    for fragment in tqdm(all_fragments):
        filename = f'../docking/{prefix} - docked to {fragment}.oeb'
        if os.path.exists(filename):
            datarecords = drconvert.RecordConvertToMols(filename)
            for molecule in datarecords:
                # Annotate which fragment this was docked to
                oechem.OESetSDData(molecule, 'fragments', fragment)

                if molecule.NumAtoms() > 1000:
                    # Store protein
                    molecule.SetTitle(fragment)
                    receptors[fragment] = oechem.OEMol(molecule)
                else:
                    # Store ligands in best docked poses
                    oechem.OEDeleteSDData(molecule, 'Number of Confs')
                    CID = molecule.GetTitle()
                    if CID not in docked_molecules:
                        docked_molecules[CID] = oechem.OEMol(molecule)
                    else:
                        if score(molecule) < score(docked_molecules[CID]):
                            docked_molecules[CID] = oechem.OEMol(molecule)

    # Sort compounds
    docked_molecules = [docked_molecules[CID] for CID in docked_molecules]
    docked_molecules.sort(key=score)

    # Write all molecules in order of increasing score
    print('Writing molecules to CSV...')
    filename = f'{prefix} - top docked.csv'
    with oechem.oemolostream(filename) as ofs:
        ofs.SetFormat(oechem.OEFormat_CSV)
        for molecule in docked_molecules:
            if molecule.NumAtoms() < 1000:
                # Only write small molecules (not proteins)
                oechem.OEWriteMolecule(ofs, molecule)

    # Write molecules
    print('Writing molecules to SDF...')
    filename = f'{prefix} - top docked.sdf'
    with oechem.oemolostream(filename) as ofs:
        for molecule in docked_molecules:
            if molecule.NumAtoms() < 1000:
                # Only write small molecules (not proteins)
                oechem.OEWriteMolecule(ofs, molecule)

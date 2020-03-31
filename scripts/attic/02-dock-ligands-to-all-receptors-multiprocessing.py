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

# Read molecules to dock
import os
prefix = 'covid_submissions_03_26_2020'
#prefix = 'mpro_fragments_03_25_2020'
molecules = read_csv_molecules(os.path.join('../molecules', prefix + '.csv'))


def dock_molecules_to_receptor(receptor_filename):
    """
    Dock the specified molecules, writing out to specified file

    Parameters
    ----------
    receptor_filename : str
        Receptor .oeb.gz filename
    fragment : str
        The fragment name to dock to

    """
    import os

    # Read the receptor
    from openeye import oechem, oedocking
    receptor = oechem.OEGraphMol()
    if not oedocking.OEReadReceptorFile(receptor, receptor_filename):
        oechem.OEThrow.Fatal("Unable to read receptor")

    if not oedocking.OEReceptorHasBoundLigand(receptor):
        raise Exception("Receptor does not have bound ligand")

    #print('Initializing receptor...')
    dockMethod = oedocking.OEDockMethod_Hybrid2
    dockResolution = oedocking.OESearchResolution_Default
    dock = oedocking.OEDock(dockMethod, dockResolution)
    success = dock.Initialize(receptor)

    # Set up Omega
    #print('Expanding conformers...')
    from openeye import oeomega
    #omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Dense)
    omegaOpts = oeomega.OEOmegaOptions()
    omega = oeomega.OEOmega(omegaOpts)
    omega.SetStrictStereo(False)

    # Dock molecules
    docked_molecules = list()
    for mol in molecules:
        dockedMol = oechem.OEGraphMol()

        # Expand conformers
        omega.Build(mol)

        # Dock molecule
        #print(f'Docking {mol.NumConfs()} conformers...')
        retCode = dock.DockMultiConformerMolecule(dockedMol, mol)
        if (retCode != oedocking.OEDockingReturnCode_Success):
            print("Docking Failed with error code " + oedocking.OEDockingReturnCodeGetName(retCode))
            continue

        # Store docking data
        sdtag = oedocking.OEDockMethodGetName(dockMethod)
        oedocking.OESetSDScore(dockedMol, dock, sdtag)
        dock.AnnotatePose(dockedMol)

        docked_molecules.append( oechem.OEMol(dockedMol) )

    return docked_molecules

def dock_molecules_to_file(fragment):
    import oechem
    basedir = 'parallel'

    # Determine reecptor filename
    receptor_filename = os.path.join(f'../receptors/Mpro-{fragment}-receptor.oeb.gz')

    # Dock molecules
    try:
        docked_molecules = dock_molecules_to_receptor(receptor_filename)
    except Exception as e:
        print(e)
        return

    # Write molecules
    if len(docked_molecules) > 0:
        output_filename = os.path.join(basedir, f'{fragment} - docked.oedb')
        with oechem.oemolostream(output_filename) as ofs:
            for docked_molecule in docked_molecules:
                oechem.OESetSDData(docked_molecule, "fragments", fragment)
                oechem.OEWriteMolecule(ofs, docked_molecule)

if __name__ == '__main__':
    # Dock the ligands
    import os
    from openeye import oechem

    # Generate list of all X-ray fragments
    fragment_molecules = read_csv_molecules(os.path.join('../molecules', 'mpro_fragments_03_25_2020.csv'))
    all_fragments = [ oechem.OEGetSDData(molecule, "fragments") for molecule in fragment_molecules ]

    # filter fragments
    all_fragments = [ fragment for fragment in all_fragments if os.path.exists(os.path.join(f'../receptors/Mpro-{fragment}-receptor.oeb.gz')) ]

    # Dock molecules in parallel
    print('Processing in parallel...')
    from multiprocessing import Pool
    from tqdm import tqdm
    with Pool() as pool:
        max_ = len(all_fragments)
        with tqdm(total=max_) as pbar:
            for i, _ in enumerate(pool.imap_unordered(dock_molecules_to_file, all_fragments)):
                pbar.update()

    # Join
    print('Joining files...')
    with oechem.oemolostream(f'{prefix} - docked.sdf') as ofs:
        for molecule in molecules:
            filename = os.path.join(basedir, f'{molecule.GetTitle()} - all-docked.oeb')
            if os.path.exists(filename):
                with oechem.oemolistream(filename) as ifs:
                    mol = oechem.OEMol()
                    while oechem.OEReadMolecule(ifs, mol):
                        oechem.OEWriteMolecule(ofs, mol)

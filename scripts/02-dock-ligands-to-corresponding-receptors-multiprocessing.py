

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

def dock_molecule(molecule, default_receptor='x0387'):
    """
    Dock the specified molecules, writing out to specified file

    Parameters
    ----------
    molecule : OEMol
        The molecule to dock
    default_receptor : str, optional, default='0387'
        The default receptor to dock to

    Returns
    -------
    all_docked_molecules : list of OEMol
        All docked molecules
    """
    import os
    import oechem

    # Make a copy of the molecule
    molecule = oechem.OEMol(molecule)

    # Extract list of corresponding receptor(s)
    fragments = list()    
    fragments = oechem.OEGetSDData(molecule, "fragments").split(',')
    fragments = [ fragment for fragment in fragments if os.path.exists(f'../receptors/Mpro-{fragment}-receptor.oeb.gz') ]

    if len(fragments) == 0:
        fragments = [default_receptor]

    # Dock them
    all_docked_molecules = list()
    for fragment in fragments:
        molecule_to_dock = oechem.OEMol(molecule)

        import os
        receptor_filename = os.path.join(f'../receptors/Mpro-{fragment}-receptor.oeb.gz')
        oechem.OESetSDData(molecule_to_dock, "fragments", fragment)

        # Enumerate reasonable protomers/tautomers
        from openeye import oequacpac
        protomer = oechem.OEMol()
        protomers = [ oechem.OEMol(protomer) for protomer in oequacpac.OEGetReasonableProtomers(molecule_to_dock) ]
        docked_molecules = dock_molecules_to_receptor(receptor_filename, protomers)
        all_docked_molecules += docked_molecules

    return all_docked_molecules

def dock_molecules_to_receptor(receptor_filename, molecules):
    """
    Dock the specified molecules, writing out to specified file

    Parameters
    ----------
    receptor_filename : str
        Receptor .oeb.gz filename
    molecules : list of openeye.oechem.OEMol
        The read molecules to dock

    Returns
    -------
    docked_molecules : list of OEMol
        All docked molecules

    """
    # Read the receptor
    from openeye import oechem, oedocking
    receptor = oechem.OEGraphMol()
    if not oedocking.OEReadReceptorFile(receptor, receptor_filename):
        oechem.OEThrow.Fatal("Unable to read receptor")

    if not oedocking.OEReceptorHasBoundLigand(receptor):
        raise Exception("Receptor does not have bound ligand")

    #print('Initializing receptor...')
    dockMethod = oedocking.OEDockMethod_Hybrid2
    dockResolution = oedocking.OESearchResolution_High
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

def dock_molecule_from_file(name):
    import oechem
    basedir = 'parallel'

    # Read input molecule
    input_filename = os.path.join(basedir, f'{name}.oeb')
    molecule = oechem.OEMol()
    with oechem.oemolistream(input_filename) as ifs:
        oechem.OEReadMolecule(ifs, molecule)

    # Dock it
    try:
        docked_molecules = dock_molecule(molecule)
    except Exception as e:
        print(e)
        return

    # Write molecule
    if len(docked_molecules) > 0:
        output_filename = os.path.join(basedir, f'{name} - docked.oeb')
        with oechem.oemolostream(output_filename) as ofs:
            for docked_molecule in docked_molecules:
                oechem.OEWriteMolecule(ofs, docked_molecule)

if __name__ == '__main__':
    # Dock the ligands
    import os
    from openeye import oechem

    prefix = 'covid_submissions_03_26_2020'
    molecules = read_csv_molecules(os.path.join('../molecules', prefix + '.csv'))

    # Write molecules independently
    print('Splitting molecules...')
    import os
    basedir = 'parallel'
    if not os.path.exists(basedir):
        os.mkdir(basedir)
    for molecule in molecules:
        filename = os.path.join(basedir, f'{molecule.GetTitle()}.oeb')
        with oechem.oemolostream(filename) as ofs:
            oechem.OEWriteMolecule(ofs, molecule)

    # Dock molecules in parallel
    print('Processing in parallel...')
    from multiprocessing import Pool
    from tqdm import tqdm
    with Pool() as pool:
        max_ = len(molecules)
        with tqdm(total=max_) as pbar:
            molecule_names = [ molecule.GetTitle() for molecule in molecules ]
            for i, _ in enumerate(pool.imap_unordered(dock_molecule_from_file, molecule_names)):
                pbar.update()

    # Join
    print('Joining files...')
    with oechem.oemolostream(f'{prefix} - docked.sdf') as ofs:
        for molecule in molecules:
            filename = os.path.join(basedir, f'{molecule.GetTitle()} - docked.oeb')
            if os.path.exists(filename):
                with oechem.oemolistream(filename) as ifs:
                    mol = oechem.OEMol()
                    while oechem.OEReadMolecule(ifs, mol):
                        oechem.OEWriteMolecule(ofs, mol)

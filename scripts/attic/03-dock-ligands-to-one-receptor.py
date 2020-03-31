molecule_sets = [
    'melatonin-and-metabolites',
    'clinical-melatonin-receptor-agonists',
    'enamine-melatonin-realspace-analogues',
    'covid_submissions_03_26_2020',
]

# Read ligands
def read_molecules(filename):
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
        while oechem.OEReadMolecule(ifs, mol):
            molecules.append(oechem.OEMol(mol))
    return molecules

def dock_molecules(receptor_filename, molecules, filename):
    """
    Dock the specified molecules, writing out to specified file

    Parameters
    ----------
    receptor_filename : str
        Receptor .oeb.gz filename
    molecules : list of openeye.oechem.OEMol
        The read molecules to dock
    filename : str
        The filename to stream docked molecules to

    """
    # Read the receptor
    print('Loading receptor...')
    from openeye import oechem, oedocking
    receptor = oechem.OEGraphMol()
    if not oedocking.OEReadReceptorFile(receptor, 'receptor.oeb.gz'):
        oechem.OEThrow.Fatal("Unable to read receptor")

    if oedocking.OEReceptorHasBoundLigand(receptor):
        print("Receptor has a bound ligand")
    else:
        print("Receptor does not have bound ligand")

    print('Initializing receptor...')
    dockMethod = oedocking.OEDockMethod_Hybrid2
    dockResolution = oedocking.OESearchResolution_High
    dock = oedocking.OEDock(dockMethod, dockResolution)
    success = dock.Initialize(receptor)

    # Set up Omega
    from openeye import oeomega
    omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Dense)
    #omegaOpts = oeomega.OEOmegaOptions()
    omega = oeomega.OEOmega(omegaOpts)
    omega.SetStrictStereo(False)

    # Dock molecules
    with oechem.oemolostream(filename) as ofs:
        from tqdm import tqdm
        for mol in tqdm(molecules):
            dockedMol = oechem.OEGraphMol()

            # Expand conformers
            omega.Build(mol)

            # Dock molecule
            retCode = dock.DockMultiConformerMolecule(dockedMol, mol)
            if (retCode != oedocking.OEDockingReturnCode_Success):
                oechem.OEThrow.Fatal("Docking Failed with error code " + oedocking.OEDockingReturnCodeGetName(retCode))

            # Store docking data
            sdtag = oedocking.OEDockMethodGetName(dockMethod)
            oedocking.OESetSDScore(dockedMol, dock, sdtag)
            dock.AnnotatePose(dockedMol)

            # Write molecule
            oechem.OEWriteMolecule(ofs, dockedMol)

# Dock the ligands
for prefix in molecule_sets:
    print(f'Docking {prefix}')
    import os
    molecules = read_molecules(os.path.join('../molecules', prefix + '.csv'))
    dock_molecules('receptor-0104.oeb.gz', molecules, prefix + '-docked.sdf')

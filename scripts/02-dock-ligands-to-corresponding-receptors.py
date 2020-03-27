

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

def dock_molecule(molecule, ofs, default_receptor='x0387'):
    """
    Dock the specified molecules, writing out to specified file

    Parameters
    ----------
    molecule : OEMol
        The molecule to dock
    ofs : oechem.oemolostream
        The filename to stream docked molecules to
    default_receptor : str, optional, default='0387'
        The default receptor to dock to
    """
    import os
    # Extract list of corresponding receptor(s)
    import oechem
    if oechem.OEHasSDData(molecule, "fragments"):
        fragments = oechem.OEGetSDData(molecule, "fragments").split(',')
        fragments = [ fragment for fragment in fragments if os.path.exists(f'../receptors/Mpro-{fragment}-receptor.oeb.gz') ]
        if len(fragments) == 0:
            fragments = [default_receptor]
        for fragment in fragments:
            import os
            receptor_filename = os.path.join(f'../receptors/Mpro-{fragment}-receptor.oeb.gz')
            oechem.OESetSDData(molecule, "fragments", fragment)

            # Enumerate reasonable protomers/tautomers
            from openeye import oequacpac
            protomer = oechem.OEMol()
            protomers = [ oechem.OEMol(protomer) for protomer in oequacpac.OEGetReasonableProtomers(molecule) ]
            dock_molecules_to_receptor(receptor_filename, protomers, ofs)

def dock_molecules_to_receptor(receptor_filename, molecules, ofs):
    """
    Dock the specified molecules, writing out to specified file

    Parameters
    ----------
    receptor_filename : str
        Receptor .oeb.gz filename
    molecules : list of openeye.oechem.OEMol
        The read molecules to dock
    ofs : oechem.oemolostream
        The filename to stream docked molecules to

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
    for mol in molecules:
        dockedMol = oechem.OEGraphMol()

        # Expand conformers
        omega.Build(mol)

        # Dock molecule
        #print(f'Docking {mol.NumConfs()} conformers...')
        retCode = dock.DockMultiConformerMolecule(dockedMol, mol)
        if (retCode != oedocking.OEDockingReturnCode_Success):
            print("Docking Failed with error code " + oedocking.OEDockingReturnCodeGetName(retCode))
            return 

        # Store docking data
        sdtag = oedocking.OEDockMethodGetName(dockMethod)
        oedocking.OESetSDScore(dockedMol, dock, sdtag)
        dock.AnnotatePose(dockedMol)

        # Write molecule
        oechem.OEWriteMolecule(ofs, dockedMol)

# Dock the ligands
prefix = 'covid_submissions_03_26_2020'
import os
from tqdm import tqdm
molecules = read_csv_molecules(os.path.join('../molecules', prefix + '.csv'))
import oechem
with oechem.oemolostream(f'{prefix} - docked.sdf') as ofs:
    for molecule in tqdm(molecules):
        try:
            dock_molecule(molecule, ofs)
        except Exception as e:
            print(e)

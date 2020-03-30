"""
Dock specified ligand to all DiamondMX Mpro structures and prepare for alchemical free energy calculations

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

def dock_molecule_to_receptor(molecule, receptor_filename):
    """
    Dock the specified molecules, writing out to specified file

    Parameters
    ----------
    molecule : oechem.OEMol
        The molecule to dock
    receptor_filename : str
        Receptor .oeb.gz filename

    Returns
    -------
    docked_molecule : openeye.oechem.OEMol
        Returns the best tautomer/protomer in docked geometry, annotated with docking score
        None is returned if no viable docked pose found
    """
    import os

    # Extract the fragment name for the receptor
    fragment = extract_fragment_from_filename(receptor_filename)

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

    # Enumerate tautomers
    from openeye import oequacpac
    tautomer_options = oequacpac.OETautomerOptions()
    tautomer_options.SetMaxTautomersGenerated(4096)
    tautomer_options.SetMaxTautomersToReturn(16)
    tautomer_options.SetCarbonHybridization(True)
    tautomer_options.SetMaxZoneSize(50)
    tautomer_options.SetApplyWarts(True)
    pKa_norm = True
    tautomers = [ oechem.OEMol(tautomer) for tautomer in oequacpac.OEGetReasonableTautomers(molecule, tautomer_options, pKa_norm) ]

    # Set up Omega
    #print('Expanding conformers...')
    from openeye import oeomega
    #omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Dense)
    omegaOpts = oeomega.OEOmegaOptions()
    omega = oeomega.OEOmega(omegaOpts)
    omega.SetStrictStereo(False)

    # Dock tautomers
    docked_molecules = list()
    from tqdm import tqdm
    for mol in tautomers:
        dockedMol = oechem.OEGraphMol()

        # Expand conformers
        omega.Build(mol)

        # Dock molecule
        retCode = dock.DockMultiConformerMolecule(dockedMol, mol)
        if (retCode != oedocking.OEDockingReturnCode_Success):
            #print("Docking Failed with error code " + oedocking.OEDockingReturnCodeGetName(retCode))
            continue

        # Store docking data
        sdtag = oedocking.OEDockMethodGetName(dockMethod)
        oedocking.OESetSDScore(dockedMol, dock, sdtag)
        oechem.OESetSDData(dockedMol, "fragments", fragment)
        dock.AnnotatePose(dockedMol)

        docked_molecules.append( dockedMol.CreateCopy() )

    if len(docked_molecules) == 0:
        return None

    # Select the best-ranked molecule and pose
    # Note that this ignores protonation state and tautomer penalties
    docked_molecules.sort(key=score)
    best_molecule = docked_molecules[0]

    return best_molecule

def extract_fragment_from_filename(filename):
    """Extract the fragment name (e.g. 'x0104') from a filename

    Parameters
    ----------
    filename : str
        The filename ('/path/to/file/Mpro-{fragment}-receptor.oeb.gz')

    Returns
    -------
    fragment : str
        The fragment name ('x####')
    """

    import re
    match = re.search('Mpro-(?P<fragment>x\d+)-receptor.oeb.gz', filename)
    fragment = match.group('fragment')
    return fragment

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
    # Parse arguments
    import argparse

    parser = argparse.ArgumentParser(description='Dock a molecule and (optionally) prepare it for alchemical free energy calculations.')
    parser.add_argument('--molecules', dest='molecules_filename', type=str, default='../molecules/covid_submissions_03_26_2020.csv',
                        help='molecules CSV file to pull from (default: ../molecules/covid_submissions_03_26_2020.csv)')
    parser.add_argument('--index', dest='molecule_index', type=int, required=True,
                        help='index of molecule to dock (0 indexed)')
    parser.add_argument('--receptors', dest='receptor_basedir', type=str, default='../receptors',
                        help='directory of receptor conformations (default: ../receptors)')
    parser.add_argument('--output', dest='output_basedir', type=str, default='.',
                        help='base directory for produced output (default: .)')
    parser.add_argument('--simulate', dest='simulate', action='store_true', default=False,
                        help='prepare for simulation in OpenMM and gromacs (default: False)')

    args = parser.parse_args()

    # Extract molecule
    molecules = read_csv_molecules(args.molecules_filename)
    print(f'{len(molecules)} molecules read')
    if not ((0 <= args.molecule_index) and (args.molecule_index < len(molecules))):
        raise Exception(f'--index <index> must be between 0 and {len(molecules)} for {args.molecules_filename}')
    molecule = molecules[args.molecule_index]

    # Replace title if there is none
    if molecule.GetTitle() == '':
        import os
        head, tail = os.path.split(args.molecules_filename)
        prefix, ext = os.path.splitext(tail)
        molecule.SetTitle(f'{prefix}-{args.molecule_index}')

    # Generate list of all X-ray fragments with prepared receptor structures
    from glob import glob
    import os
    receptor_filenames = glob(os.path.join(args.receptor_basedir, 'Mpro-*-receptor.oeb.gz'))

    # Dock molecules in parallel
    print(f'Docking {molecule.GetTitle()} to all receptors...')
    from tqdm import tqdm
    docked_molecules = list()
    for receptor_filename in tqdm(receptor_filenames):
        docked_molecule = dock_molecule_to_receptor(molecule, receptor_filename)
        if docked_molecule is not None:
            docked_molecules.append(docked_molecule)

    # Extract top pose
    docked_molecules.sort(key=score)
    docked_molecule = docked_molecules[0]

    # Write out top pose
    print('Writing out best pose...')
    import os
    from openeye import oechem, oedocking
    os.makedirs(args.output_basedir, exist_ok=True)
    # Write molecule as CSV
    output_filename = os.path.join(args.output_basedir, f'{molecule.GetTitle()} - docked.csv')
    with oechem.oemolostream(output_filename) as ofs:
        oechem.OEWriteMolecule(ofs, docked_molecule)
    # Write molecule as SDF
    output_filename = os.path.join(args.output_basedir, f'{molecule.GetTitle()} - ligand.sdf')
    with oechem.oemolostream(output_filename) as ofs:
        oechem.OEWriteMolecule(ofs, docked_molecule)
    # Write molecule as mol2
    output_filename = os.path.join(args.output_basedir, f'{molecule.GetTitle()} - ligand.mol2')
    with oechem.oemolostream(output_filename) as ofs:
        oechem.OEWriteMolecule(ofs, docked_molecule)
    # Read receptor
    fragment = oechem.OEGetSDData(docked_molecule, 'fragments')
    receptor_filename = os.path.join(args.receptor_basedir, f'Mpro-{fragment}-receptor.oeb.gz')
    from openeye import oechem, oedocking
    receptor = oechem.OEGraphMol()
    if not oedocking.OEReadReceptorFile(receptor, receptor_filename):
        oechem.OEThrow.Fatal("Unable to read receptor")
    # Write receptor
    output_filename = os.path.join(args.output_basedir, f'{molecule.GetTitle()} - protein.pdb')
    with oechem.oemolostream(output_filename) as ofs:
        oechem.OEWriteMolecule(ofs, receptor)
    # Write joined PDB with ligand and receptor
    output_filename = os.path.join(args.output_basedir, f'{molecule.GetTitle()} - complex.pdb')
    with oechem.oemolostream(output_filename) as ofs:
        oechem.OEWriteMolecule(ofs, receptor)
        oechem.OEWriteMolecule(ofs, docked_molecule)

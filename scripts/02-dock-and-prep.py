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
        Receptor to dock to

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

def prepare_simulation(molecule, output_directory):
    """
    Prepare simulation systems

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
       The molecule to set up
    output_directory : str
       The output directory
    """
    # Parameters
    from simtk import unit, openmm
    water_model = 'tip3p'
    solvent_padding = 10.0 * unit.angstrom
    box_size = openmm.vec3.Vec3(3.4,3.4,3.4)*unit.nanometers
    ionic_strength = 100 * unit.millimolar # 100
    pressure = 1.0 * unit.atmospheres
    collision_rate = 91.0 / unit.picoseconds
    temperature = 300.0 * unit.kelvin
    timestep = 4.0 * unit.femtoseconds
    nsteps_per_iteration = 250
    iterations = 100

    protein_forcefield = 'amber14/protein.ff14SB.xml'
    small_molecule_forcefield = 'openff-1.1.0'
    #small_molecule_forcefield = 'gaff-2.11' # only if you really like atomtypes
    solvation_forcefield = 'amber14/tip3p.xml'

    # Create SystemGenerators
    import os
    from simtk.openmm import app
    from openforcefield.topology import Molecule
    off_molecule = Molecule.from_openeye(molecule, allow_undefined_stereo=True)
    print(off_molecule)
    barostat = openmm.MonteCarloBarostat(pressure, temperature)
    cache = os.path.join(output_directory, f'{molecule.GetTitle()}.json')
    common_kwargs = {'removeCMMotion': False, 'ewaldErrorTolerance': 5e-04,
                     'nonbondedMethod': app.PME, 'hydrogenMass': 3.0*unit.amu}    
    unconstrained_kwargs = {'constraints': None, 'rigidWater': False}
    constrained_kwargs = {'constraints': app.HBonds, 'rigidWater': True}
    forcefields = [protein_forcefield, solvation_forcefield]
    from openmmforcefields.generators import SystemGenerator
    parmed_system_generator = SystemGenerator(forcefields=forcefields, 
                                              molecules=[off_molecule], small_molecule_forcefield=small_molecule_forcefield, cache=cache,
                                              barostat=barostat, 
                                              forcefield_kwargs={**common_kwargs, **unconstrained_kwargs})
    openmm_system_generator = SystemGenerator(forcefields=forcefields, 
                                              molecules=[off_molecule], small_molecule_forcefield=small_molecule_forcefield, cache=cache,
                                              barostat=barostat, 
                                              forcefield_kwargs={**common_kwargs, **constrained_kwargs})

    # Prepare phases
    import os
    print(f'Setting up simulation for {molecule.GetTitle()}...')
    for phase in ['ligand', 'complex']:
        phase_name = f'{molecule.GetTitle()} - {phase}'
        phase_prefix = os.path.join(output_directory, phase_name)
        print(phase_name)

        if os.path.exists(phase_prefix+'.gro') and os.path.exists(phase_prefix+'.top'):
            continue
        
        # Read the unsolvated system into an OpenMM Topology
        pdbfile = app.PDBFile(phase_prefix+'.pdb')
        topology, positions = pdbfile.topology, pdbfile.positions

        # Add solvent
        print('Adding solvent...')
        modeller = app.Modeller(topology, positions)
        if phase == 'ligand':
            kwargs = {'boxSize' : box_size}
        else:
            kwargs = {'padding' : solvent_padding}
        modeller.addSolvent(openmm_system_generator.forcefield, model='tip3p', ionicStrength=ionic_strength, **kwargs)

        # Create an OpenMM system
        system = openmm_system_generator.create_system(modeller.topology)
        
        # Create OpenM Context
        platform = openmm.Platform.getPlatformByName('CUDA')
        platform.setPropertyDefaultValue('Precision', 'mixed')
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(modeller.positions)

        # Report initial potential energy
        state = context.getState(getEnergy=True)
        print(f'{molecule.GetTitle()} {phase} : Initial potential energy is {state.getPotentialEnergy()/unit.kilocalories_per_mole:.3f} kcal/mol')

        # Minimize
        print('Minimizing...')
        openmm.LocalEnergyMinimizer.minimize(context)

        # Equilibrate
        print('Equilibrating...')
        from tqdm import tqdm
        for iteration in tqdm(range(iterations)):
            integrator.step(nsteps_per_iteration)

        # Retrieve state
        state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
        system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
        modeller.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
        print(f'{molecule.GetTitle()} {phase} : Final potential energy is {state.getPotentialEnergy()/unit.kilocalories_per_mole:.3f} kcal/mol')

        # Save as OpenMM
        print('Saving as OpenMM...')
        import gzip
        with gzip.open(phase_prefix+'.integrator.xml.gz', 'wt') as f:
            f.write(openmm.XmlSerializer.serialize(integrator))
        with gzip.open(phase_prefix+'.state.xml.gz','wt') as f:
            f.write(openmm.XmlSerializer.serialize(state))
        with gzip.open(phase_prefix+'.system.xml.gz','wt') as f:
            f.write(openmm.XmlSerializer.serialize(system))
        
        # Convert to gromacs via ParmEd
        print('Saving as gromacs...')
        import parmed
        parmed_system = parmed_system_generator.create_system(modeller.topology)
        #parmed_system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
        structure = parmed.openmm.load_topology(modeller.topology, parmed_system, xyz=state.getPositions(asNumpy=True))
        structure.save(phase_prefix+'.gro', overwrite=True)
        structure.save(phase_prefix+'.top', overwrite=True)

def ensemble_dock(molecule, fragments_to_dock_to):
    """Perform ensemble docking on all fragment X-ray structures
    
    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The molecule to dock
    fragments_to_dock_to : list of str
        List of fragments to dock to

    """
    import os
    from tqdm import tqdm

    # Dock molecule to all receptors
    print(f'Docking {molecule.GetTitle()} to {len(fragments_to_dock_to)} fragment structures...')
    docked_molecules = list()
    for fragment in tqdm(fragments_to_dock_to):
        receptor_filename = os.path.join(args.receptor_basedir, f'Mpro-{fragment}-receptor.oeb.gz')
        if not os.path.exists(receptor_filename):
            # Skip receptors that are not set up
            continue
        docked_molecule = dock_molecule_to_receptor(molecule, receptor_filename)
        if docked_molecule is not None:
            docked_molecules.append(docked_molecule)

    # Extract top pose
    docked_molecules.sort(key=score)
    docked_molecule = docked_molecules[0].CreateCopy()

    # Fix title
    docked_molecule.SetTitle(molecule.GetTitle())

    # Populate scores for all complexes
    from openeye import oechem
    for score_molecule in docked_molecules:
        fragment = oechem.OEGetSDData(score_molecule, 'fragments')
        oechem.OESetSDData(docked_molecule, f'Mpro-{fragment}_dock', str(score(score_molecule)))

    # Populate site info
    fragment = oechem.OEGetSDData(docked_molecule, 'fragments')    
    if fragment in active_site_fragments:
        oechem.OESetSDData(docked_molecule, 'site', 'active-noncovalent')
    elif fragment in covalent_active_site_fragments:
        oechem.OESetSDData(docked_molecule, 'site', 'active-covalent')
    elif fragment in dimer_interface_fragments:
        oechem.OESetSDData(docked_molecule, 'site', 'dimer-interface')

    return docked_molecule

def set_serial(molecule, chainid, first_serial):
    import oechem
    if not oechem.OEHasResidues(molecule):
        oechem.OEPerceiveResidues(molecule, oechem.OEPreserveResInfo_All)
    for atom in molecule.GetAtoms():
        residue = oechem.OEAtomGetResidue(atom)
        residue.SetExtChainID(chainid)
        serial_number = residue.GetSerialNumber()
        residue.SetSerialNumber(serial_number + first_serial)

def transfer_data(molecule, source_directory):
    cmd = f'rsync -avz --include="*/" --include="{molecule.GetTitle()}*" --exclude="*" {source_directory} tug27224@owlsnest.hpc.temple.edu:work'
    import subprocess
    output = subprocess.getoutput(cmd)
    print(output)

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
    parser.add_argument('--output', dest='output_basedir', type=str, default='docked',
                        help='base directory for produced output (default: docked/)')
    parser.add_argument('--simulate', dest='simulate', action='store_true', default=False,
                        help='prepare for simulation in OpenMM and gromacs (default: False)')
    parser.add_argument('--transfer', dest='transfer', action='store_true', default=False,
                        help='transfer simulation data (default: False)')

    args = parser.parse_args()

    active_site_fragments = ['x0072', 'x0104', 'x0107', 'x0161', 'x0195', 'x0305', 'x0354', 'x0387', 'x0395', 'x0397', 'x0426', 'x0434', 'x0540',
        'x0678', 'x0874', 'x0946', 'x0967', 'x0991', 'x0995', 'x1077', 'x1093', 'x1249']
    covalent_active_site_fragments = ['x0689', 'x0691', 'x0692', 'x0705','x0708', 'x0731', 'x0734', 'x0736', 'x0749', 'x0752', 'x0755', 'x0759',
    'x0769', 'x0770', 'x0771', 'x0774', 'x0786', 'x0820', 'x0830', 'x0831', 'x0978', 'x0981', 'x1308', 'x1311', 'x1334', 'x1336', 'x1348',
    'x1351', 'x1358', 'x1374', 'x1375', 'x1380', 'x1382', 'x1384', 'x1385', 'x1386', 'x1392', 'x1402', 'x1412', 'x1418', 'x1425', 'x1458',
    'x1478', 'x1493']
    dimer_interface_fragments = ['x0887', 'x1187']
    fragments_to_dock_to = active_site_fragments + covalent_active_site_fragments + dimer_interface_fragments

    # Make output directory
    import os
    os.makedirs(args.output_basedir, exist_ok=True)

    # Extract molecule
    molecules = read_csv_molecules(args.molecules_filename)
    print(f'{len(molecules)} molecules read')
    if not ((0 <= args.molecule_index) and (args.molecule_index < len(molecules))):
        raise Exception(f'--index <index> must be between 0 and {len(molecules)} for {args.molecules_filename}')
    molecule = molecules[args.molecule_index]

    # Replace title if there is none
    import os
    if molecule.GetTitle() == '':
        head, tail = os.path.split(args.molecules_filename)
        prefix, ext = os.path.splitext(tail)
        molecule.SetTitle(f'{prefix}-{args.molecule_index}')

    # Dock if the molecule has not already been docked
    from openeye import oechem
    sdf_filename = os.path.join(args.output_basedir, f'{molecule.GetTitle()} - ligand.sdf')
    if not os.path.exists(sdf_filename):
        # Dock the molecule
        docked_molecule = ensemble_dock(molecule, fragments_to_dock_to)
    else:
        # Read the molecule
        with oechem.oemolistream(sdf_filename) as ifs:            
            docked_molecule = oechem.OEGraphMol()
            oechem.OEReadMolecule(ifs, docked_molecule)

    import os
    from openeye import oechem, oedocking
        
    # Write molecule as CSV with cleared SD tags
    output_filename = os.path.join(args.output_basedir, f'{molecule.GetTitle()} - docked.csv')
    if not os.path.exists(output_filename):
        docked_molecule_clean = docked_molecule.CreateCopy()
        for sdpair in oechem.OEGetSDDataPairs(docked_molecule_clean):
            if sdpair.GetTag() not in ['Hybrid2', 'fragments', 'site']:
                oechem.OEDeleteSDData(docked_molecule_clean, sdpair.GetTag())
        with oechem.oemolostream(output_filename) as ofs:
            oechem.OEWriteMolecule(ofs, docked_molecule_clean)

    # Write molecule as SDF
    output_filename = os.path.join(args.output_basedir, f'{molecule.GetTitle()} - ligand.sdf')
    if not os.path.exists(output_filename):
        with oechem.oemolostream(output_filename) as ofs:
            oechem.OEWriteMolecule(ofs, docked_molecule)

    # Write molecule as mol2
    output_filename = os.path.join(args.output_basedir, f'{molecule.GetTitle()} - ligand.mol2')
    if not os.path.exists(output_filename):
        with oechem.oemolostream(output_filename) as ofs:
            oechem.OEWriteMolecule(ofs, docked_molecule)

    # Read receptor
    fragment = oechem.OEGetSDData(docked_molecule, 'fragments')
    receptor_filename = os.path.join(args.receptor_basedir, f'Mpro-{fragment}-receptor.oeb.gz')
    print(f'Reading receptor from {receptor_filename}')
    from openeye import oechem, oedocking
    receptor = oechem.OEGraphMol()
    if not oedocking.OEReadReceptorFile(receptor, receptor_filename):
        oechem.OEThrow.Fatal("Unable to read receptor")

    # Set chains
    #set_serial(docked_molecule, 'A', 0)

    # Write receptor
    output_filename = os.path.join(args.output_basedir, f'{molecule.GetTitle()} - protein.pdb')
    if not os.path.exists(output_filename):
        with oechem.oemolostream(output_filename) as ofs:
            oechem.OEWriteMolecule(ofs, receptor)

    # Write joined PDB with ligand and receptor
    output_filename = os.path.join(args.output_basedir, f'{molecule.GetTitle()} - complex.pdb')
    if not os.path.exists(output_filename):
        with oechem.oemolostream(output_filename) as ofs:
            oechem.OEClearResidues(docked_molecule)
            
            #oechem.OEWriteMolecule(ofs, docked_molecule)

            oechem.OETriposAtomNames(docked_molecule)
            oechem.OEWritePDBFile(ofs, docked_molecule, oechem.OEOFlavor_PDB_Default | oechem.OEOFlavor_PDB_BONDS)

            set_serial(receptor, 'X', docked_molecule.NumAtoms()+1)
            #oechem.OEWritePDBFile(ofs, receptor, oechem.OEOFlavor_PDB_Default | oechem.OEOFlavor_PDB_BONDS)
            
            oechem.OEWriteMolecule(ofs, receptor)

    # Write PDB of just ligand
    output_filename = os.path.join(args.output_basedir, f'{molecule.GetTitle()} - ligand.pdb')
    if not os.path.exists(output_filename):
        with oechem.oemolostream(output_filename) as ofs:
            oechem.OEClearResidues(docked_molecule)
            
            #oechem.OEWriteMolecule(ofs, docked_molecule)

            oechem.OETriposAtomNames(docked_molecule)
            oechem.OEWritePDBFile(ofs, docked_molecule, oechem.OEOFlavor_PDB_Default | oechem.OEOFlavor_PDB_BONDS)

    # Prepare simulation
    if args.simulate:
        prepare_simulation(docked_molecule, args.output_basedir)
        
    if args.transfer:
        transfer_data(docked_molecule, args.output_basedir)

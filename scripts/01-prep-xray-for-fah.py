"""
Prepare all X-ray structures for FAH

# Projects

    13430 : apo Mpro monomer His41(0) Cys145(0)
    13431 : apo Mpro monomer His41(+) Cys145(-)
    13432 : holo Mpro monomer His41(0) Cys145(0)
    13433 : holo Mpro monomer His41(+) Cys145(-)
    13434 : apo Mpro dimer His41(0) Cys145(0)
    13435 : apo Mpro dimer His41(+) Cys145(-)
    13436 : holo Mpro dimer His41(0) Cys145(0)
    13437 : holo Mpro dimer His41(+) Cys145(-)

Each RUN corresponds to a different fragment structure

# Manifest

`../structures/metadata.csv` : master index of fragment IDs and RUNs
```
,crystal_name,RealCrystalName,smiles,new_smiles,alternate_name,site_name,pdb_entry
1,Mpro-1q2w-2020-04-Bonanno_0,Mpro-1q2w-2020-04-Bonanno,C[C@H](O)CC(C)(C)O,NA,NA,Mpro-SARS1,1Q2W
3,Mpro-1wof-2020-04-Yang_0,Mpro-1wof-2020-04-Yang,CCOC(O)CC[C@H](C[C@@H]1CCNC1O)N[C@H](O)[C@H](CC(C)C)NC(O)[C@@H](NC(O)[C@H](C)NC(O)C1CC(C)ON1)C(C)C,NA,NA,Mpro-SARS1,1WOF
4,Mpro-2a5i-2020-04-Lee_0,Mpro-2a5i-2020-04-Lee,CCOC(O)[C@@H](O)C[C@@H](O)N(CCC(N)O)NC(O)[C@H](CC1CCCCC1)N[C@H](O)[C@H](CC(C)C)N[C@H](O)OCC1CCCCC1,NA,NA,Mpro-SARS1,2A5I
...
```

First column is used to identify RUN:
* `RUN0` is skipped
* `RUN1` is Mpro-1q2w-2020-04-Bonanno_0
* `RUN2` is skipped
* `RUN3` is Mpro-1wof-2020-04-Yang_0
...

"""

def setup_fah_run(destination_path, protein_pdb_filename, oemol=None, cache=None, restrain_rmsd=False):
    """
    Prepare simulation

    Parameters
    ----------
    destination_path : str
        The path to the RUN to be created
    protein_pdb_filename : str
        Path to protein PDB file
    oemol : openeye.oechem.OEMol, optional, default=None
        The molecule to parameterize, with SDData attached
        If None, don't include the small molecule
    restrain_rmsd : bool, optional, default=False
        If True, restrain RMSD during first equilibration phase
    """
    # Parameters
    from simtk import unit, openmm
    protein_forcefield = 'amber14/protein.ff14SB.xml'
    solvent_forcefield = 'amber14/tip3p.xml'
    small_molecule_forcefield='openff-1.2.0'
    water_model = 'tip3p'
    solvent_padding = 10.0 * unit.angstrom
    ionic_strength = 70 * unit.millimolar # assay buffer: 20 mM HEPES pH 7.3, 1 mM TCEP, 50 mM NaCl, 0.01% Tween-20, 10% glycerol
    pressure = 1.0 * unit.atmospheres
    collision_rate = 1.0 / unit.picoseconds
    temperature = 300.0 * unit.kelvin
    timestep = 4.0 * unit.femtoseconds
    iterations = 1000 # 1 ns equilibration
    nsteps_per_iteration = 250

    # Prepare phases
    import os
    system_xml_filename = os.path.join( destination_path, 'system.xml.bz2')
    integrator_xml_filename = os.path.join(destination_path, 'integrator.xml.bz2')
    state_xml_filename = os.path.join(destination_path, 'state.xml.bz2')

    # Check if we can skip setup
    openmm_files_exist = os.path.exists(system_xml_filename) and os.path.exists(state_xml_filename) and os.path.exists(integrator_xml_filename)
    if openmm_files_exist:
        return

    # Create barostat
    barostat = openmm.MonteCarloBarostat(pressure, temperature)

    # Create RUN directory if it does not yet exist
    os.makedirs(destination_path, exist_ok=True)

    # Load any molecule(s)
    molecule = None
    if oemol is not None:
        from openforcefield.topology import Molecule
        molecule = Molecule.from_openeye(oemol, allow_undefined_stereo=True)
        molecule.name = 'MOL' # Ensure residue is MOL
        print([res for res in molecule.to_topology().to_openmm().residues()])

    # Create SystemGenerator
    import os
    from simtk.openmm import app
    forcefield_kwargs = {'removeCMMotion': False, 'hydrogenMass': 3.0*unit.amu, 'constraints': app.HBonds, 'rigidWater': True}
    periodic_kwargs = {'nonbondedMethod': app.PME, 'ewaldErrorTolerance': 2.5e-04}
    forcefields = [protein_forcefield, solvent_forcefield]
    from openmmforcefields.generators import SystemGenerator
    openmm_system_generator = SystemGenerator(
                                forcefields=forcefields,
                                molecules=molecule, small_molecule_forcefield=small_molecule_forcefield, cache=cache,
                                barostat=barostat,
                                forcefield_kwargs=forcefield_kwargs, periodic_forcefield_kwargs=periodic_kwargs)

    # Read protein
    print(f'Reading protein from {protein_pdb_filename}...')
    pdbfile = app.PDBFile(protein_pdb_filename)
    modeller = app.Modeller(pdbfile.topology, pdbfile.positions)

    if oemol is not None:
        # Add small molecule to the system
        modeller.add(molecule.to_topology().to_openmm(), molecule.conformers[0])
        # DEBUG : Check residue name
        with open(os.path.join(destination_path, 'initial-complex.pdb'), 'wt') as outfile:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

    # Add solvent
    print('Adding solvent...')
    kwargs = {'padding' : solvent_padding}
    modeller.addSolvent(openmm_system_generator.forcefield, model='tip3p', ionicStrength=ionic_strength, **kwargs)

    # Create an OpenMM system
    print('Creating OpenMM system...')
    system = openmm_system_generator.create_system(modeller.topology)

    # Add a virtual bond between protein and ligand to make sure they are not imaged separately
    if oemol is not None:
        import mdtraj as md
        mdtop = md.Topology.from_openmm(modeller.topology) # excludes solvent and ions
        for res in mdtop.residues:
            print(res)
        protein_atom_indices = mdtop.select('(protein and name CA)') # protein CA atoms
        ligand_atom_indices = mdtop.select('((resname MOL) and (mass > 1))') # ligand heavy atoms
        protein_atom_index = int(protein_atom_indices[0])
        ligand_atom_index = int(ligand_atom_indices[0])
        print(f'Creating virtual bond between atoms {protein_atom_index} and {ligand_atom_index}')
        force = openmm.CustomBondForce('0')
        force.addBond(protein_atom_index, ligand_atom_index, [])
        system.addForce(force)

    # Add RMSD restraints if requested
    if restrain_rmsd:
        print('Adding RMSD restraint...')
        kB = unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB 
        kT = kB * temperature     
        import mdtraj as md
        mdtop = md.Topology.from_openmm(pdbfile.topology) # excludes solvent and ions
        #heavy_atom_indices = mdtop.select('mass > 1') # heavy solute atoms
        rmsd_atom_indices = mdtop.select('(protein and (name CA)) or ((resname MOL) and (mass > 1))') # CA atoms and ligand heavy atoms
        rmsd_atom_indices = [ int(index) for index in rmsd_atom_indices ]
        custom_cv_force = openmm.CustomCVForce('(K_RMSD/2)*RMSD^2')
        custom_cv_force.addGlobalParameter('K_RMSD', kT / unit.angstrom**2)
        rmsd_force = openmm.RMSDForce(modeller.positions, rmsd_atom_indices)
        custom_cv_force.addCollectiveVariable('RMSD', rmsd_force)
        force_index = system.addForce(custom_cv_force)

    # Create OpenM Context
    platform = openmm.Platform.getPlatformByName('OpenCL')
    platform.setPropertyDefaultValue('Precision', 'mixed')
    from openmmtools import integrators
    integrator = integrators.LangevinIntegrator(temperature, collision_rate, timestep)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(modeller.positions)

    # Report initial potential energy
    state = context.getState(getEnergy=True)
    print(f'Initial potential energy is {state.getPotentialEnergy()/unit.kilocalories_per_mole:.3f} kcal/mol')

    # Store snapshots in MDTraj trajectory to examine RMSD
    import mdtraj as md
    import numpy as np
    mdtop = md.Topology.from_openmm(pdbfile.topology)
    atom_indices = mdtop.select('all') # all solute atoms
    protein_atom_indices = mdtop.select('protein and (mass > 1)') # heavy solute atoms
    if oemol is not None:
        ligand_atom_indices = mdtop.select('(resname MOL) and (mass > 1)') # ligand heavy atoms
    trajectory = md.Trajectory(np.zeros([iterations+1,len(atom_indices),3], np.float32), mdtop)
    trajectory.xyz[0,:,:] = context.getState(getPositions=True).getPositions(asNumpy=True)[atom_indices] / unit.nanometers

    # Minimize
    print('Minimizing...')
    openmm.LocalEnergyMinimizer.minimize(context)

    # Equilibrate (with RMSD restraint if needed)
    import numpy as np
    from rich.progress import track
    import time
    initial_time = time.time()
    for iteration in track(range(iterations), 'Equilibrating...'):
        integrator.step(nsteps_per_iteration)
        trajectory.xyz[iteration+1,:,:] = context.getState(getPositions=True).getPositions(asNumpy=True)[atom_indices] / unit.nanometers
    elapsed_time = (time.time() - initial_time) * unit.seconds
    ns_per_day = (context.getState().getTime() / elapsed_time) / (unit.nanoseconds / unit.day)
    print(f'Performance: {ns_per_day:8.3f} ns/day')
    
    if restrain_rmsd:
        # Disable RMSD restraint
        context.setParameter('K_RMSD', 0.0)
        
        print('Minimizing...')        
        openmm.LocalEnergyMinimizer.minimize(context)        

        for iteration in track(range(iterations), 'Equilibrating without RMSD restraint...'):
            integrator.step(nsteps_per_iteration)

    # Retrieve state
    state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
    system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    modeller.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
    print(f'Final potential energy is {state.getPotentialEnergy()/unit.kilocalories_per_mole:.3f} kcal/mol')

    # Equilibrate again if we restrained the RMSD
    if restrain_rmsd:
        print('Removing RMSD restraint from system...')
        system.removeForce(force_index)

    #if oemol is not None:
    #    # Check final RMSD
    #    print('checking RMSD...')
    #    trajectory.superpose(trajectory, atom_indices=protein_atom_indices)
    #    protein_rmsd = md.rmsd(trajectory, trajectory[-1], atom_indices=protein_atom_indices)[-1] * 10 # Angstroms
    #    oechem.OESetSDData(oemol, 'equil_protein_rmsd', f'{protein_rmsd:.2f} A')
    #    ligand_rmsd = md.rmsd(trajectory, trajectory[-1], atom_indices=ligand_atom_indices)[-1] * 10 # Angstroms
    #    oechem.OESetSDData(oemol, 'equil_ligand_rmsd', f'{ligand_rmsd:.2f} A')
    #    print('RMSD after equilibration: protein {protein_rmsd:8.2f} A | ligand {ligand_rmsd:8.3f} A')

    # Save as OpenMM
    print('Exporting for OpenMM FAH simulation...')
    import bz2
    with bz2.open(integrator_xml_filename, 'wt') as f:
        f.write(openmm.XmlSerializer.serialize(integrator))
    with bz2.open(state_xml_filename,'wt') as f:
        f.write(openmm.XmlSerializer.serialize(state))
    with bz2.open(system_xml_filename,'wt') as f:
        f.write(openmm.XmlSerializer.serialize(system))
    with bz2.open(os.path.join(destination_path, 'equilibrated-all.pdb.bz2'), 'wt') as f:
        app.PDBFile.writeFile(modeller.topology, state.getPositions(), f)
    with open(os.path.join(destination_path, 'equilibrated-solute.pdb'), 'wt') as f:
        import mdtraj
        mdtraj_topology = mdtraj.Topology.from_openmm(modeller.topology)
        mdtraj_trajectory = mdtraj.Trajectory([state.getPositions(asNumpy=True) / unit.nanometers], mdtraj_topology)
        selection = mdtraj_topology.select('not water')
        mdtraj_trajectory = mdtraj_trajectory.atom_slice(selection)
        app.PDBFile.writeFile(mdtraj_trajectory.topology.to_openmm(), mdtraj_trajectory.openmm_positions(0), f)
    if oemol is not None:
        # Write molecule as SDF, SMILES, and mol2
        for extension in ['sdf', 'mol2', 'smi', 'csv']:
            filename = os.path.join(destination_path, f'molecule.{extension}')
            with oechem.oemolostream(filename) as ofs:
                oechem.OEWriteMolecule(ofs, oemol)

    # Clean up
    del context, integrator

if __name__ == '__main__':
    # Parse arguments
    import argparse

    parser = argparse.ArgumentParser(description='Prepare the specified RUN for FAH by preparing all X-ray structure variants of a specific fragment')
    parser.add_argument('--receptors', dest='receptors_path', type=str, default='../receptors',
                        help='directory of receptor conformations (default: ../receptors)')
    parser.add_argument('--metadata', dest='metadata_filename', type=str, default='../fah-xray/fah-metadata.csv',
                        help='metadata (default: ../fah-xray/fah-metadata.csv)')
    parser.add_argument('--run', dest='run', type=str, required=True,
                        help='RUN index to prepare (zero-indexes first column contents in Mpro.zip metadata.csv)')
    parser.add_argument('--output', dest='output_path', type=str, default='projects',
                        help='base directory for produced output (default: projects/)')

    args = parser.parse_args()

    # Read DiamondMX/XChem structure medatadata
    import os, csv
    from collections import OrderedDict
    metadata = OrderedDict()
    with open(args.metadata_filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            run = row['RUN']
            metadata[run] = row

    # Extract relevant metadata
    run = f'RUN{args.run}'
    if run not in metadata:
        raise Exception(f'{run} not found in metadata.csv')
    print(f'Preparing {run}')
    metadata = metadata[run]

    # Extract crystal_name
    crystal_name = metadata['crystal_name']

    # Read molecule in the appropriate protonation state
    from openeye import oechem
    oemol = oechem.OEMol()
    molecule_filename = os.path.join(args.receptors_path, 'monomer', f'{crystal_name}_bound-ligand.mol2')
    if not os.path.exists(molecule_filename):
        raise Exception(f'{molecule_filename} does not exist')
    with oechem.oemolistream(molecule_filename) as ifs:
        oechem.OEReadMolecule(ifs, oemol)
    # Rename the molecule
    title = metadata['alternate_name']
    print(f'Setting title to {title}')
    oemol.SetTitle(title)

    # Remove dummy atoms
    for atom in oemol.GetAtoms():
        if atom.GetName().startswith('Du'):
            print('Removing dummy atom.')
            oemol.DeleteAtom(atom)

    # Attach all structure metadata to the molecule
    for key in metadata:
        oechem.OESetSDData(oemol, key, metadata[key])

    # Set up all variants
    # TODO: Generalize this to just use just one protonation state
    # and maybe constant-pH simulations?
    import tempfile
    import traceback, sys
    with tempfile.TemporaryDirectory() as tmpdir:
        cache = os.path.join(tmpdir, 'cache.json')

        def prepare_variant(project, run, crystal_name, biounit, dyad_state, oemol):
            assert biounit in ['monomer', 'dimer']

            try:
                print('')
                print(f'PROJ{project}')

                if dyad_state == 'His41(0) Cys145(0)':
                    protein_pdb_filename = os.path.join(args.receptors_path, biounit, f'{crystal_name}_bound-protein.pdb')
                elif dyad_state == 'His41(+) Cys145(-)':
                    protein_pdb_filename = os.path.join(args.receptors_path, biounit, f'{crystal_name}_bound-protein.pdb')
                else:
                    raise Exception("dyad_state must be one of ['His41(0) Cys145(0)', 'His41(+) Cys145(-)']")
                destination_path = os.path.join(args.output_path, project, 'RUNS', f'RUN{run}')
                setup_fah_run(destination_path, protein_pdb_filename, oemol=oemol, cache=cache)
                print('')
            except Exception as e:
                traceback.print_exc(file=sys.stdout)
                print(e)

        prepare_variant('13430', args.run, crystal_name, 'monomer', 'His41(0) Cys145(0)', None)
        prepare_variant('13431', args.run, crystal_name, 'monomer', 'His41(+) Cys145(-)', None)
        prepare_variant('13432', args.run, crystal_name, 'monomer', 'His41(0) Cys145(0)', oemol)
        prepare_variant('13433', args.run, crystal_name, 'monomer', 'His41(+) Cys145(-)', oemol)
        prepare_variant('13434', args.run, crystal_name, 'dimer',   'His41(0) Cys145(0)', None)
        prepare_variant('13435', args.run, crystal_name, 'dimer',   'His41(+) Cys145(-)', None)
        prepare_variant('13436', args.run, crystal_name, 'dimer',   'His41(0) Cys145(0)', oemol)
        prepare_variant('13437', args.run, crystal_name, 'dimer',   'His41(+) Cys145(-)', oemol)

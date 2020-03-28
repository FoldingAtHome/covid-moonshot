#!/usr/bin/env python

from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField
from simtk import openmm, unit
from simtk.openmm import app
from openmmforcefields.generators import SystemGenerator
import pandas as pd
import parmed
import os, sys, subprocess

water_model = 'tip3p'
solvent_padding = 10.0 * unit.angstrom
box_size = openmm.vec3.Vec3(3.4,3.4,3.4)*unit.nanometers
ionic_strength = 100 * unit.millimolar # 100
pressure = 1.0 * unit.atmospheres
collision_rate = 91.0 / unit.picoseconds
temperature = 300.0 * unit.kelvin
timestep = 2.0 * unit.femtoseconds
nsteps_equil = 50000 # test

protein_forcefield = 'amber14/protein.ff14SB.xml'
small_molecule_forcefield = 'openff-1.1.0'
#small_molecule_forcefield = 'gaff-2.11' # only if you really like atomtypes
solvation_forcefield = 'amber14/tip3p.xml'

dataset = 'MOON'
opennmm_write_cutoff = 1000 # only write XML files for first n ligands (bc they big)

def prepare_RL_system():
    RL_complex_structure = ligand_structure + receptor_structure

    barostat = openmm.MonteCarloBarostat(pressure, temperature)
    forcefield_kwargs = {'removeCMMotion': False, 'ewaldErrorTolerance': 5e-04,
        'nonbondedMethod': app.PME, 'constraints': None, 'rigidWater': False}
    system_generator = SystemGenerator(forcefields=[protein_forcefield,solvation_forcefield],
        barostat=barostat, forcefield_kwargs=forcefield_kwargs, molecules=[ligand],
        small_molecule_forcefield=small_molecule_forcefield)

    modeller = app.Modeller(RL_complex_structure.topology, RL_complex_structure.positions)
    modeller.addSolvent(system_generator.forcefield, model='tip3p',
        padding=solvent_padding, ionicStrength=ionic_strength) # padding=solvent_padding for RL, boxSize=box_size for L

    system = system_generator.create_system(modeller.topology)
    solvated_structure = parmed.openmm.load_topology(modeller.topology,
        system, xyz=modeller.positions)

    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)

    # minimize and equilibrate
    print('Minimizing...')
    platform = openmm.Platform.getPlatformByName('CUDA')
    platform.setPropertyDefaultValue('CudaDeviceIndex', sys.argv[3])
    context = openmm.Context(system, integrator)
    context.setPositions(modeller.positions)
    openmm.LocalEnergyMinimizer.minimize(context)

    print('Equilibrating...')
    integrator.step(nsteps_equil)

    print('Saving RL GMX files...')
    state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
    system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    parmed_system = parmed.openmm.load_topology(modeller.topology,
        system, xyz=state.getPositions(asNumpy=True))
    parmed_system.save(f'{RL_output_prefix}/conf.gro', overwrite=True)
    parmed_system.save(f'{RL_output_prefix}/topol.top', overwrite=True)

    if ligand_ndx <= openmm_write_cutoff:
        with open(f'{RL_output_prefix}/integrator.xml', 'w') as f:
            f.write(openmm.XmlSerializer.serialize(integrator))
        with open(f'{RL_output_prefix}/state.xml','w') as f:
            f.write(openmm.XmlSerializer.serialize(state))
        with open(f'{RL_output_prefix}/system.xml','w') as f:
            f.write(openmm.XmlSerializer.serialize(system))

def prepare_L_system():
    L_complex_structure = ligand_structure

    barostat = openmm.MonteCarloBarostat(pressure, temperature)
    forcefield_kwargs = {'removeCMMotion': False, 'ewaldErrorTolerance': 5e-04,
        'nonbondedMethod': app.PME, 'constraints': None, 'rigidWater': False}
    system_generator = SystemGenerator(forcefields=[protein_forcefield,solvation_forcefield],
        barostat=barostat, forcefield_kwargs=forcefield_kwargs, molecules=[ligand],
        small_molecule_forcefield=small_molecule_forcefield)

    modeller = app.Modeller(L_complex_structure.topology, L_complex_structure.positions)
    modeller.addSolvent(system_generator.forcefield, model='tip3p',
        boxSize=box_size, ionicStrength=ionic_strength) # padding=solvent_padding for RL, boxSize=box_size for L

    system = system_generator.create_system(modeller.topology)
    solvated_structure = parmed.openmm.load_topology(modeller.topology,
        system, xyz=modeller.positions)

    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)

    # minimize and equilibrate
    print('Minimizing...')
    platform = openmm.Platform.getPlatformByName('CUDA')
    platform.setPropertyDefaultValue('CudaDeviceIndex', sys.argv[3])
    context = openmm.Context(system, integrator)
    context.setPositions(modeller.positions)
    openmm.LocalEnergyMinimizer.minimize(context)

    print('Equilibrating...')
    integrator.step(nsteps_equil)

    print('Saving GMX files...')
    state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
    system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    parmed_system = parmed.openmm.load_topology(modeller.topology,
        system, xyz=state.getPositions(asNumpy=True))
    parmed_system.save(f'{L_output_prefix}/conf.gro', overwrite=True)
    parmed_system.save(f'{L_output_prefix}/topol.top', overwrite=True)

    if ligand_ndx <= openmm_write_cutoff:
        with open(f'{L_output_prefix}/integrator.xml', 'w') as f:
            f.write(openmm.XmlSerializer.serialize(integrator))
        with open(f'{L_output_prefix}/state.xml','w') as f:
            f.write(openmm.XmlSerializer.serialize(state))
        with open(f'{L_output_prefix}/system.xml','w') as f:
            f.write(openmm.XmlSerializer.serialize(system))


for ligand_ndx in range(int(sys.argv[1]),int(sys.argv[2])):

    print(f'Processing LIG{ligand_ndx}...')
    receptor_file = f'protein-0387.pdb'
    RL_output_prefix = f'MOON_RL/LIG{ligand_ndx}'
    L_output_prefix = f'MOON_L/LIG{ligand_ndx}'
    ligand_file = f'{L_output_prefix}/LIG{ligand_ndx}.sdf'
    ligand = Molecule.from_file(ligand_file)
    receptor = app.PDBFile(receptor_file)
    receptor_structure = parmed.load_file(receptor_file)
    ligand_structure = parmed.load_file(f'{L_output_prefix}/LIG{ligand_ndx}.pdb')

    # prepare RL system
    if not os.path.exists(f'{RL_output_prefix}/conf.gro'):
        print(f'Processing RL System for LIG{ligand_ndx}')
        try: # catch bad ligands
            prepare_RL_system()
        except Exception as e:
            print(f'Exception occured for RL LIG{ligand_ndx}: {e}')
            continue

    # prepare L system
    if not os.path.exists(f'{L_output_prefix}/conf.gro'):
        print(f'Processing L System for LIG{ligand_ndx}')
        try: # catch bad ligands
            prepare_L_system()
        except Exception as e:
            print(f'Exception occured for L LIG{ligand_ndx}: {e}')
            continue

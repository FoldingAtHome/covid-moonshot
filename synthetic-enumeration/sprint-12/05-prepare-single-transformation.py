import yaml 
import sys
import itertools
import os
import json

yaml_filename = 'sprint-11.yaml'

def run_relative_perturbation(compound_series, transformation_index, microstate_sdf_filename, tidy=True, rmsd_restraint=True):
    """
    Set up a relative transformation
    """
    fah_projects = compound_series.metadata.fah_projects
    transformation = compound_series.transformations[transformation_index]

    from openeye import oechem
    with oechem.oemolistream(microstate_sdf_filename) as ifs:
        microstate_ids = [ mol.GetTitle() for mol in ifs.GetOEGraphMols() ]

    # Decode assembly and charge states from microstate_sdf_filename prefix
    assembly_state = 'monomer' if 'monomer' in compound_series.metadata.description else 'dimer'
    suffix = '-thiolate' if 'charged' in compound_series.metadata.description else ''
    protein_pdb_filename = f'receptors/{assembly_state}/Mpro-{transformation.xchem_fragment_id}_bound-protein{suffix}.pdb'
    print(protein_pdb_filename)

    outputdir = f'RUN{transformation.run_id}'

    print(f'Starting relative calculation')
    print(transformation)
    
    # rewrite yaml file
    with open(yaml_filename, "r") as yaml_file:
        options = yaml.load(yaml_file, Loader=yaml.FullLoader)
    options['old_ligand_index'] = microstate_ids.index(transformation.initial_microstate.microstate_id)
    options['new_ligand_index'] = microstate_ids.index(transformation.final_microstate.microstate_id)
    options['trajectory_directory'] = outputdir
    options['protein_pdb'] = protein_pdb_filename
    options['ligand_file'] = microstate_sdf_filename
    options['use_given_geometries'] = True
    options['rmsd_restraint'] = rmsd_restraint
    options['num_equilibration_iterations'] = 2000 # 2 ns
    options['num_equilibration_steps_per_iteration'] = 250

    # Metadata
    microstate_smiles = { microstate.microstate_id : microstate.smiles for compound in compound_series.compounds for microstate in compound.microstates } 

    options['initial_microstate'] = transformation.initial_microstate.microstate_id
    options['initial_compound'] = transformation.initial_microstate.compound_id
    options['initial_smiles'] = microstate_smiles[transformation.initial_microstate.microstate_id]

    options['final_microstate'] = transformation.final_microstate.microstate_id
    options['final_compound'] = transformation.final_microstate.compound_id
    options['final_smiles'] = microstate_smiles[transformation.final_microstate.microstate_id]

    # Write options YAML file
    print(options)
    for project in fah_projects.solvent_phase, fah_projects.complex_phase:
        yamldir = f'{project}/RUNS/{outputdir}'
        new_yaml = os.path.join(yamldir, f'perses.yaml')
        os.makedirs(yamldir, exist_ok=True)
        with open(new_yaml, 'w') as outfile:
            yaml.dump(options, outfile)
    
    # run the simulation
    os.system(f'perses-fah {new_yaml}')

    print('Relative calculation of ligand {} to {} complete'.format(transformation.initial_microstate.microstate_id, transformation.final_microstate.microstate_id))

    if tidy:
        os.remove(new_yaml)

    return


# work out which ligand pair to run
            
def load_json(filename):
    """
    Load a JSON file that may be .bz2 or .gz compressed

    Parameters
    ----------
    filename : str
        JSON filename

    Returns
    -------
    json : list or dict
        Contents of the JSON file in Python object model

    """
    if '.bz2' in filename:
        import bz2
        with bz2.open(filename, 'rt') as infile:
            return json.load(infile)
    elif '.gz' in filename:
        import gzip
        with gzip.open(filename, 'rt') as infile:
            return json.load(infile)
    else:
        with open(filename, 'rt') as infile:
            return json.load(infile)

# Sprint 
import os
from glob import glob
json_filenames = glob('json/*.json')

# DEBUG: Force order
json_filenames = [
    'json/sprint-11-2021-12-26-P1800_0A-dimer-neutral-unrestrained.json',
    'json/sprint-11-2021-12-26-P1800_0A-dimer-neutral-restrained.json',
    'json/sprint-11-2021-12-26-P1800_0A-dimer-charged-unrestrained.json',
    'json/sprint-11-2021-12-26-P1800_0A-dimer-charged-restrained.json',
]

transformation_index = int(sys.argv[1])

# Select appropriate run index from entire set of JSON filenames
for json_filename in json_filenames:
    print(json_filename)
    head, tail = os.path.split(json_filename)
    prefix, ext = os.path.splitext(tail)
    microstate_sdf_filename = f'docked/{prefix}.sdf' # TODO: Encapsulate this in JSON file as source_sdf_filename and source_molecule_index?
    json_data = load_json(json_filename)
    from fah_xchem.schema import CompoundSeries
    compound_series = CompoundSeries.parse_obj(json_data)
    print(prefix)

    if transformation_index < len(compound_series.transformations):
        break
    else:
        transformation_index -= len(compound_series.transformations)

rmsd_restraint = False
if compound_series.metadata.receptor_variant['rmsd_restraint'] == 'restrained':
    rmsd_restraint = True
        
run_relative_perturbation(compound_series, transformation_index, microstate_sdf_filename, rmsd_restraint=rmsd_restraint)

import yaml 
import sys
import itertools
import os
import json

yaml_filename = 'sprint-7.yaml'

def run_relative_perturbation(compound_series, transformation_index, microstate_sdf_filename, tidy=True, rmsd_restraint=True):
    """
    Set up a relative transformation
    """
    transformation = compound_series.transformations[transformation_index]

    from openeye import oechem
    with oechem.oemolistream(microstate_sdf_filename) as ifs:
        microstate_ids = [ mol.GetTitle() for mol in ifs.GetOEGraphMols() ]

    # Decode assembly and charge states from microstate_sdf_filename prefix
    assembly_state = 'monomer' if 'monomer' in compound_series.metadata.description else 'dimer'
    suffix = '-thiolate' if 'charged' in compound_series.metadata.description else ''
    protein_pdb_filename = f'receptors/{assembly_state}/Mpro-{transformation.xchem_fragment_id}-protein{suffix}.pdb'
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
    options['num_equilibration_steps_per_iteration'] = 500

    # Write options YAML file
    print(options)
    yamldir = 'yaml'
    new_yaml = os.path.join(yamldir, f'fah_{outputdir}.yaml')
    os.makedirs(yamldir, exist_ok=True)
    with open(new_yaml, 'w') as outfile:
        yaml.dump(options, outfile)
    
    # run the simulation
    os.system(f'perses-fah {new_yaml}')

    print('Relative calcluation of ligand {} to {} complete'.format(transformation.initial_microstate.microstate_id, transformation.final_microstate.microstate_id))

    if tidy:
        os.remove(new_yaml)

    return


# work out which ligand pair to run

def load_json(filename):
    """
    Load a JSON file that may be .bz2 or .gz compressed
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

# Sprint 5
import os
from glob import glob
json_filenames = glob('json/*.json')

# DEBUG: Force order
json_filenames = [
    'json/sprint-7-2021-05-11-stereofilter-P0033-dimer-neutral-restrained.json',
    'json/sprint-7-2021-05-11-stereofilter-P0033-dimer-neutral-unrestrained.json',
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

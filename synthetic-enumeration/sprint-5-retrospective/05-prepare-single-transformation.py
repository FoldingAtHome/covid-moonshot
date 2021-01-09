import yaml 
import sys
import itertools
import os
import json

yaml_filename = 'sprint-5.yaml'

def run_relative_perturbation(compound_series, transformation_index, microstate_sdf_filename, tidy=True):
    """
    Set up a relative transformation
    """
    transformation = compound_series.transformations[transformation_index]

    from openeye import oechem
    with oechem.oemolistream(microstate_sdf_filename) as ifs:
        microstate_ids = [ mol.GetTitle() for mol in ifs.GetOEGraphMols() ]

    # TODO : Decode assembly and charge states from microstate_sdf_filename prefix
    protein_pdb_filename = f'../../receptors/monomer/Mpro-{transformation.xchem_fragment_id}_0A_bound-protein.pdb'

    outputdir = f'RUN{transformation.run_id}'

    print(f'Starting relative calcluation')
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
    options['rmsd_restraint'] = True

    # Write options YAML file
    print(options)
    yamldir = 'yaml-monomer'
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
json_filenames = ['json/sprint-5-retrospective-x11498-monomer-neutral.json'] # DEBUG
# TODO : select appropriate  run index from entire set of JSON filenames
for json_filename in json_filenames:
    head, tail = os.path.split(json_filename)
    prefix, ext = os.path.splitext(tail)
    # TODO Fix SDF filename
    prefix = 'sprint-5-retrospective-microstates-x11498-monomer-neutral'
    microstate_sdf_filename = f'docked/{prefix}.sdf' # TODO: Encapsulate this in JSON file as source_sdf_filename and source_molecule_index?
    json_data = load_json(json_filename)
    from fah_xchem.schema import CompoundSeries
    compound_series = CompoundSeries.parse_obj(json_data)

    # Process 0-indexed transformation index
    transformation_index = int(sys.argv[1])
    run_relative_perturbation(compound_series, transformation_index, microstate_sdf_filename)

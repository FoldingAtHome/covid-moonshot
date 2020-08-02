import yaml 
import sys
import itertools
import os
import json

def run_relative_perturbation(run, tidy=True):
    print(f'run details : {run}')
    ligA = run['start']
    ligB = run['end']
    ligandfile = run['ligand']
    proteinfile = run['protein']
    outputdir = run['directory']
    index = run['JOBID']
    ff = run['ff']
    print(f'Starting relative calcluation of ligand {ligA} to {ligB} for {run["target"]} with forcefield {ff}')
    new_yaml = f'fah_{outputdir}.yaml'
    
    # rewrite yaml file
    with open('backtesting.yaml', "r") as yaml_file:
        options = yaml.load(yaml_file, Loader=yaml.FullLoader)
    options['old_ligand_index'] = ligA
    options['new_ligand_index'] = ligB
    options['trajectory_directory'] = outputdir 
    options['protein_pdb'] = f'{proteinfile}'
    options['ligand_file'] = f'{ligandfile}'
    options['small_molecule_forcefield'] = ff 
    with open(new_yaml, 'w') as outfile:
        yaml.dump(options, outfile)
    
    # run the simulation
    os.system(f'perses-fah {new_yaml}')

    print('Relative calcluation of ligand {} to {} complete'.format(ligA,ligB))

    if tidy:
        os.remove(new_yaml)

    return

# work out which ligand pair to run
#series = '2020-07-28.json'
#series = '2020-07-29-retrospective-aminopyridines.json'
series = '2020-08-02-retrospective-aminopyridines-matt.json'
with open(series, 'r') as f:
    data = json.load(f)
this_run = data[sys.argv[1]]

run_relative_perturbation(this_run)

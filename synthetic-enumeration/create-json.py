import numpy as np
import json
import math
import itertools
from openeye import oechem

ff = 'openff-1.2.0'
ligand_dict = {
    # Sprint 1

    #'primary_amine_enumeration' : 'primary_amine_enumeration_for_chodera_lab_FEP-permuted-conformers-x10789.sdf',
    #'boronic_ester_enumeration' : 'boronic_ester_enumeration_for_chodera_lab_FEP-permuted-conformers-x10789.sdf',
    #'retrospective-aminopyridines' : 'activity-data-2020-07-29-conformers-x10789.sdf',
    #'retrospective-aminopyridines-matt' : 'aminopyridine_compounds_for_FEP_benchmarking-conformers-x10789.sdf',
    #'retrospective-aminopyridines-matt-dockscores' : 'aminopyridine_compounds_for_FEP_benchmarking-dockscores-x10789.sdf',
    #'fastgrant-table1' : 'fastgrant-table1-dockscores-x10789.sdf',
    # 'RAL-THA-6b94ceba' : 'RAL-THA-6b94ceba-dockscores-x10789.sdf',

    # Sprint 3
    #'nucleophilic-displacement-x10876' : '2020-08-20-benzotriazoles-dockscores-x10876.sdf',
    #'nucleophilic-displacement-x10871' : '2020-08-20-benzotriazoles-dockscores-x10871.sdf',
    ##'nucleophilic-displacement-x10820' : '2020-08-20-benzotriazoles-dockscores-x10820.sdf',

    # retrospective
    'benzotriazoles-retrospective-x10876' : '2020-09-01-benzotriazoles-retrospective-dockscores-x10876.sdf',
    'benzotriazoles-retrospective-x10871' : '2020-09-01-benzotriazoles-retrospective-dockscores-x10871.sdf',
    'benzotriazoles-retrospective-x10820' : '2020-09-01-benzotriazoles-retrospective-dockscores-x10820.sdf',
}

# Sprint 2
#json_filename = '2020-08-20-benzotriazoles.json'
json_filename = '2020-09-01-benzotriazoles-retrospective.json'

# SMARTS for common core scaffold
smarts = 'c1ccc(NC(=O)[C,N]n2nnc3ccccc32)cc1'

receptors = [
    '../receptors/monomer/Mpro-x2646_0_bound-protein.pdb',
    '../receptors/monomer/Mpro-x2646_0_bound-protein-thiolate.pdb',
    '../receptors/monomer/Mpro-x10876_0_bound-protein.pdb',
    '../receptors/monomer/Mpro-x10876_0_bound-protein-thiolate.pdb',
    '../receptors/monomer/Mpro-x10871_0_bound-protein.pdb',
    '../receptors/monomer/Mpro-x10871_0_bound-protein-thiolate.pdb',
    '../receptors/monomer/Mpro-x10820_0_bound-protein.pdb',
    '../receptors/monomer/Mpro-x10820_0_bound-protein-thiolate.pdb',
]

def get_moldb(series):
    lig_file_name = ligand_dict[series]
    moldb = oechem.OEMolDatabase()
    moldb.Open(lig_file_name)    
    return moldb

def add_run(master_dict, i, mol_i, j, mol_j, starting_index=0, **kwargs):
    index = len(master_dict) + starting_index

    smiles_flag = oechem.OESMILESFlag_Canonical | oechem.OESMILESFlag_ISOMERIC

    master_dict[index] = {
        'JOBID':index,
        'directory':f'RUN{index}',
        
        'start':i,
        'start_title':mol_i.GetTitle(),
        'start_smiles':oechem.OECreateSmiString(mol_i, smiles_flag),

        'end':j,
        'end_title':mol_j.GetTitle(),
        'end_smiles':oechem.OECreateSmiString(mol_j,  smiles_flag),                        

        **kwargs
    }
    if oechem.OEHasSDData(mol_j,'f_avg_pIC50'):
        master_dict[index]['start_pIC50'] = oechem.OEGetSDData(mol_j, 'f_avg_pIC50')
    if oechem.OEHasSDData(mol_i,'f_avg_pIC50'):
        master_dict[index]['end_pIC50'] = oechem.OEGetSDData(mol_i, 'f_avg_pIC50')

def add_direction(master_dict, direction, ligand_dict, reference_molecule_index=0, starting_index=0):
    i = reference_molecule_index
    mol_i = oechem.OEGraphMol()
    mol_j = oechem.OEGraphMol()

    for series in ligand_dict:
        print(series)

        moldb = get_moldb(series)
        num_mols = moldb.NumMols()
        moldb.GetMolecule(mol_i, reference_molecule_index)

        for protein in receptors:
            print(protein)

            kwargs = { 
                'series':series,
                'direction':direction,
                'smarts':smarts,
                'target':'SARS-CoV-2 Mpro',
                'protein':protein,
                'ligand':ligand_dict[series],
                'ff':ff,
            }

            for j in range(1, num_mols):
                # Get the molecule
                moldb.GetMolecule(mol_j, j)

                if direction == 'forwards':
                    add_run(master_dict, i, mol_i, j, mol_j, starting_index=starting_index, **kwargs)

                if direction == 'backwards':                    
                    add_run(master_dict, j, mol_j, i, mol_i, starting_index=starting_index, **kwargs)

# Sprint 3
starting_index = 2688 # retrospective start
master_dict = dict()
add_direction(master_dict, 'backwards', ligand_dict, starting_index=starting_index)
add_direction(master_dict, 'forwards', ligand_dict, starting_index=starting_index)

with open(json_filename, "w") as f:
    json.dump(master_dict, f, sort_keys=True, indent=4, separators=(',', ': '))

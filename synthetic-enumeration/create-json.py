import numpy as np
import json
import math
import itertools

ff = 'openff-1.2.0'
# protein_dict = {'Series5':'protein.pdb'}
ligand_dict = {
    'primary_amine_enumeration' : 'primary_amine_enumeration_for_chodera_lab_FEP-permuted-conformers.sdf',
#    'boronic_ester_enumeration' : 'boronic_ester_enumeration_for_chodera_lab_FEP-permuted-conformers.sdf',
}
receptors = ['receptors/monomer/Mpro-x2646_0_bound-protein.pdb']
index = 0 # starting index

# Count ligands
print('Counting ligands...')
def count_molecules(lig_file_name):
    moldb = oechem.OEMolDatabase()
    moldb.Open(lig_file_name)
    num_mols = moldb.NumMols()
    return num_mols

n_ligands = dict()
from openeye import oechem 
for series in ligand_dict:
    n_ligands[series] = count_molecules(ligand_dict[series])
print(n_ligands)

i = 0 # reference molecule
master_dict = dict()
for series in ligand_dict:
    for protein in receptors:
        for j in range(1, n_ligands[series]):
            # Only include forward directions
            master_dict[index] = {'target':'backtesting','start':i,'end':j,
                                  'protein':protein,'ligand':ligand_dict[series],'ff':ff,
                                           'directory':f'RUN{index-1}','JOBID':index}
            index += 1

with open(f"2020-07-24.json", "w") as f:
    json.dump(master_dict, f, sort_keys=True, indent=4, separators=(',', ': '))

import numpy as np
import json
import math
import itertools
from openeye import oechem

ff = 'openff-1.2.0'
# protein_dict = {'Series5':'protein.pdb'}
ligand_dict = {
    'primary_amine_enumeration' : 'primary_amine_enumeration_for_chodera_lab_FEP-permuted-conformers-x10789.sdf',
    'boronic_ester_enumeration' : 'boronic_ester_enumeration_for_chodera_lab_FEP-permuted-conformers-x10789.sdf',
}
receptors = [
    '../receptors/monomer/Mpro-x2646_0_bound-protein.pdb',
    '../receptors/monomer/Mpro-x2646_0_bound-protein-thiolate.pdb'
]

#index = 1686 # starting index
index = 4664 # starting index

mol_i = oechem.OEGraphMol()
mol_j = oechem.OEGraphMol()

i = 0 # reference molecule
master_dict = dict()
smiles_flag = oechem.OESMILESFlag_Canonical | oechem.OESMILESFlag_ISOMERIC
for series in ligand_dict:
    print(series)

    # Open ligands as OEMolDatabase
    lig_file_name = ligand_dict[series]
    moldb = oechem.OEMolDatabase()
    moldb.Open(lig_file_name)
    num_mols = moldb.NumMols()
    
    moldb.GetMolecule(mol_i, i)

    for protein in receptors:
        print(protein)

        for j in range(1, num_mols):
            # Get the molecule
            moldb.GetMolecule(mol_j, j)

            # Only include backward directions
            master_dict[index] = {
                'JOBID':index,
                'directory':f'RUN{index-1}',

                'target':'SARS-CoV-2 Mpro',

                'start':j,
                'start_title':mol_j.GetTitle(),
                'start_smiles':oechem.OECreateSmiString(mol_j, smiles_flag),

                'end':i,
                'end_title':mol_i.GetTitle(),
                'end_smiles':oechem.OECreateSmiString(mol_i, smiles_flag),

                'protein':protein,
                'ligand':ligand_dict[series],
                'ff':ff,
            }
            index += 1

with open(f"2020-07-28.json", "w") as f:
    json.dump(master_dict, f, sort_keys=True, indent=4, separators=(',', ': '))

import numpy as np
import json
import math
import itertools
from openeye import oechem

ff = 'openff-1.2.0'
ligand_dict = {
    #'primary_amine_enumeration' : 'primary_amine_enumeration_for_chodera_lab_FEP-permuted-conformers-x10789.sdf',
    #'boronic_ester_enumeration' : 'boronic_ester_enumeration_for_chodera_lab_FEP-permuted-conformers-x10789.sdf',
    #'retrospective-aminopyridines' : 'activity-data-2020-07-29-conformers-x10789.sdf',
    #'retrospective-aminopyridines-matt' : 'aminopyridine_compounds_for_FEP_benchmarking-conformers-x10789.sdf',
    #'retrospective-aminopyridines-matt-dockscores' : 'aminopyridine_compounds_for_FEP_benchmarking-dockscores-x10789.sdf',
    #'fastgrant-table1' : 'fastgrant-table1-dockscores-x10789.sdf',
    'RAL-THA-6b94ceba' : 'RAL-THA-6b94ceba-dockscores-x10789.sdf',
}
receptors = [
    '../receptors/monomer/Mpro-x2646_0_bound-protein.pdb',
    '../receptors/monomer/Mpro-x2646_0_bound-protein-thiolate.pdb'
]

#index = 1686 # starting index
#index = 4664 # starting index
#index = 7642 # starting index for retrospective-aminopyridines
#index = 7802 # starting index for retrospective-aminopyridines-matt
#index = 8086 # starting index for retrospective-aminopyridines-matt
#index = 8374 # starting index for fastgrant-table1
index = 8410 # starting index for RAL-THA-6b94ceba

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

            # Add forwards
            master_dict[index] = {
                'JOBID':index,
                'directory':f'RUN{index-1}',
                'series':series,

                'target':'SARS-CoV-2 Mpro',

                'start':i,
                'start_title':mol_i.GetTitle(),
                'start_smiles':oechem.OECreateSmiString(mol_i, smiles_flag),

                'end':j,
                'end_title':mol_j.GetTitle(),
                'end_smiles':oechem.OECreateSmiString(mol_j, smiles_flag),

                'protein':protein,
                'ligand':ligand_dict[series],
                'ff':ff,
            }
            if oechem.OEHasSDData(mol_i,'f_avg_pIC50'):
                master_dict[index]['start_pIC50'] = oechem.OEGetSDData(mol_i, 'f_avg_pIC50')
            if oechem.OEHasSDData(mol_j,'f_avg_pIC50'):
                master_dict[index]['end_pIC50'] = oechem.OEGetSDData(mol_j, 'f_avg_pIC50')
            index += 1

            # Add backwards
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
            if oechem.OEHasSDData(mol_j,'f_avg_pIC50'):
                master_dict[index]['start_pIC50'] = oechem.OEGetSDData(mol_j, 'f_avg_pIC50')
            if oechem.OEHasSDData(mol_i,'f_avg_pIC50'):
                master_dict[index]['end_pIC50'] = oechem.OEGetSDData(mol_i, 'f_avg_pIC50')
            index += 1

with open(f"2020-08-12-RAL-THA-6b94ceba1.json", "w") as f:
    json.dump(master_dict, f, sort_keys=True, indent=4, separators=(',', ': '))

#!/bin/env python

"""
Extract retrospective subset of docked SDF files that have experimental data and are enantiomerically resolved
"""

import numpy as np

input_sdf_filename = 'docked/sprint-10-2021-07-26-x10959-dimer-neutral.sdf'
output_prefix = 'sprint-10-2021-07-26-x10959-dimer-neutral'
output_dir = 'retrospective-for-orion/'

# Read all molecules from the SDF file
molecules = list()
from openeye import oechem
with oechem.oemolistream(input_sdf_filename) as ifs:
    for mol in ifs.GetOEGraphMols():
        molecules.append( oechem.OEGraphMol(mol) )
print(f'Loaded data for {len(molecules)} molecules.')

def has_scaffold(mol):
    """Return True if this molecule has desired scaffold"""
    from openeye import oechem
    mol = oechem.OEGraphMol(mol)
    scaffold_smiles = 'O=C(Cc1cccc(Cl)c1)N'
    ss = oechem.OESubSearch(scaffold_smiles)
    oechem.OEPrepareSearch(mol, ss)
    if ss.SingleMatch(mol):
        return True
    else:
        return False

def has_ic50(mol):
    """Return True if this molecule has fluorescence IC50 data"""
    from openeye import oechem
    if not oechem.OEHasSDData(mol, 'f_avg_IC50'):
        return False

    try:
        if oechem.OEHasSDData(mol, 'f_avg_IC50'):
            IC50 = oechem.OEGetSDData(mol, 'f_avg_IC50')
            IC50 = float(IC50)
            if IC50 < 99.0:
                return True
            else:
                return False
        else:
            return False
    except Exception as e:
        return False

# Only retain those with experimental data
molecules = [ molecule for molecule in molecules if has_ic50(molecule) ]
print(f'{len(molecules)} remain after retaining only those with experimental data')

# Only retain those with only one stereoisomer
import re
molecule_names = set([ molecule.GetTitle() for molecule in molecules ])
molecules = [ molecule for molecule in molecules if (re.sub(r'_1$', r'_2', molecule.GetTitle()) not in molecule_names) ]
print(f'{len(molecules)} remain after retaining only those with single enantiomers')

# Only retain those with scaffold
molecules = [ molecule for molecule in molecules if has_scaffold(molecule) ]
print(f'{len(molecules)} remain after retaining only those with scaffold')

# Write molecules
import os
os.makedirs(output_dir, exist_ok=True)
for extension in ['smi', 'sdf', 'mol2', 'pdb']:
    output_filename = f'{output_dir}/{output_prefix}.{extension}'
    with oechem.oemolostream(output_filename) as ofs:
        for index, mol in enumerate(molecules):
            oechem.OEWriteMolecule(ofs, oechem.OEGraphMol(mol))
        print(f'{len(molecules)} written to {output_filename}')

# Write experimental data
output_filename = f'{output_dir}/experimental-data'
with open(output_filename, 'wt') as outfile:
    for index, mol in enumerate(molecules):
        IC50 = oechem.OEGetSDData(mol, 'f_avg_IC50')
        IC50 = float(IC50) * 1.0e-6
        DeltaG_in_kcal_per_mol = np.log(IC50) * 0.596
        dDeltaG_in_kcal_per_mol = 0.3 # rough guess; fix later
        outfile.write(f'{mol.GetTitle()} {DeltaG_in_kcal_per_mol:.3f} {dDeltaG_in_kcal_per_mol:.3f} kcal/mol\n')
print(f'{len(molecules)} written to {output_filename}')

# Write transformations
output_filename = f'{output_dir}/transformations'
reference_molecule_title = 'ADA-UCB-6c2cb422-1_1'
with open(output_filename, 'wt') as outfile:
    outfile.write(f'; Target Mpro - Sprint 10 retrospective - num edges: {len(molecules)-1}\n')
    for index, mol in enumerate(molecules):
        outfile.write(f'{reference_molecule_title} >> {mol.GetTitle()}\n')
print(f'{len(molecules)-1} transformations written to {output_filename}')

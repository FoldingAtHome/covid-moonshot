#!/bin/env python

"""
Extract retrospective subset of docked SDF files that have experimental data and are enantiomerically resolved
"""

input_sdf_filename = 'docked/sprint-9-2021-06-14-x10959-dimer-neutral.sdf'
output_sdf_filename = 'retrospective-for-siremol/sprint-9-2021-06-14-x10959-dimer-neutral-retrospective-docked.sdf'

# Read all molecules from the SDF file
molecules = list()
from openeye import oechem 
with oechem.oemolistream(input_sdf_filename) as ifs:
    for mol in ifs.GetOEGraphMols():
        molecules.append( oechem.OEGraphMol(mol) )
print(f'Loaded data for {len(molecules)} molecules.')

def has_ic50(mol):
    """Return True if this molecule has fluorescence IC50 data"""
    from openeye import oechem
    if not oechem.OEHasSDData(mol, 'f_avg_pIC50'):
        return False

    try:
        if oechem.OEHasSDData(mol, 'f_avg_IC50'):
            IC50 = oechem.OEGetSDData(mol, 'f_avg_IC50')
            IC50 = float(IC50)
            return True
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

# Write molecules
with oechem.oemolostream(output_sdf_filename) as ofs:
    for mol in molecules:
        oechem.OEWriteMolecule(ofs, mol)
print(f'{len(molecules)} written to {output_sdf_filename}')

# Write molecules as separate mol2 files
for index, mol in enumerate(molecules):
    with oechem.oemolostream(f'retrospective-for-siremol/molecules/{mol.GetTitle()}.sdf') as ofs:
        oechem.OEWriteMolecule(ofs, oechem.OEGraphMol(mol))
    with oechem.oemolostream(f'retrospective-for-siremol/molecules/{mol.GetTitle()}.mol2') as ofs:
        oechem.OEWriteMolecule(ofs, oechem.OEGraphMol(mol))
    with oechem.oemolostream(f'retrospective-for-siremol/molecules/{mol.GetTitle()}.pdb') as ofs:
        oechem.OEWriteMolecule(ofs, oechem.OEGraphMol(mol))


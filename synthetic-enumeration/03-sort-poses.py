from openeye import oechem

"""
Sort poses by docking score
"""

prefix = 'nucleophilic_displacement_enumeration_for_FEP'
fragment = 'x10789'
unsorted_filename = f'{prefix}-permuted-dockscores-{fragment}.sdf'
sorted_filename =  f'{prefix}-sorted-{fragment}.sdf'

mols = list()
print('Reading molecules...')
with oechem.oemolistream(unsorted_filename) as ifs:
    for mol in ifs.GetOEGraphMols():
        mols.append(oechem.OEGraphMol(mol))

print('Pruning molecules...')
mols = [ mol for mol in mols if mol.NumAtoms() > 1 ]

def score(mol):
    text = oechem.OEGetSDData(mol, 'docking_score')
    try:
        score = float(text)
        return score
    except Exception as e:
        return 0

print('Sorting molecules...')
mols = sorted(mols, key=score)

print('Writing molecules...')
with oechem.oemolostream(sorted_filename) as ofs:
    for mol in mols:
        oechem.OEWriteMolecule(ofs, mol)

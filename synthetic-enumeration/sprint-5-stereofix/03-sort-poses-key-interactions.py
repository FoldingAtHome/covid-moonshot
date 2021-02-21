from openeye import oechem

"""
Sort poses by Tim Dudgen NumKeyInteractions score, and then interleave by intermediate
"""

# TODO: Automatically list these with -microstates-*.sdf pattern
fragments = [
    'x11498',
    'x12073'
]

prefix = 'sprint-5'

def compare(mol1, mol2):
    clash_score_1 = float(oechem.OEGetSDData(mol1, 'clash_score'))
    clash_score_2 = float(oechem.OEGetSDData(mol2, 'clash_score'))

    num_interactions_1 = int(oechem.OEGetSDData(mol1, 'NumTotalInteractions'))
    num_interactions_2 = int(oechem.OEGetSDData(mol2, 'NumTotalInteractions'))

    dock_score_1 = float(oechem.OEGetSDData(mol1, 'docking_score'))
    dock_score_2 = float(oechem.OEGetSDData(mol2, 'docking_score'))

    if (clash_score_1 == 0) and (clash_score_2 == 0):
        if num_interactions_1 != num_interactions_2:
            return num_interactions_2 - num_interactions_1
        else:
            return dock_score_1 - dock_score_2
    elif (clash_score_1 == 0) or (clash_score_2 == 0):
        return clash_score_1 - clash_score_2
    else:
        return dock_score_1 - dock_score_2

for fragment in fragments:
    # Input filename
    unsorted_filename = f'docked/{prefix}-microstates-{fragment}.sdf'
    key_interactions_filename = f'docked/{prefix}-microstates-{fragment}-inters.sdf'

    # Output filename
    sorted_filename = f'docked/{prefix}-microstates-{fragment}-sorted.sdf'

    # Read number of key interactions
    num_key_interactions = dict()
    print('Reading number of key interactions')
    with oechem.oemolistream(key_interactions_filename) as ifs:
        for mol in ifs.GetOEGraphMols():
            num_key_interactions[mol.GetTitle()] = oechem.OEGetSDData(mol, 'NumTotalInteractions')

    # Read molecules
    mols = list()
    print('Reading molecules...')
    with oechem.oemolistream(unsorted_filename) as ifs:
        for mol in ifs.GetOEGraphMols():
            oechem.OESetSDData(mol, 'NumTotalInteractions', num_key_interactions[mol.GetTitle()])
            mols.append(oechem.OEGraphMol(mol))
    print(f'{len(mols)} molecules')

    print('Pruning empty molecules...')
    mols = [ mol for mol in mols if mol.NumAtoms() > 1 ]
    print(f'{len(mols)} molecules')

    print('Sorting molecules by number of key interactions score...')
    from functools import cmp_to_key
    mols = sorted(mols[1:], key=cmp_to_key(compare))
    print(f'{len(mols)} molecules')

    # Interleave molecules by intermediate category
    print('Interleaving molecules by intermediate...')
    intermediates = set([oechem.OEGetSDData(mol, 'intermediate') for mol in mols])
    mols_by_label = { label : list() for label in intermediates }
    for mol in mols:
        label = oechem.OEGetSDData(mol, 'intermediate')
        mols_by_label[label].append(mol)
    nmols_to_assign = len(mols)
    mols = list()
    while (nmols_to_assign > 0):
        for label in mols_by_label.keys():
            if len(mols_by_label[label]) > 0:
                mols.append(mols_by_label[label].pop(0))
                nmols_to_assign -= 1
    print(f'{len(mols)} molecules')

    # Move important molecules to the beginning
    print('Moving important molecules to the beginning...')
    import re
    pattern = '\w{3}-\w{3}-[\w\d]{8}-\d+'
    important_molecules = [ mol for mol in mols if re.search(pattern, mol.GetTitle()) is not None ]
    unimportant_molecules = [ mol for mol in mols if re.search(pattern, mol.GetTitle()) is None ]
    mols = important_molecules + unimportant_molecules
    print(f'{len(mols)} molecules')

    print('Writing molecules...')
    with oechem.oemolostream(sorted_filename) as ofs:
        for mol in mols:
            oechem.OEWriteMolecule(ofs, mol)

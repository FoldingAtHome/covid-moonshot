from openeye import oechem

"""
Sort poses by docking score, and then interleave by intermediate
"""

# TODO: Automatically list these with -microstates-*.sdf pattern
fragments = [
    'x11498',
    'x12073'
]

prefix = 'sprint-5'

def score(mol):
    text = oechem.OEGetSDData(mol, 'docking_score')
    try:
        score = float(text)
        return score
    except Exception as e:
        return 0

for fragment in fragments:
    # Input filename
    unsorted_filename = f'docked/{prefix}-microstates-{fragment}.sdf'

    # Output filename
    sorted_filename = f'docked/{prefix}-microstates-{fragment}-sorted.sdf'

    # Read molecules
    mols = list()
    print('Reading molecules...')
    with oechem.oemolistream(unsorted_filename) as ifs:
        for mol in ifs.GetOEGraphMols():
            mols.append(oechem.OEGraphMol(mol))
    print(f'{len(mols)} molecules')

    print('Pruning empty molecules...')
    mols = [ mol for mol in mols if mol.NumAtoms() > 1 ]
    print(f'{len(mols)} molecules')

    print('Sorting molecules by docking score...')
    mols = sorted(mols[1:], key=score)
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

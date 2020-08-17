from openeye import oechem

winners = [
    'EN300-60314',
    'EN300-6734624',
    'EN300-20814457',
    'EN300-1723947',
    'EN300-2515954',
    'EN300-295395',
    'EN300-6734624',
    'EN300-124496',
    'EN300-212829',
    'EN300-26624333',
    'EN300-365771',
    'EN300-316592',
    'EN300-298506',
    'EN300-1086686',
    'EN300-106778',
    'EN300-1425849',
    'EN300-1704613',
    'EN300-299518',
    ]

mols = list()
for filename in [
    'primary_amine_enumeration_for_chodera_lab_FEP-permuted-conformers-x10789.sdf',
    'boronic_ester_enumeration_for_chodera_lab_FEP-permuted-conformers-x10789.sdf'
]:
    with oechem.oemolistream(filename) as ifs:
        for mol in ifs.GetOEGraphMols():
            if mol.GetTitle() in winners:
                mols.append(oechem.OEGraphMol(mol))

print(f'{len(mols)} molecules read')

print('Writing molecules...')
with oechem.oemolostream('sprint1-winners.mol2') as ofs:
    for mol in mols:
        oechem.OEWriteMolecule(ofs, mol)

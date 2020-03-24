#fragid = '0104'
fragid = '0387'

complex_pdb_filename = f'../diamond-structures/Mpro full XChem screen - pdbs - active site non-covalent - ver-2020-03-20/Mpro-x{fragid}.pdb'
#complex_pdb_filename = '../diamond-structures/5r7z.pdb'

def read_pdb_file(pdb_file):
    print(f'Reading receptor from {pdb_file}...')

    from openeye import oechem
    ifs = oechem.oemolistream()
    ifs.SetFlavor(oechem.OEFormat_PDB, oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)  # noqa

    if not ifs.open(pdb_file):
        oechem.OEThrow.Fatal("Unable to open %s for reading." % pdb_file)

    temp_mol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ifs, temp_mol):
        oechem.OEThrow.Fatal("Unable to read molecule from %s." % pdb_file)
    ifs.close()

    fact = oechem.OEAltLocationFactory(temp_mol)
    mol = oechem.OEGraphMol()
    mol.Clear()
    fact.MakePrimaryAltMol(mol)
    return (mol)

# Read the receptor and identify design units
from openeye import oespruce, oechem
complex = read_pdb_file(complex_pdb_filename)
print('Identifying design units...')
design_units = list(oespruce.OEMakeDesignUnits(complex))
for i, design_unit in enumerate(design_units):
    filename = f'DU_{i}.oedu'
    print(f'Writing design unit {i} to {filename}')
    oechem.OEWriteDesignUnit(filename, design_unit)
if len(design_units) > 1:
    print('More than one design unit found---using first one')
    design_unit = design_units[0]

# Prepare the receptor
print('Preparing receptor...')
from openeye import oedocking
protein = oechem.OEGraphMol()
design_unit.GetProtein(protein)
ligand = oechem.OEGraphMol()
design_unit.GetLigand(ligand)
receptor = oechem.OEGraphMol()
oedocking.OEMakeReceptor(receptor, protein, ligand)
oedocking.OEWriteReceptorFile(receptor, f'receptor-{fragid}.oeb.gz')

with oechem.oemolostream(f'protein-{fragid}.pdb') as ofs:
    oechem.OEWriteMolecule(ofs, protein)
with oechem.oemolostream(f'ligand-{fragid}.mol2') as ofs:
    oechem.OEWriteMolecule(ofs, ligand)

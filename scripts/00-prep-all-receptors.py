"""
Prepare all receptors for OEDocking

"""

import glob, os
source_pdb_files = glob.glob("../diamond-structures/Mpro_All_PDBs - ver 2020-03-24/*.pdb")

def read_pdb_file(pdb_file):
    print(f'Reading receptor from {pdb_file}...')

    from openeye import oechem
    ifs = oechem.oemolistream()
    ifs.SetFlavor(oechem.OEFormat_PDB, oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)  # noqa

    if not ifs.open(pdb_file):
        oechem.OEThrow.Fatal("Unable to open %s for reading." % pdb_file)

    mol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ifs, mol):
        oechem.OEThrow.Fatal("Unable to read molecule from %s." % pdb_file)
    ifs.close()

    return (mol)

def prepare_receptor(complex_pdb_filename, prefix):
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
    oedocking.OEWriteReceptorFile(receptor, f'{prefix}-receptor.oeb.gz')

    with oechem.oemolostream(f'{prefix}-protein.pdb') as ofs:
        oechem.OEWriteMolecule(ofs, protein)
    with oechem.oemolostream(f'{prefix}-ligand.mol2') as ofs:
        oechem.OEWriteMolecule(ofs, ligand)
    with oechem.oemolostream(f'{prefix}-ligand.pdb') as ofs:
        oechem.OEWriteMolecule(ofs, ligand)
    with oechem.oemolostream(f'{prefix}-ligand.sdf') as ofs:
        oechem.OEWriteMolecule(ofs, ligand)

# Process all receptors
for complex_pdb_filename in source_pdb_files:
    basepath, filename = os.path.split(complex_pdb_filename)
    prefix, extension = os.path.splitext(filename)
    import os
    if not os.path.exists(f'{prefix}-receptor.oeb.gz'):
        try:
            print(prefix)
            prepare_receptor(complex_pdb_filename, prefix)
        except Exception as e:
            print(e)

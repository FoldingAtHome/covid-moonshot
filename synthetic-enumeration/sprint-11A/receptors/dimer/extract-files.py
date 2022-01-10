from openeye import oechem
import copy

prefix = 'Mpro-P1800_0A_bound'

dus = list()
du = oechem.OEDesignUnit()
filename = f'{prefix}.oedu'
ifs = oechem.oeifstream(filename)
while oechem.OEReadDesignUnit(ifs, du):
    break

oemol = oechem.OEGraphMol()

du.GetComponents(oemol, oechem.OEDesignUnitComponents_Protein | oechem.OEDesignUnitComponents_Ligand | oechem.OEDesignUnitComponents_Solvent)
with oechem.oemolostream(f'{prefix}-complex.pdb') as ofs:
    oechem.OEWriteMolecule(ofs, oemol)

du.GetComponents(oemol, oechem.OEDesignUnitComponents_Protein | oechem.OEDesignUnitComponents_Solvent)
with oechem.oemolostream(f'{prefix}-protein.pdb') as ofs:
    oechem.OEWriteMolecule(ofs, oemol)

du.GetComponents(oemol, oechem.OEDesignUnitComponents_Ligand)
with oechem.oemolostream(f'{prefix}-ligand.pdb') as ofs:
    oechem.OEWriteMolecule(ofs, oemol)

du.GetComponents(oemol, oechem.OEDesignUnitComponents_Ligand)
with oechem.oemolostream(f'{prefix}-ligand.sdf') as ofs:
    oechem.OEWriteMolecule(ofs, oemol)
with oechem.oemolostream(f'{prefix}-ligand.mol2') as ofs:
    oechem.OEWriteMolecule(ofs, oemol)

from openeye import oedocking
recOpts = oedocking.OEMakeReceptorOptions()
result = oedocking.OEMakeReceptor(du, recOpts)
#ofs = oechem.oeofstream(f'{prefix}-receptor.oeb.gz')
#oechem.OEWriteDesignUnit(ofs, du)
protein = oechem.OEGraphMol()
du.GetComponents(protein, oechem.OEDesignUnitComponents_Protein | oechem.OEDesignUnitComponents_Solvent) 
ligand = oechem.OEGraphMol()
du.GetComponents(ligand, oechem.OEDesignUnitComponents_Ligand) 
receptor = oechem.OEGraphMol()
oedocking.OEMakeReceptor(receptor, protein, ligand, False)
oedocking.OEWriteReceptorFile(oechem.OEGraphMol(receptor), f'{prefix}-receptor.oeb.gz')


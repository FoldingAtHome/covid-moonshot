"""
Check all receptors can be set up in OpenMM
"""

from simtk import openmm, unit
from simtk.openmm import app

protein_forcefield = 'amber14/protein.ff14SB.xml'
forcefield = app.ForceField(protein_forcefield)

import glob
files = glob.glob('../receptors/Mpro-*-protein.pdb')
for filename in files:
    print(filename)

    # Filter out UNK
    pdbfile_contents = [ line for line in open(filename, 'r') if 'UNK' not in line ]
    from io import StringIO
    infile = StringIO(''.join(pdbfile_contents))

    pdbfile = app.PDBFile(infile)
    try:
        system = forcefield.createSystem(pdbfile.topology)
    except Exception as e:
        print(e)

# Read compound IDs
cids = set()
import csv
with open('merged-compounds-2020-04-01.csv', 'r') as csvfile:
    spamreader = csv.reader(csvfile)
    for row in spamreader:
        cid = row[1]
        cids.add(cid)

from openeye import oechem
mol = oechem.OEGraphMol()
with oechem.oemolistream('../moonshot-submissions/covid_submissions_03_31_2020-docked.csv') as ifs:
    with oechem.oemolostream('merged-compounds-2020-04-01-docked.csv') as ofs:
        while oechem.OEReadMolecule(ifs, mol):
            if mol.GetTitle() in cids:
                print(mol.GetTitle())
                oechem.OEWriteMolecule(ofs, mol)


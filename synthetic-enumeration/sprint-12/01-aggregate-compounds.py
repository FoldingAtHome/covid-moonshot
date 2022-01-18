from openeye import oechem
from rich.progress import track

"""
Prepare list of annotated compounds for docking

* Aggregate all compound designs
* Annotate by intermediate
* Sort compounds by size
* Interleave compounds by intermediate

"""

mols = list()

# TODO: Auto-download and timestamp
submissions_csv_filename = 'submissions/submissions-2022-01-17.csv.bz2'

# Read all submitted designs: Compounds with the key substructure will be retained
print('Reading submitted designs...')
# Drop columns that cause trouble for OpenEye
import pandas as pd
drop_columns = ['Submission Rationale', 'Submission Notes']
df = pd.read_csv(submissions_csv_filename, dtype=str)
df.drop(columns=drop_columns, inplace=True)
print(f'{len(df)} records read')
import tempfile
with tempfile.NamedTemporaryFile(suffix='.csv') as csv_file:
    df.to_csv(csv_file.name, header=True, index=False)
    # Read file
    with oechem.oemolistream(csv_file.name) as ifs:
        mol = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, mol):
            # Clear SD tags
            oechem.OEClearSDData(mol)
            # Store the molecule
            mols.append(mol.CreateCopy())
print(f'{len(mols)} molecules read')

# Aggregate all compound designs
print('reading compound designs...')
source_filenames = [
    # Filtered synthetic designs
    #'virtual-synthetic-library/reductive-aminations.csv',
]
for source_filename in source_filenames:
    print(source_filename)
    with oechem.oemolistream(source_filename) as ifs:
        mol = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, mol):
            # Clear SD tags
            oechem.OEClearSDData(mol)
            # Store the molecule
            mols.append(mol.CreateCopy())
print(f'{len(mols)} molecules read')

# Annotate molecules with SMARTS labels
print('Annotating SMARTS labels...')
import csv
labels_filename = 'annotations/annotations.csv' # list of labels for various SMARTS patterns
smarts_labels = dict()
with open(labels_filename, 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    for row in csvreader:
        smarts = row[0]
        label = row[1]
        if smarts[0] != '#': # allow comments
            smarts_labels[smarts] = label
# Label the molecules
for smarts, label in smarts_labels.items():
    ss = oechem.OESubSearch(smarts)
    for mol in track(mols, description=label):
        oechem.OEPrepareSearch(mol, ss)
        if ss.SingleMatch(mol):
            oechem.OESetSDData(mol, 'scaffold', label)
# Discard molecules without labels
mols = [ mol for mol in mols if oechem.OEHasSDData(mol, 'scaffold') ]
print(f'{len(mols)} molecules remain after discarding molecules that do not match scaffold')

# Sort based on molecular weight
print(f'Sorting molecules by number of atoms...')
mols.sort(key=lambda mol : mol.NumAtoms())

# Filter the number of heavy atoms
#n_heavy_max = 35
n_heavy_max = 50
mols = [ mol for mol in mols if oechem.OECount(mol, oechem.OEIsHeavy()) <= n_heavy_max ]
print(f'{len(mols)} molecule remain after filtering atoms with more than {n_heavy_max} heavy atoms')

# Write molecules
output_filename = 'sorted/sprint-12.csv'
with oechem.oemolostream(output_filename) as ofs:
    for mol in track(mols, description='Writing molecules...'):
        oechem.OEWriteMolecule(ofs, mol)

# Generate PDF
pdf_filename = output_filename.replace('.csv', '.pdf')
print(f'Generating PDF as {pdf_filename}')
cmd = f"mols2pdf.py -in {output_filename} -out {pdf_filename}"
import os
os.system(cmd)

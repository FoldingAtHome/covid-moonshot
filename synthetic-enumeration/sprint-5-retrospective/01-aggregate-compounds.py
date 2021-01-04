from openeye import oechem
from rich.progress import track

"""
Prepare list of annotated compounds for docking

* Aggregate all compound designs
* Annotate by intermediate
* Sort compounds by size
* Interleave compounds by intermediate

"""

# Latest submissions downloaded from PostEra site
# TODO: Auto-download and timestamp
submissions_csv_filename = 'activity-data/activity-data-2020-12-30.csv'
#submissions_csv_filename = 'submissions/submissions-2020-12-30.csv'

# Aggregate all compound designs
source_filenames = [
    # Filtered synthetic designs
    #'filtered/filtered_alkyl_halide_bb_ether_synthesis_for_FEP.smi',
    #'filtered/filtered_amine_bb_amide_couplings_for_FEP.smi',
    #'filtered/filtered_cooh_bb_amide_couplings_for_FEP.smi',
]
mols = list()
print('Reading compound designs...')
for source_filename in source_filenames:
    with open(source_filename, 'rt') as infile:
        for line in track(infile.readlines(), description=f'Reading {source_filename}'):
            smiles, title, _, demerits, reason_for_demerits = line.split()
            mol = oechem.OEGraphMol()
            oechem.OESmilesToMol(mol, smiles)
            mol.SetTitle(title)
            oechem.OESetSDData(mol, 'demerits', demerits)
            oechem.OESetSDData(mol, 'reason_for_demerits', reason_for_demerits)
            mols.append(mol)
print(f'{len(mols)} molecules read')

# Read all submitted designs: Compounds with the key substructure will be retained
print('Reading submitted designs...')
# Drop columns that cause trouble for OpenEye
import pandas as pd
df = pd.read_csv(submissions_csv_filename)
#drop_columns = ['Submission Rationale', 'Submission Notes']
#df.drop(columns=drop_columns, inplace=True)
import tempfile
with tempfile.NamedTemporaryFile(suffix='.csv') as csv_file:
    df.to_csv(csv_file.name, header=True, index=False)
    # Read file
    tags_to_retain = ['f_avg_IC50']
    with oechem.oemolistream(csv_file.name) as ifs:
        mol = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, mol):
            sdtags = { tag : oechem.OEGetSDData(mol, tag) for tag in tags_to_retain }
            # Clear SD tags
            oechem.OEClearSDData(mol)
            # Repopulate
            [ oechem.OESetSDData(mol, tag, value) for tag, value in sdtags.items() ]
            # Store the molecule
            mols.append(mol.CreateCopy())
print(f'{len(mols)} molecules read')

# Annotate molecules with SMARTS labels
print('Annotating SMARTS labels...')
import csv
labels_filename = 'annotations/benzopyran_annotations.csv' # list of labels for various SMARTS patterns
smarts_labels = dict()
with open(labels_filename, 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    for row in csvreader:
        smarts = row[0]
        label = row[1]
        smarts_labels[smarts] = label
# Label the molecules
for smarts, label in smarts_labels.items():
    ss = oechem.OESubSearch(smarts)
    for mol in track(mols, description=label):
        oechem.OEPrepareSearch(mol, ss)
        if ss.SingleMatch(mol):
            oechem.OESetSDData(mol, 'intermediate', label)
# Discard molecules without labels
mols = [ mol for mol in mols if oechem.OEHasSDData(mol, 'intermediate') ]
print(f'{len(mols)} molecules remain after discarding unlabeled molecules')

# Sort molecules by IC50
print('Sorting molecules by fluorescence IC50...')
def ic50(mol):
    value = oechem.OEGetSDData(mol, 'f_avg_IC50')
    try:
        return float(value)
    except Exception as e:
        return 100

mols.sort(key=ic50)

# Write molecules
output_filename = 'sorted/sprint-5-retrospective.csv'
with oechem.oemolostream(output_filename) as ofs:
    for mol in track(mols, description='Writing molecules...'):
        oechem.OEWriteMolecule(ofs, mol)

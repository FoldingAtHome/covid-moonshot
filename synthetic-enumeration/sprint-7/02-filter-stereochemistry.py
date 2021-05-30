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

print('reading compound designs...')
source_filenames = [
    'sorted/sprint-7-2021-05-11.csv'
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

# Expand uncertain stereochemistry
print('Enumerating undefiend stereochemistry...')
from openeye import oeomega
omegaOpts = oeomega.OEOmegaOptions()
omegaOpts.SetMaxConfs(1)
omegaOpts.SetWarts(True)
omegaOpts.SetStrictStereo(False)
omega = oeomega.OEOmega(omegaOpts)
flipperOpts = oeomega.OEFlipperOptions()
flipperOpts.SetEnumSpecifiedStereo(False)
enumerated_mols = list()
for mol in track(mols):
    for enantiomer in oeomega.OEFlipper(mol, flipperOpts):
        enantiomer = oechem.OEMol(enantiomer)
        ret_code = omega.Build(enantiomer)
        if ret_code == oeomega.OEOmegaReturnCode_Success:
            enumerated_mols.append(enantiomer)
        else:
            oechem.OEThrow.Warning("%s: %s" %
                (enantiomer.GetTitle(), oeomega.OEGetOmegaError(ret_code)))
mols = enumerated_mols
print(f'{len(mols)} molecules after enumeration')

# Annotate molecules with SMARTS labels
print('Annotating SMARTS labels...')
import csv
labels_filename = 'annotations/annotations-stereo.csv' # list of labels for various SMARTS patterns
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

# Write molecules
output_filename = 'sorted/sprint-7-2021-05-11-stereofilter.csv'
with oechem.oemolostream(output_filename) as ofs:
    for mol in track(mols, description='Writing molecules...'):
        oechem.OEWriteMolecule(ofs, mol)

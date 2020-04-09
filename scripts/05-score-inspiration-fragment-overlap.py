"""
Score overlap of docked ligands with X-ray fragments

"""

# Make a list of all fragment sites to dock to
noncovalent_active_site_fragments = ['x0072', 'x0104', 'x0107', 'x0161', 'x0195', 'x0305', 'x0354', 'x0387', 'x0395', 'x0397', 'x0426', 'x0434', 'x0540',
    'x0678', 'x0874', 'x0946', 'x0967', 'x0991', 'x0995', 'x1077', 'x1093', 'x1249']
covalent_active_site_fragments = ['x0689', 'x0691', 'x0692', 'x0705','x0708', 'x0731', 'x0734', 'x0736', 'x0749', 'x0752', 'x0755', 'x0759',
'x0769', 'x0770', 'x0771', 'x0774', 'x0786', 'x0820', 'x0830', 'x0831', 'x0978', 'x0981', 'x1308', 'x1311', 'x1334', 'x1336', 'x1348',
'x1351', 'x1358', 'x1374', 'x1375', 'x1380', 'x1382', 'x1384', 'x1385', 'x1386', 'x1392', 'x1402', 'x1412', 'x1418', 'x1425', 'x1458',
'x1478', 'x1493']
active_site_fragments = noncovalent_active_site_fragments + covalent_active_site_fragments
dimer_interface_fragments = ['x0887', 'x1187']
all_fragments = noncovalent_active_site_fragments + covalent_active_site_fragments + dimer_interface_fragments

if __name__ == '__main__':
    from openeye import oechem, oeshape
    import csv
    from tqdm import tqdm
    from glob import glob
    import os, re

    # Parse arguments
    import argparse

    parser = argparse.ArgumentParser(description='Aggregate results and coalesce Folding@Home RUNs.')
    parser.add_argument('--docked', dest='docked_molecules', type=str, default='covid_submissions_03_31_2020-docked.sdf',
                        help='structures of docked molecules (default: covid_submissions_03_31_2020-docked)')
    parser.add_argument('--output', dest='output_prefix', type=str, default='covid_submissions_03_31_2020-overlap',
                        help='output aggregated CSV file (default: covid_submissions_03_31_2020-overlap.csv)')
    parser.add_argument('--receptors', dest='receptor_basedir', type=str, default='../receptors/monomer',
                        help='directory of receptor conformations (default: ../receptors/monomer)')
    parser.add_argument('--clean', dest='clean', action='store_true', default=False,
                        help='if specified, will only store minimal information for each molecule (default: False)')
    parser.add_argument('--sort', dest='sort', action='store_true', default=False,
                        help='if specified, will sort according to overlap (default: False)')
    parser.add_argument('--covalent', dest='covalent', action='store_true', default=False,
                        help='if specified, will only consider those with `covalent_warhead=True` (default: False)')

    args = parser.parse_args()

    # Read the docked molecules as CSV
    print('Loading molecules...')
    docked_molecules = list()
    with oechem.oemolistream(args.docked_molecules) as ifs:
        docked_molecules = list()
        molecule = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, molecule):
            oechem.OESuppressHydrogens(molecule)
            docked_molecules.append( molecule.CreateCopy() )
    print(f'{len(docked_molecules)} read')

    if args.covalent:
        print('Only filtering covalent fragments')
        docked_molecules = [molecule for molecule in docked_molecules if oechem.OEGetSDData(molecule, 'covalent_warhead')=='TRUE']
        print(f'{len(docked_molecules)} remain after filtering')

    print('Loading fragments')
    filenames = glob(os.path.join(args.receptor_basedir, 'Mpro-x*-ligand.mol2'))
    fragments = dict()
    for filename in tqdm(filenames):
        with oechem.oemolistream(filename) as ifs:
            fragment = oechem.OEGraphMol()
            while oechem.OEReadMolecule(ifs, fragment):
                oechem.OESuppressHydrogens(fragment)
                match = re.search('Mpro-(?P<fragment_name>x\d\d\d\d)-ligand.mol2', filename)
                fragment_name = match.group('fragment_name')
                # Set fragment name as title
                fragment.SetTitle(fragment_name)

                # Store it
                fragments[fragment_name] = fragment.CreateCopy()
    print(f'{len(fragments)} fragments loaded.')

    # Get appropriate function to calculate analytic shape
    prep = oeshape.OEOverlapPrep()
    result = oeshape.OEOverlapResults()
    fragment_func = oeshape.OEExactColorFunc()
    print('Computing overlap scores...')
    for molecule in tqdm(docked_molecules):
        # Compute overlap with merged query
        refmol = molecule.CreateCopy()
        prep.Prep(refmol)
        fragment_func.SetupRef(refmol)
        volume = oeshape.OECalcVolume(molecule)

        # Compute overlaps
        overlapping_fragments = list()
        fragment_overlap_scores = dict()
        score = 0.0
        inspiration_fragment_names = oechem.OEGetSDData(molecule, 'fragments').split(',')
        for fragment_name in inspiration_fragment_names:
            if fragment_name not in fragments:
                continue
            fragment = fragments[fragment_name]
            fitmol = fragment.CreateCopy()
            prep.Prep(fitmol)
            fragment_func.Overlap(fitmol, result)
            # Compute overlap (fraction of the fragment covered)
            fragment_overlap = result.GetRefTverskyCombo()
            # Store fragment
            fragment_overlap_scores[fragment_name] = fragment_overlap
            overlapping_fragments.append(fragment_name)
            score += fragment_overlap

        overlapping_fragments.sort(key=lambda fragment_name : -fragment_overlap_scores[fragment_name])

        # Attach scores
        oechem.OESetSDData(molecule, 'overlapping_fragments', ','.join(overlapping_fragments))
        oechem.OESetSDData(molecule, 'overlap_score', str(score))
        oechem.OESetSDData(molecule, 'volume', str(volume))

    def overlap_score(molecule):
        overlap_score = float(oechem.OEGetSDData(molecule, 'overlap_score'))
        return -overlap_score

    if args.sort:
        print('Sorting by overlap score...')
        docked_molecules.sort(key=overlap_score)

    # Write molecules
    print('Writing molecules...')
    for ext in ['sdf', 'csv', 'pdb', 'mol2']:
        filename = f'{args.output_prefix}.{ext}'
        with oechem.oemolostream(filename) as ofs:
            print(filename)
            for molecule in tqdm(docked_molecules):
                if args.clean:
                    for sdpair in oechem.OEGetSDDataPairs(molecule):
                        if sdpair.GetTag() not in ['Hybrid2', 'fragments', 'site', 'number_of_overlapping_fragments', 'overlapping_fragments', 'overlap_score', 'volume', 'covalent_distance_min', 'covalent_distance_mean', 'covalent_distance_stderr']:
                            oechem.OEDeleteSDData(molecule, sdpair.GetTag())
                oechem.OEWriteMolecule(ofs, molecule)

    # Writing molecules with overlapping fragments for inspection
    print('Writing molecules with overlapping fragments for inspection...')
    basedir = f'{args.output_prefix}'
    os.makedirs(basedir, exist_ok=True)
    for index, molecule in tqdm(enumerate(docked_molecules)):
        filename = os.path.join(basedir, f'{index:05d}-{molecule.GetTitle()}.pdb')
        with oechem.oemolostream(filename) as ofs:
            oechem.OEWriteMolecule(ofs, molecule)
            overlapping_fragments = oechem.OEGetSDData(molecule, 'overlapping_fragments').split(',')
            for fragment in overlapping_fragments:
                if fragment != '':
                    oechem.OEWriteMolecule(ofs, fragments[fragment])

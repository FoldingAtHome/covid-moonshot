"""
Score overlap of docked ligands with X-ray fragments

"""

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

    # Create merged query from fragments
    from openeye import oegrid, oeshape
    prep = oeshape.OEOverlapPrep()
    query = oeshape.OEShapeQuery()
    query_molecule = oechem.OEGraphMol()
    for fragment_name, fragment in fragments.items():
        fragment_copy = fragment.CreateCopy()
        prep.Prep(fragment_copy)
        for atom in oeshape.OEGetColorAtoms(fragment_copy):
            coords = oechem.OEFloatArray(3)
            fragment_copy.GetCoords(atom, coords)
            prefactor = 1.0
            width = 1.0
            gauss = oegrid.OEGaussian(prefactor, width, coords, oeshape.OEGetColorType(atom))
            query.AddColorGaussian(gauss)
            #query_molecule.NewAtom(atom.GetAtomicNum())

    query.SetMolecule(fragment_copy) # TODO: Can this be any molecule?
    merged_func = oeshape.OEOverlapFunc()
    merged_func.SetupRef(query)

    # Get appropriate function to calculate analytic shape
    result = oeshape.OEOverlapResults()
    fragment_func = oeshape.OEOverlapFunc()
    print('Computing overlap scores...')
    OVERLAP_THRESHOLD = 0.60
    for molecule in tqdm(docked_molecules):
        # Compute overlap with merged query
        fitmol = molecule.CreateCopy()
        prep.Prep(fitmol)
        merged_func.Overlap(fitmol, result)
        score = result.GetRefTverskyCombo()

        # Compute overlaps
        fragment_func.SetupRef(molecule)
        overlapping_fragments = list()
        for fragment_name, fragment in fragments.items():
            #overlap, volume = compute_fragment_overlap(molecule, fragment)
            fragment_func.Overlap(fragment, result)
            # Compute overlap (fraction of the fragment covered)
            fragment_overlap = result.GetFitTverskyCombo()
            # Store fragment
            if fragment_overlap > OVERLAP_THRESHOLD:
                overlapping_fragments.append(fragment_name)

        n_overlapping_fragments = len(overlapping_fragments)

        # Attach scores
        oechem.OESetSDData(molecule, 'number_of_overlapping_fragments', str(n_overlapping_fragments))
        oechem.OESetSDData(molecule, 'overlapping_fragments', ','.join(overlapping_fragments))
        oechem.OESetSDData(molecule, 'overlap_score', str(score))

    print('Sorting...')
    def overlap_score(molecule):
        overlap_score = float(oechem.OEGetSDData(molecule, 'overlap_score'))
        return -overlap_score

    def nfragments_score(molecule):
        overlapping_fragments = oechem.OEGetSDData(molecule, 'overlapping_fragments')
        if overlapping_fragments == '':
            return 0
        else:
            return -len(overlapping_fragments.split(','))

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
                        if sdpair.GetTag() not in ['Hybrid2', 'fragments', 'site', 'number_of_overlapping_fragments', 'overlapping_fragments', 'overlap_score']:
                            oechem.OEDeleteSDData(molecule, sdpair.GetTag())
                oechem.OEWriteMolecule(ofs, molecule)

    # Writing molecules with overlapping fragments for inspection
    print('Writing molecules with overlapping fragments for inspection...')
    basedir = f'{args.output_prefix}'
    os.makedirs(basedir, exist_ok=True)
    for index, molecule in tqdm(enumerate(docked_molecules)):
        if nfragments_score(molecule) != 0:
            filename = os.path.join(basedir, f'{index:05d}-{molecule.GetTitle()}.pdb')
            with oechem.oemolostream(filename) as ofs:
                oechem.OEWriteMolecule(ofs, molecule)
                overlapping_fragments = oechem.OEGetSDData(molecule, 'overlapping_fragments').split(',')
                for fragment in overlapping_fragments:
                    if fragment != '':
                        oechem.OEWriteMolecule(ofs, fragments[fragment])

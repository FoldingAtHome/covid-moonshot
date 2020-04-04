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
    parser.add_argument('--docked', dest='docked_basedir', type=str, default='covid_submissions_03_31_2020',
                        help='base directory for docked molecules (default: covid_submissions_03_31_2020/)')
    parser.add_argument('--output', dest='output_filename', type=str, default='covid_submissions_03_31_2020-overlap.csv',
                        help='output aggregated CSV file (default: covid_submissions_03_31_2020-overlap.csv)')
    parser.add_argument('--receptors', dest='receptor_basedir', type=str, default='../receptors/monomer',
                        help='directory of receptor conformations (default: ../receptors/monomer)')

    args = parser.parse_args()

    # Read the docked molecules as CSV
    print('Loading molecules...')
    docked_filenames = glob(f'{args.docked_basedir}/docking/*ligand.sdf')
    docked_molecules = list()
    molecule = oechem.OEGraphMol()
    for filename in tqdm(docked_filenames):
        with oechem.oemolistream(filename) as ifs:
            molecules = list()
            while oechem.OEReadMolecule(ifs, molecule):
                molecules.append( molecule.CreateCopy() )
        docked_molecules += molecules
    print(f'{len(docked_molecules)} read')

    print('Loading fragments')
    filenames = glob(os.path.join(args.receptor_basedir, 'Mpro-x*-ligand.mol2'))
    fragments = dict()
    for filename in tqdm(filenames):
        with oechem.oemolistream(filename) as ifs:
            fragment = oechem.OEGraphMol()
            while oechem.OEReadMolecule(ifs, fragment):
                match = re.search('Mpro-(?P<fragment_name>x\d\d\d\d)-ligand.mol2', filename)
                fragment_name = match.group('fragment_name')
                # Set fragment name as title
                fragment.SetTitle(fragment_name)
                # Store it
                fragments[fragment_name] = fragment.CreateCopy()

    # Compute distinct fragments
    compute_fragment_basis = False
    if compute_fragment_basis:
        # TODO: Use PCCA+ when we have time to implement this
        print('Extract a fragment basis...')
        import numpy as np
        from sklearn.cluster import SpectralClustering
        fragment_names = list(fragments.keys())
        nfragments = len(fragment_names)
        n_clusters = 30
        affinity_matrix = np.ones([nfragments, nfragments], np.float32)
        for i in range(nfragments):
            shapeFunc = oeshape.OEAnalyticShapeFunc()
            shapeFunc.SetupRef(fragments[fragment_names[i]])
            result = oeshape.OEOverlapResults()
            for j in range(i+1,nfragments):
                shapeFunc.Overlap(fragments[fragment_names[j]], result)
                overlap = result.GetTanimoto()
                affinity_matrix[i,j] = overlap
                affinity_matrix[j,i] = overlap
        clustering = SpectralClustering(n_clusters=n_clusters, affinity='precomputed').fit(affinity_matrix)
        unique_fragment_names = [ index for index in range(n_clusters) ]
        for fragment_index, cluster_index in enumerate(clustering.labels_):
            unique_fragment_names[cluster_index] = fragment_names[fragment_index]
        fragments = { fragment_name : fragments[fragment_name] for fragment_name in unique_fragment_names }

    print('Computing overlap scores...')
    OVERLAP_THRESHOLD = 0.4
    molecules = list()
    directories = glob('Files/x*')
    max_overlap = 0.0
    for directory in tqdm(directories):
        _, docked_fragment = os.path.split(directory)
        with oechem.oemolistream(os.path.join(directory, 'poses.mol2')) as ifs:
            molecule = oechem.OEGraphMol()
            index = 1
            while oechem.OEReadMolecule(ifs, molecule):
                # Correct molecule title to CID
                dockvalent_id, _ = molecule.GetTitle().split()
                covalent_id, _ = dockvalent_id.split('_')
                refmol = molecule_index[covalent_id]
                molecule.SetTitle(refmol.GetTitle())
                oechem.OESetSDData(molecule, 'dockvalent_id', dockvalent_id)

                # Get appropriate function to calculate analytic shape
                shapeFunc = oeshape.OEAnalyticShapeFunc()
                shapeFunc.SetupRef(molecule)

                # Compute overlaps
                score = 0.0
                result = oeshape.OEOverlapResults()
                overlapping_fragments = list()
                for fragment_name, fragment in fragments.items():
                    #overlap, volume = compute_fragment_overlap(molecule, fragment)
                    shapeFunc.Overlap(fragment, result)
                    # Compute overlap
                    overlap = result.GetTanimoto()
                    if overlap > max_overlap:
                        max_overlap = overlap
                    # Compute volume
                    volume = oeshape.OECalcVolume(fragment)
                    # Accumulate score
                    score += volume * overlap
                    # Store fragment
                    if overlap > OVERLAP_THRESHOLD:
                        overlapping_fragments.append(fragment_name)
                    # DEBUG
                    #print(f'  {fragment_name} {overlap:8.3f} {volume:10.3f}')

                    # Correct for fragment-fragment overlap


                n_overlapping_fragments = len(overlapping_fragments)
                #print(f'{directory:16} {index:6} {molecule.GetTitle():18} {covalent_id:8} {score:10.3f} {overlapping_fragments}')

                # Attach scores
                oechem.OESetSDData(molecule, 'number_of_overlapping_fragments', str(n_overlapping_fragments))
                oechem.OESetSDData(molecule, 'overlapping_fragments', ','.join(overlapping_fragments))
                oechem.OESetSDData(molecule, 'overlap_score', str(score))
                oechem.OESetSDData(molecule, 'docked_fragment', docked_fragment)

                # Store a copy
                molecules.append( molecule.CreateCopy() )

                index += 1

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

    #molecules.sort(key=overlap_score)
    molecules.sort(key=overlap_score)
    print(f'max overlap: {max_overlap}')

    # Write molecules
    print('Writing molecules...')
    for ext in ['sdf', 'mol2', 'pdb', 'csv']:
        filename = f'covalent-docking-overlap.{ext}'
        with oechem.oemolostream(filename) as ofs:
            print(filename)
            for molecule in tqdm(molecules):
                oechem.OEWriteMolecule(ofs, molecule)

    # Writing molecules with overlapping fragments for inspection
    print('Writing molecules with overlapping fragments for inspection...')
    basedir = 'docked-with-overlapping-fragments'
    os.makedirs(basedir, exist_ok=True)
    for index, molecule in tqdm(enumerate(molecules)):
        if nfragments_score(molecule) != 0:
            filename = os.path.join(basedir, f'{index:05d}-{molecule.GetTitle()}-{dockvalent_id}.pdb')
            with oechem.oemolostream(filename) as ofs:
                oechem.OEWriteMolecule(ofs, molecule)
                overlapping_fragments = oechem.OEGetSDData(molecule, 'overlapping_fragments').split(',')
                for fragment in overlapping_fragments:
                    if fragment != '':
                        oechem.OEWriteMolecule(ofs, fragments[fragment])

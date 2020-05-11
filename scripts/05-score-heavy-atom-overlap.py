"""
Score overlap of docked ligands with X-ray fragments

This script computes the following:
* A continuous `overlap_score` that roughly measures the number of heavy atoms for the docked molecule that overlap with any fragment heavy atoms
* A list of `overlapping_fragments` that contains a minimal set of distinct fragments spanned by the docked molecule
* The `volume` of the docked molecule (in A^3)

This should be invoked in `moonshot-compounds/` with:

$ PREFIX=covid_submissions_all_info
$ rm -rf covid_submissions_all_info-docked-overlap
$ python -i ../scripts/05-score-heavy-atom-overlap.py --docked $PREFIX-docked.sdf --output $PREFIX-docked-overlap --sort

This will produce several sets of files sorted from most overlapping fragments to fewest
* `$PREFIX-docked-overlap.csv` : summary CSV file
* `$PREFIX-docked-overlap.sdf` : annotated SDF file with original structures
* `$PREFIX-docked-overlap.pdb` : same as input structures
* `$PREFIX-docked-overlap.mol2` : same as input structures
* `$PREFIX-docked-overlap/` : directory with files numbered in order of most overlapping fragments to fewest
  * `/00000-MAK-UNK-9e4a73aa-2.pdb` - PDB file containing the docked molecule (first molecule) and overlapping fragments (subsequent molecules)
  * ...

Conda package requirements (see `environment.yml`):
* numpy
* scipy
* tqdm
* openeye-toolkits

"""

def greedy_set_cover(list_of_subsets, min_increment=1):
    """At each step, choose the subset with the greatest number of uncovered elements

    This is essentially as good as you can do
        https://en.wikipedia.org/wiki/Set_cover_problem#Greedy_algorithm

    @author Josh Fass

    Parameters
    ----------
    list_of_subsets : list of set
        List of subsets
    min_increment : int, optional, default=1
        minimum acceptable increment

    Returns
    -------
    included_indices : list of int
        List of indices of selected subsets

    """
    # which elements have we included so far
    running_union = set()

    target = set.union(*list_of_subsets)

    # which subsets have we included so far
    included_indices = []
    remaining_indices = list(range(len(list_of_subsets)))

    # at each step, greedily add whichever subset would add the most elements to running_union
    improvable = True
    while (len(remaining_indices) > 0) and (improvable):
        weights = [len(running_union.union(list_of_subsets[i])) for i in remaining_indices]
        if max(weights) == len(target):
            improvable = False

        if not max(weights) >= len(running_union) + min_increment:
            improvable = False
            continue

        i = remaining_indices[np.argmax(weights)]
        running_union.update(list_of_subsets[i])
        included_indices.append(i)
        remaining_indices.remove(i)

    # return the indices of the selected subsets
    return included_indices

if __name__ == '__main__':
    from openeye import oechem, oeshape
    import csv
    from tqdm import tqdm
    from glob import glob
    import os, re
    import numpy as np
    import scipy

    # Parse arguments
    import argparse

    parser = argparse.ArgumentParser(description='Score fragment poses by how many heavy atoms overlap with unique fragment heavy atoms')
    parser.add_argument('--docked', dest='docked_molecules', type=str, default='covid_submissions_all_info-docked.sdf',
                        help='structures of docked molecules (default: covid_submissions_all_info-docked)')
    parser.add_argument('--receptors', dest='receptor_basedir', type=str, default='../receptors/monomer',
                        help='directory of receptor conformations (default: ../receptors/monomer)')
    parser.add_argument('--output', dest='output_prefix', type=str, default='covid_submissions_all_info-docked-overlap',
                        help='output aggregated CSV file (default: covid_submissions_all_info-docked-overlap.csv)')
    parser.add_argument('--clean', dest='clean', action='store_true', default=False,
                        help='if specified, will only store minimal information for each molecule (default: False)')
    parser.add_argument('--sort', dest='sort', action='store_true', default=False,
                        help='if specified, will sort according to overlap (default: False)')

    args = parser.parse_args()

    # Read the docked molecules as CSV
    print('Loading molecules and suppressing hydrogens...')
    docked_molecules = list()
    with oechem.oemolistream(args.docked_molecules) as ifs:
        docked_molecules = list()
        molecule = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, molecule):
            oechem.OESuppressHydrogens(molecule)
            docked_molecules.append( molecule.CreateCopy() )
    print(f'{len(docked_molecules)} read')

    # Read fragments
    print('Loading fragments and suppressing hydrogens...')
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
    nfragments = len(fragments)
    print(f'{nfragments} fragments loaded.')

    # Extract fragment heavy atom positions
    print('Extracting fragment positions...')
    fragment_positions = dict()
    heavy_atom_positions = list() # concatenated positions to be clustered
    fragment_indices = list() # fragment_indices[atom_index] is the fragment index of atom atom_index
    fragment_names = list()
    fragment_index = 0
    for fragment_name, fragment in tqdm(fragments.items()):
        # Extract positions
        positions = np.zeros([fragment.NumAtoms(),3], np.float32)
        for index, atom in enumerate(fragment.GetAtoms()):
            coords = oechem.OEFloatArray(3)
            fragment.GetCoords(atom, coords)
            positions[index,:] = coords[:]
            # Store concatenated positions and fragment identities
            heavy_atom_positions.append(coords)
            fragment_indices.append(fragment_index)

        fragment_names.append(fragment_name)
        fragment_index += 1
        # Store positions
        fragment_positions[fragment_name] = positions

    heavy_atom_positions = np.array(heavy_atom_positions)
    nheavy = heavy_atom_positions.shape[0]
    fragment_indices = np.array(fragment_indices)

    # Cluster heavy atom positions
    print('Clustering heavy atoms...')
    from scipy.cluster.hierarchy import fclusterdata
    CUTOFF = 0.60 # Angstroms : TODO: Change this to be twice the stddev distance
    cluster_labels = fclusterdata(heavy_atom_positions, t=CUTOFF, criterion='distance')
    nclusters = cluster_labels.max()
    print(f'There are {nclusters} clusters')
    bins = np.arange(nclusters+1)+0.5
    cluster_counts, bins = np.histogram(cluster_labels, bins=bins)
    weight_j = np.array([ 1.0 / cluster_counts[cluster-1] for cluster in cluster_labels ])

    # Write out cluster centers
    print('Writing clustered-heavy-atoms.pdb')
    with open('clustered-heavy-atoms.pdb', 'wt') as outfile:
        cluster_exemplars = np.zeros([nclusters,3], np.float32)
        for index in range(nheavy):
            cluster_index = cluster_labels[index] - 1
            cluster_exemplars[cluster_index,:] = heavy_atom_positions[index,:]

        for cluster_index in range(nclusters):
            #       ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
            x, y, z = cluster_exemplars[cluster_index,:]
            line = f'ATOM  {index+1:5}  C   XXX X {1:4d} {x:8.3f}{y:8.3f}{z:8.3f}\n'
            outfile.write(line)

    # Get appropriate function to calculate analytic shape
    print('Computing overlap scores...')
    for molecule in tqdm(docked_molecules):
        # Extract positions
        positions = np.zeros([molecule.NumAtoms(),3], np.float32)
        for index, atom in enumerate(molecule.GetAtoms()):
            coords = oechem.OEFloatArray(3)
            molecule.GetCoords(atom, coords)
            positions[index,:] = coords[:]

        # Compute volume
        volume = oeshape.OECalcVolume(molecule)

        # Compute overlap score as the total number of overlapping heavy atoms
        from scipy.spatial.distance import cdist
        distances_ij = cdist(positions, heavy_atom_positions, 'euclidean')
        SIGMA = CUTOFF / 2.0
        phi_ij = np.exp(-0.5 * (distances_ij / SIGMA)**2)
        score = np.einsum('ij,j->', phi_ij, weight_j)

        # Compute set coverings for all fragments
        list_of_subsets = list()
        for fragment_index in range(nfragments):
            # Determine which atoms in the docked molecule this fragment overlaps with
            heavy_atom_indices = np.where(fragment_indices==fragment_index)
            overlapping_atoms = { index for index in range(molecule.NumAtoms()) if (distances_ij[index,heavy_atom_indices].min() < CUTOFF) }
            list_of_subsets.append(overlapping_atoms)
        # Identify the smallest number of fragments that provides the maximum possible covering of docked molecule atoms
        overlapping_fragment_indices = greedy_set_cover(list_of_subsets, min_increment=4)
        overlapping_fragments = [ fragment_names[fragment_index] for fragment_index in overlapping_fragment_indices ]
        n_overlapping_fragments = len(overlapping_fragments)

        # Attach scores
        oechem.OESetSDData(molecule, 'number_of_overlapping_fragments', str(n_overlapping_fragments))
        oechem.OESetSDData(molecule, 'overlapping_fragments', ','.join(overlapping_fragments))
        oechem.OESetSDData(molecule, 'overlap_score', str(score))
        oechem.OESetSDData(molecule, 'volume', str(volume))


    def overlap_score(molecule):
        overlap_score = float(oechem.OEGetSDData(molecule, 'overlap_score'))
        return -overlap_score

    def nfragments_score(molecule):
        overlapping_fragments = oechem.OEGetSDData(molecule, 'overlapping_fragments')
        if overlapping_fragments == '':
            return 0
        else:
            return -len(overlapping_fragments.split(','))

    if args.sort:
        print('Sorting by overlap score...')
        #docked_molecules.sort(key=overlap_score)
        docked_molecules.sort(key=nfragments_score)


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

"""
Combine all docking results

Prepare SDF file in format for Fragalysis upload Spec 1.2:
https://discuss.postera.ai/t/providing-computed-poses-for-others-to-look-at/1155/8?u=johnchodera

"""

# Read ligands
def read_csv_molecules(filename):
    """Read molecules from the specified path

    Parameters
    ----------
    filename : str
        File from which molecules are to be read

    Returns
    -------
    molecules : list of openeye.oechem.OEMol
        The read molecules
    """

    from openeye import oechem
    mol = oechem.OEMol()
    molecules = list()
    with oechem.oemolistream(filename) as ifs:
        while oechem.OEReadCSVFile(ifs, mol):
            molecules.append(oechem.OEMol(mol))
    return molecules

def read_original_smiles(filename):
    """Read mapping of original SMILES for each compound name

    Parameters
    ----------
    filename : str
        File from which molecules are to be read

    Returns
    -------
    original_smiles : dict
        original_smiles[compound_id] is the original SMILES
    """
    original_smiles = dict()
    with open(filename, 'r') as infile:
        for line in infile:
            elements = line.split(',')
            smiles = elements[0]
            cid = elements[1]
            original_smiles[cid] = smiles

    return original_smiles

def score(molecule, field='Hybrid2'):
    """Return the docking score"""
    from openeye import oechem
    value = oechem.OEGetSDData(molecule, field)
    return float(value)

if __name__ == '__main__':
    # Parse arguments
    import argparse

    parser = argparse.ArgumentParser(description='Aggregate docking results')
    parser.add_argument('--molecules', dest='molecules_filename', type=str, default='covid_submissions_all_info.csv',
                        help='molecules CSV file to pull from (default: covid_submissions_all_info-2020-04-23.csv)')
    parser.add_argument('--docked', dest='docked_basedir', type=str, default='covid_submissions_all_info-docked',
                        help='base directory for docked molecules (default: covid_submissions_all_info-docked/)')
    parser.add_argument('--output', dest='output_filename', type=str, default='covid_submissions_all_info-docked.csv',
                        help='output aggregated CSV file (default: covid_submissions_all_info.csv)')
    parser.add_argument('--clean', dest='clean', action='store_true', default=False,
                        help='if specified, will only store minimal information for each molecule (default: False)')
    parser.add_argument('--fragalysis', dest='fragalysis', default=False,
                        help='if specified, URL to use in formatting data for fragalysis; must be SDF output file (default: False)')

    args = parser.parse_args()

    import os
    from openeye import oechem
    from tqdm import tqdm
    from glob import glob

    # Read the docked molecules from SDF
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

    # Sort molecules
    docked_molecules.sort(key=score)

    # Read original SMILES to annotate molecules
    if args.molecules_filename:
        original_smiles = read_original_smiles(args.molecules_filename)

    # Clean up SD tags, if requested
    if args.clean:
        print('Cleaning SD tags...')
        allowed_tags = ['Hybrid2', 'docked_fragment', 'fragments', 'covalent_distance_min', 'covalent_distance_mean', 'covalent_distance_stddev']
        for molecule in tqdm(docked_molecules):
            for sdpair in oechem.OEGetSDDataPairs(molecule):
                tag = sdpair.GetTag()
                value = sdpair.GetValue()
                if tag not in allowed_tags:
                    oechem.OEDeleteSDData(molecule, sdpair.GetTag())

    # Translate SD tags, if requested
    if args.fragalysis:
        print('Cleaning SD tags for fragalysis...')
        descriptions = {
            'Chemgauss4' : 'The [Chemgauss4 docking score](https://docs.eyesopen.com/toolkits/cpp/dockingtk/scoring.html#section-scoring-chemgauss4) of the best scoring pose and protonation/tautomeric state to any Mpro structure; more negative is better, and fragments score ~ -10',
            'CYS145-warhead dist minimum (A)' : 'the minimum distance (in Angstroms) between CYS145 SG and any covalent warhead heavy atom during 10 ns simulation, where distances < 4A indicate the warhead is potentially well-positioned for irreversible binding (covalent inhibitors only)',
            'CYS145-warhead dist mean (A)' : 'the mean distance (in Angstroms) between CYS145 SG and any covalent warhead heavy atom during 10 ns simulation, where distances < 4A indicate the warhead is potentially well-positioned for irreversible binding (covalent inhibitors only)',
            'CYS145-warhead dist stddev (A)' : 'is the standard deviation (in Angstroms) of the distance between CYS145 SG and any covalent warhead heavy atom during 10 ns simulation, where large values may indicate the warhead samples a variety of distances from CYS145 (covalent inhibitors only)',        
            'ref_mols' : 'a comma separated list of the fragments that inspired the design of the new molecule (codes as they appear in fragalysis - e.g. x0104_0,x0692_0)',
            'ref_pdb' : 'The name of the fragment (and corresponding Mpro fragment structure) with the best scoring hybrid docking pose',
            'original SMILES' : 'the original SMILES of the compound before any computation was carried out',            
        }

        for molecule in tqdm(docked_molecules):
            for sdpair in oechem.OEGetSDDataPairs(molecule):
                tag = sdpair.GetTag()
                value = sdpair.GetValue()
                # Remove hydrogens
                oechem.OESuppressHydrogens(molecule, True)
                # Remap SD tags
                if tag == 'Hybrid2':
                    oechem.OESetSDData(molecule, 'Chemgauss4', value)
                    oechem.OEDeleteSDData(molecule, tag)
                elif tag == 'fragments':
                    value = ','.join([ f'{fragment}_0' for fragment in value.split(',') ])
                    oechem.OESetSDData(molecule, 'ref_mols', value)
                    oechem.OEDeleteSDData(molecule, tag)
                elif tag == 'docked_fragment':
                    oechem.OESetSDData(molecule, 'ref_pdb', value + '_0')
                    oechem.OEDeleteSDData(molecule, tag)
                elif tag == 'covalent_distance_min':
                    oechem.OESetSDData(molecule, 'CYS145-warhead dist minimum (A)', value)
                    oechem.OEDeleteSDData(molecule, tag)
                elif tag == 'covalent_distance_mean':
                    oechem.OESetSDData(molecule, 'CYS145-warhead dist mean (A)', value)
                    oechem.OEDeleteSDData(molecule, tag)
                elif tag == 'covalent_distance_stddev':
                    oechem.OESetSDData(molecule, 'CYS145-warhead dist stddev (A)', value)
                    oechem.OEDeleteSDData(molecule, tag)
                # Add original SMILES
                oechem.OESetSDData(molecule, 'original SMILES', original_smiles[molecule.GetTitle()])

        # Add initial blank molecule (that includes distances)
        import copy
        from datetime import datetime    
        # Find a molecule that includes distances, if present
        def has_dist(mol):
            for sdpair in oechem.OEGetSDDataPairs(mol):
                tag = sdpair.GetTag()
                if 'dist' in tag:
                    return True
            return False
        molecule = molecules[0].CreateCopy()
        for molecule in molecules:
            if has_dist(molecule):
                molecule = molecule.CreateCopy()
        # Add descriptions
        for sdpair in oechem.OEGetSDDataPairs(molecule):
            tag = sdpair.GetTag()
            value = sdpair.GetValue()
            oechem.OESetSDData(molecule, tag, descriptions[tag])
        # Add other fields
        molecule.SetTitle('ver_1.2')
        oechem.OESetSDData(molecule, 'ref_url', args.fragalysis)
        oechem.OESetSDData(molecule, 'submitter_name', 'John D. Chodera')
        oechem.OESetSDData(molecule, 'submitter_email', 'john.chodera@choderalab.org')
        oechem.OESetSDData(molecule, 'submitter_institution', 'MSKCC')        
        oechem.OESetSDData(molecule, 'generation_date', datetime.today().strftime('%Y-%m-%d'))
        oechem.OESetSDData(molecule, 'method', 'ensemble-hybrid-oedocking')        
        docked_molecules.insert(0, molecule) # make it first molecule

    # Write sorted molecules
    print(f'Writing sorted molecules to {args.output_filename}')
    with oechem.oemolostream(args.output_filename) as ofs:        
        for molecule in tqdm(docked_molecules):
            oechem.OEWriteMolecule(ofs, molecule)



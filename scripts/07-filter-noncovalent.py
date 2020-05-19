"""
Generate ROC plot from knwon actives

"""

covalent_warhead_precedence = [
    'acrylamide', 'acrylamide_adduct', 'chloroacetamide', 'chloroacetamide_adduct', 'vinylsulfonamide', 'vinylsulfonamide_adduct', 'nitrile', 'nitrile_adduct',
]

covalent_warhead_smarts = {
    'acrylamide' : '[C;H2:1]=[C;H1]C(N)=O',
    'acrylamide_adduct' : 'NC(C[C:1]S)=O',
    'chloroacetamide' : 'Cl[C;H2:1]C(N)=O',
    'chloroacetamide_adduct' : 'S[C:1]C(N)=O',
    'vinylsulfonamide' : 'NS(=O)([C;H1]=[C;H2:1])=O',
    'vinylsulfonamide_adduct' : 'NS(=O)(C[C:1]S)=O',
    #'nitrile' : 'N#[C:1]-[*]',
    #'nitrile_adduct' : 'C-S-[C:1](=N)',
    }

def get_covalent_warhead_atom(molecule, covalent_warhead_type):
    """
    Get tagged atom index in provided tagged SMARTS string, or None if no match found.

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The molecule to search
    covalent_warhead : str
        Covalent warhead name

    Returns
    -------
    index : int or None
        The atom index in molecule of the covalent atom, or None if SMARTS does not match

    """
    smarts = covalent_warhead_smarts[covalent_warhead_type]
    qmol = oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, smarts):
        raise ValueError(f"Error parsing SMARTS '{smarts}'")
    substructure_search = oechem.OESubSearch(qmol)
    substructure_search.SetMaxMatches(1)
    matches = list()
    for match in substructure_search.Match(molecule):
        # Compile list of atom indices that match the pattern tags
        for matched_atom in match.GetAtoms():
            if(matched_atom.pattern.GetMapIdx()==1):
                return matched_atom.target.GetIdx()

    return None

def find_warheads(molecule):
    """Find covalent warheads

    Parameters
    ----------
    molecule : OEMol
      The molecule

    Returns
    -------
    warheads_found : dict
       warheads_found[warhead_type] corresponds to atom index of heavy atom that forms covalent adduct

    """
    # Identify highest priority warhead found in compound
    #print('Finding warhead atoms...')
    warheads_found = dict()
    for warhead_type in covalent_warhead_smarts.keys():
        warhead_atom_index = get_covalent_warhead_atom(molecule, warhead_type)
        if warhead_atom_index is not None:
            warheads_found[warhead_type] = warhead_atom_index
            #print(f'* Covalent warhead atom found: {warhead_type} {warhead_atom_index}')

    return warheads_found

def score(molecule, field='Hybrid2'):
    """Return the docking score"""
    from openeye import oechem
    value = oechem.OEGetSDData(molecule, field)
    return float(value)

def is_active(molecule):
    """
    """
    from openeye import oechem
    value = oechem.OEGetSDData(molecule, 'IC50 (uM)')
    try:
        IC50 = float(value)
        if (IC50 < 99.0):
            return True
    except:
        pass

    return False

def compute_ROC(docked_molecules):
    """
    Parameters
    ----------
    docked_molecules : list of oechem.OEMol
        The docked molecules

    Returns
    -------
    fraction : np.array with shape [nmolecules]
        fraction[n] is fraction of library screened at molecule index n
    fpr : np.array with shape [nmolecules]
        fpr[n] is the false positive rate at molecule index n
    tpr : np.array with shape [nmolecules]
        tpr[n] is the true positive rate at molecule index n
    """
    # Create ROC curve
    import numpy as np
    actives = np.array([is_active(molecule) for molecule in docked_molecules])
    import sklearn.metrics
    scores = np.array([score(molecule) for molecule in docked_molecules])
    fpr, tpr, thresholds = sklearn.metrics.roc_curve(actives, -scores)
    AUC = sklearn.metrics.roc_auc_score(actives, -scores)

    return fpr, tpr, AUC

def bootstrap_AUC(docked_molecules, CI=0.95, nreplicates=200):
    """
    Bootstrap the AUC

    Parameters
    ----------
    docked_molecules : list of openeye.oechem.OEMol
        The docked molecules with 'Hybrid2' scores and 'IC50 (uM)' fields
    CI : float, optional, default=0.95
        Confidence interval range
    nreplicates : int, optional, default=500
        Number of replicates

    Returns
    -------
    AUC_low, AUC_high : float
        AUC low and high ends of 95% CI
    """
    import numpy as np
    AUC_n = np.zeros([nreplicates], np.float32)
    nmolecules = len(docked_molecules)
    successful_replicates = 0
    for replicate in tqdm(range(nreplicates)):
        indices = np.random.choice(np.arange(nmolecules), nmolecules, replace=True)
        sample = [ docked_molecules[index] for index in indices ]
        try:
            sample_fpr, sample_tpr, sample_AUC = compute_ROC(sample)
            AUC_n[successful_replicates] = sample_AUC
            successful_replicates += 1
        except Exception as e:
            pass

    nreplicates = successful_replicates
    AUC_n.sort()
    AUC_low = AUC_n[int(nreplicates*(1-CI)/2)]
    AUC_high = AUC_n[int(nreplicates*(1-(1-CI)/2))]
    AUC_std = AUC_n.std()

    return AUC_low, AUC_high, AUC_std

if __name__ == '__main__':
    from openeye import oechem, oeshape
    import csv
    from tqdm import tqdm
    from glob import glob
    import os, re

    # Parse arguments
    import argparse

    parser = argparse.ArgumentParser(description='Filter out noncovalent compounds.')
    parser.add_argument('--docked', dest='docked_molecules', type=str, default='compound-set_ensemble-hybrid-oedocking.sdf',
                        help='docked assayed molecules (default: compound-set_ensemble-hybrid-oedocking.sdf)')
    parser.add_argument('--output', dest='output_filename', type=str, default='compound-set_ensemble-hybrid-oedocking-noncovalent.sdf',
                        help='noncovalent compounds (default: compound-set_ensemble-hybrid-oedocking-noncovalent.sdf)')
    args = parser.parse_args()

    # Read the docked molecules
    print('Loading molecules...')
    docked_molecules = list()
    with oechem.oemolistream(args.docked_molecules) as ifs:
        molecule = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, molecule):
            docked_molecules.append( molecule.CreateCopy() )
    nmolecules = len(docked_molecules)
    print(f'{nmolecules} read')

    # Prepare filtered subsets
    molecules = {
        'covalent': list(),
        'noncovalent' : list(),
        'all' : list()
    }
    for molecule in docked_molecules:
        # Determine if covalent
        covalent_warheads = find_warheads(molecule)
        is_covalent = True if (len(covalent_warheads) > 0) else False
        oechem.OESetSDData(molecule, 'is_covalent', str(is_covalent))
        # Add to subsets
        molecules['all'].append(molecule)
        if is_covalent:
            molecules['covalent'].append(molecule)
        else:
            molecules['noncovalent'].append(molecule)

    # Write sorted molecules
    print(f'Writing sorted molecules to {args.output_filename}')
    with oechem.oemolostream(args.output_filename) as ofs:
        for molecule in tqdm(molecules['noncovalent']):
            oechem.OEWriteMolecule(ofs, molecule)

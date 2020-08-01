#!/usr/bin/env python

"""
Generate poses for relative free energy calculations using fragment structures

"""
from openeye import oechem
import numpy as np


def GetFragments(mol, minbonds, maxbonds):
    from openeye import oegraphsim
    frags = []
    fptype = oegraphsim.OEGetFPType("Tree,ver=2.0.0,size=4096,bonds=%d-%d,atype=AtmNum,btype=Order"
                                    % (minbonds, maxbonds))

    for abset in oegraphsim.OEGetFPCoverage(mol, fptype, True):
        fragatompred = oechem.OEIsAtomMember(abset.GetAtoms())

        frag = oechem.OEGraphMol()
        adjustHCount = True
        oechem.OESubsetMol(frag, mol, fragatompred, adjustHCount)
        oechem.OEFindRingAtomsAndBonds(frag)
        frags.append(oechem.OEGraphMol(frag))

    return frags


def GetCommonFragments(mollist, frags,
                       atomexpr=oechem.OEExprOpts_DefaultAtoms,
                       bondexpr=oechem.OEExprOpts_DefaultBonds):

    corefrags = []

    from rich.progress import track
    for frag in track(frags, description='Finding common fragments'):

        ss = oechem.OESubSearch(frag, atomexpr, bondexpr)
        if not ss.IsValid():
            print('Is not valid')
            continue

        validcore = True
        for mol in mollist:
            oechem.OEPrepareSearch(mol, ss)
            validcore = ss.SingleMatch(mol)
            if not validcore:
                break

        if validcore:
            corefrags.append(frag)

    return corefrags


def GetCoreFragment(refmol, mols,
                    minbonds=3, maxbonds=200,
                    atomexpr=oechem.OEExprOpts_DefaultAtoms,
                    bondexpr=oechem.OEExprOpts_DefaultBonds):

    #print("Number of molecules = %d" % len(mols))

    frags = GetFragments(refmol, minbonds, maxbonds)
    if len(frags) == 0:
        oechem.OEThrow.Error("No fragment is enumerated with bonds %d-%d!" % (minbonds, maxbonds))

    #print("Number of fragments = %d" % len(frags))
    frags = frags[:100] # DEBUG


    commonfrags = GetCommonFragments(mols, frags, atomexpr, bondexpr)
    if len(commonfrags) == 0:
        oechem.OEThrow.Error("No common fragment is found!")

    #print("Number of common fragments = %d" % len(commonfrags))

    core = None
    for frag in commonfrags:
        if core is None or GetFragmentScore(core) < GetFragmentScore(frag):
            core = frag

    return core

def GetFragmentScore(mol):

    score = 0.0
    score += 2.0 * oechem.OECount(mol, oechem.OEAtomIsInRing())
    score += 1.0 * oechem.OECount(mol, oechem.OENotAtom(oechem.OEAtomIsInRing()))

    return score

def compute_common_substructure(refmol, newmol):
    """Determine the common heavy-atom substructure between reference molecule and new molecule.

    Parameters
    ----------
    refmol : openeye.oechem.OEMol
        The reference molecule
    newmol : openeye.oechem.OEMol
        The new molecule

    Returns
    -------
    submol : openeye.oechem.OEMol
        Common substructure containing coordinates from refmol; only one is returned
    """
    # Set up MCSS
    from openeye import oechem
    mcss = oechem.OEMCSSearch(oechem.OEMCSType_Exhaustive)
    mcss.Init(refmol, atomexpr, bondexpr)
    mcss.SetMCSFunc(oechem.OEMCSMaxAtomsCompleteCycles())
    unique = True
    # Extract the first match
    match = [m for m in mcss.Match(newmol, unique)][0]
    # Create a new molecule containing just the matched substructure
    oemol = oechem.OEMol()
    oemol.SetTitle(refmol.GetTitle() + ':' + newmol.GetTitle())
    atoms = list()
    for matchpair in match.GetAtoms():
        atom = oemol.NewAtom(matchpair.pattern)
        atoms.append(atom)
    for matchpair in match.GetBonds():
        bond = oemol.NewBond(atoms[matchpair.pattern.GetBgnIdx()], atoms[matchpair.pattern.GetEndIdx()], matchpair.pattern.GetOrder())

    return submol

def expand_stereochemistry(mols):
    """Expand stereochemistry when uncertain

    Parameters
    ----------
    mols : openeye.oechem.OEGraphMol
        Molecules to be expanded

    Returns
    -------
    expanded_mols : openeye.oechem.OEMol
        Expanded molecules
    """
    expanded_mols = list()

    from openeye import oechem, oeomega
    omegaOpts = oeomega.OEOmegaOptions()
    omega = oeomega.OEOmega(omegaOpts)
    maxcenters = 12
    forceFlip = False
    enumNitrogen = True
    warts = False
    for mol in mols:
        for enantiomer in oeomega.OEFlipper(mol, maxcenters, forceFlip, enumNitrogen, warts):
            enantiomer = oechem.OEMol(enantiomer)
            expanded_mols.append(enantiomer)

    return expanded_mols

class BumpCheck:
    def __init__(self, prot_mol, cutoff=2.0):
        self.near_nbr = oechem.OENearestNbrs(prot_mol, cutoff)
        self.cutoff = cutoff

    def count(self, lig_mol):
        bump_count = 0
        for nb in self.near_nbr.GetNbrs(lig_mol):
            if (not nb.GetBgn().IsHydrogen()) and (not nb.GetEnd().IsHydrogen()):
                bump_count += np.exp(-0.5 * (nb.GetDist() / self.cutoff)**2)
        return bump_count

def generate_restricted_conformers(receptor, fixmol, mol):
    """
    Generate and select a conformer of the specified molecule using the fixedmol

    Parameters
    ----------
    receptor : openeye.oechem.OEGraphMol
       Receptor to score poses against
    fixmol : openeye.oechem.OEGraphMol
       The reference fixed part of the molecule (common core)
    mol : openeye.oechem.OEGraphMol
       The molecule to enumerate conformers for
    """
    from openeye import oechem, oeomega

    # Create an Omega instance
    #omegaOpts = oeomega.OEOmegaOptions()
    omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Dense)

    # Set the fixed reference molecule
    omegaFixOpts = oeomega.OEConfFixOptions()
    omegaFixOpts.SetFixMaxMatch(10) # allow multiple MCSS matches
    omegaFixOpts.SetFixDeleteH(True) # only use heavy atoms
    omegaFixOpts.SetFixMol(fixmol)
    omegaOpts.SetConfFixOptions(omegaFixOpts)

    molBuilderOpts = oeomega.OEMolBuilderOptions()
    molBuilderOpts.SetStrictAtomTypes(False) # don't give up if MMFF types are not found
    omegaOpts.SetMolBuilderOptions(molBuilderOpts)

    omegaOpts.SetWarts(False) # expand molecule title
    omegaOpts.SetStrictStereo(True) # set strict stereochemistry
    omegaOpts.SetIncludeInput(False) # don't include input
    #omegaOpts.SetMaxConfs(50000) # generate lots of conformers
    #omegaOpts.SetEnergyWindow(10.0) # allow high energies
    omega = oeomega.OEOmega(omegaOpts)

    from openeye import oequacpac
    if not oequacpac.OEGetReasonableProtomer(mol):
        print('No reasonable protomer found')
        return None

    mol = oechem.OEMol(mol) # multi-conformer molecule

    ret_code = omega.Build(mol)
    if (mol.GetDimension() != 3) or (ret_code != oeomega.OEOmegaReturnCode_Success):
        print(f'Omega failure: {mol.GetDimension()} and {oeomega.OEGetOmegaError(ret_code)}')
        return None

    # Extract poses
    poses = [ pose for pose in mol.GetConfs() ]

    # Score clashes
    bump_check = BumpCheck(receptor)
    clash_scores = [ bump_check.count(pose) for pose in poses ]

    # Score docking poses
    from openeye import oedocking
    score = oedocking.OEScore(oedocking.OEScoreType_Chemgauss4)
    score.Initialize(receptor)
    docking_scores = [ score.ScoreLigand(pose) for pose in poses ]

    # Select the best docking score
    import numpy as np
    pose_index = np.argmin(docking_scores)
    mol.SetActive(poses[pose_index])
    oechem.OESetSDData(mol, 'clash_score', str(clash_scores[pose_index]))
    oechem.OESetSDData(mol, 'docking_score', str(docking_scores[pose_index]))

    # Convert to single-conformer molecule
    mol = oechem.OEGraphMol(mol)

    return mol

# TODO: import this from https://github.com/postera-ai/COVID_moonshot_submissions/blob/master/lib/utils.py
def get_series(smi):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import Descriptors
    series_SMARTS_dict = {
        # "3-aminopyridine": "[R1][C,N;R0;!$(NC(=O)CN)]C(=O)[C,N;R0;!$(NC(=O)CN)][c]1cnccc1",
        "3-aminopyridine-like": "[R1]!@[C,N]C(=O)[C,N]!@[R1]",
        "3-aminopyridine-strict": "c1ccncc1NC(=O)!@[R1]",
        "Ugi": "[c,C:1][C](=[O])[N]([c,C,#1:2])[C]([c,C,#1:3])([c,C,#1:4])[C](=[O])[NH1][c,C:5]",
        "quinolones": "NC(=O)c1cc(=O)[nH]c2ccccc12",
        "piperazine-chloroacetamide": "O=C(CCl)N1CCNCC1",
    }

    def check_if_smi_in_series(
        smi, SMARTS, MW_cutoff=550, num_atoms_cutoff=70, num_rings_cutoff=10
    ):
        mol = Chem.MolFromSmiles(smi)
        MW = Chem.Descriptors.MolWt(mol)
        num_heavy_atoms = mol.GetNumHeavyAtoms()
        num_rings = Chem.rdMolDescriptors.CalcNumRings(mol)
        patt = Chem.MolFromSmarts(SMARTS)
        if (
            (
                len(
                    Chem.AddHs(Chem.MolFromSmiles(smi)).GetSubstructMatches(
                        patt
                    )
                )
                > 0
            )
            and (MW <= MW_cutoff)
            and (num_heavy_atoms <= num_atoms_cutoff)
            and (num_rings <= num_rings_cutoff)
        ):
            return True
        else:
            return False

    for series in series_SMARTS_dict:
        series_SMARTS = series_SMARTS_dict[series]
        if series == "3-amonipyridine-like":
            if check_if_smi_in_series(
                smi,
                series_SMARTS,
                MW_cutoff=410,
                num_rings_cutoff=3,
                num_atoms_cutoff=28,
            ):
                return series
        else:
            if check_if_smi_in_series(smi, series_SMARTS):
                return series
    return None

def generate_poses(fragment, prefix, fragment_title=None, filter_series=None):
    """
    Parameters
    ----------
    fragment_title : str, optional, default=None
        If specified, replace the reference X-ray ligand title with this name
    filter_series : str, optional, default=None
        If specified, filter out only molecules matching this series name from the loaded molecules
    """
    # Read receptor
    print('Reading receptor...')
    from openeye import oechem
    receptor = oechem.OEGraphMol()
    receptor_filename = f'../receptors/monomer/Mpro-{fragment}_0_bound-receptor.oeb.gz'
    #with oechem.oemolistream(receptor_filename) as infile:
    #    oechem.OEReadMolecule(infile, receptor)
    from openeye import oedocking
    oedocking.OEReadReceptorFile(receptor, receptor_filename)
    print(f'  Receptor has {receptor.NumAtoms()} atoms.')

    # Read target molecules
    print('Reading target molecules...')
    from openeye import oechem
    targets_filename = prefix + '.csv'
    target_molecules = list()
    with oechem.oemolistream(targets_filename) as ifs:
        for mol in ifs.GetOEGraphMols():
            target_molecules.append( oechem.OEGraphMol(mol) )
    if len(target_molecules) == 0:
        raise Exception('No target molecules specified; check filename!')
    print(f'  There are {len(target_molecules)} target molecules')

    # Read reference molecule with coordinates
    refmol_filename = f'../receptors/monomer/Mpro-{fragment}_0_bound-ligand.mol2'
    refmol = None
    with oechem.oemolistream(refmol_filename) as ifs:
        for mol in ifs.GetOEGraphMols():
            refmol = mol
            break
    if refmol is None:
        raise Exception(f'Could not read {refmol_filename}')
    print(f'Reference molecule has {refmol.NumAtoms()} atoms')
    if fragment_title is not None:
        refmol.SetTitle(fragment_title)
        # Copy data from target molecles if present
        for mol in target_molecules:
            if mol.GetTitle() == fragment_title:
                print(f'{fragment_title} found in target_molecules; copying SDData')
                oechem.OECopySDData(refmol, mol)
                break

    if filter_series is not None:
        print(f'Filtering out series {filter_series}...')
        target_molecules = [ mol for mol in target_molecules if get_series(oechem.OECreateSmiString(mol)) == filter_series ]
        print(f'  There are {len(target_molecules)} target molecules')
        with oechem.oemolostream(f'filtered.mol2') as ofs:
            for mol in target_molecules:
                oechem.OEWriteMolecule(ofs, oechem.OEGraphMol(mol))

    # Expand uncertain stereochemistry
    print('Expanding uncertain stereochemistry...')
    target_molecules = expand_stereochemistry(target_molecules)
    print(f'  There are {len(target_molecules)} target molecules')

    # Get core fragment
    print('Identifying core fragment...')
    from openeye import oechem
    core_fragment = GetCoreFragment(refmol, target_molecules)
    oechem.OESuppressHydrogens(core_fragment)
    print(f'  Core fragment has {core_fragment.NumAtoms()} heavy atoms')

    # Write core fragment (without modifying it)
    with oechem.oemolostream(f'{prefix}-core-{fragment}.mol2') as ofs:
        oechem.OEWriteMolecule(ofs, oechem.OEGraphMol(core_fragment))

    # Expand conformers
    with oechem.oemolostream(prefix + f'-conformers-{fragment}.sdf') as ofs:
        # Write reference molecule copy
        refmol_copy = oechem.OEGraphMol(refmol)
        oechem.OESetSDData(refmol_copy, 'clash_score', '0.0')
        oechem.OEWriteMolecule(ofs, refmol_copy)

        #from rich.progress import track
        #for mol in track(target_molecules, f'Generating poses for {len(target_molecules)} target molecules'):
        from tqdm import tqdm
        for mol in tqdm(target_molecules):
            pose = generate_restricted_conformers(receptor, core_fragment, mol)
            if pose is not None:
                oechem.OEWriteMolecule(ofs, pose)

if __name__ == '__main__':
    #fragment = 'x2646' # TRY-UNI-714a760b-6 (the main amino pyridine core)
    #fragment = 'x10789' # TRY-UNI-2eddb1ff-7 (beta-lactam an chloride)

    # TODO: Figure out why this single-threading workaround is needed to prevent issues
    from openeye import oechem
    #oechem.OESetMemPoolMode(oechem.OEMemPoolMode_SingleThreaded |
    #                        oechem.OEMemPoolMode_UnboundedCache)

    for fragment in ['x10789']:
        for prefix in [
                'nucleophilic_displacement_enumeration_for_FEP-permuted',
                #'activity-data-2020-07-29',
                #'primary_amine_enumeration_for_chodera_lab_FEP-permuted',
                #'boronic_ester_enumeration_for_chodera_lab_FEP-permuted',
        ]:
            #generate_poses(fragment, prefix, fragment_title='TRY-UNI-2eddb1ff-7', filter_series="3-aminopyridine-strict")
            generate_poses(fragment, prefix, fragment_title='TRY-UNI-2eddb1ff-7')

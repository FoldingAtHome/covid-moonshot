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

    print("Number of molecules = %d" % len(mols))

    frags = GetFragments(refmol, minbonds, maxbonds)
    if len(frags) == 0:
        oechem.OEThrow.Error("No fragment is enumerated with bonds %d-%d!" % (minbonds, maxbonds))

    print("Number of fragments = %d" % len(frags))

    commonfrags = GetCommonFragments(mols, frags, atomexpr, bondexpr)
    if len(commonfrags) == 0:
        oechem.OEThrow.Error("No common fragment is found!")

    print("Number of common fragments = %d" % len(commonfrags))

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
    # Generate a dense sampling of conformers
    omegaOpts = oeomega.OEOmegaOptions()
    #omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Dense)
    # Set the fixed reference molecule
    omegaFixOpts = oeomega.OEConfFixOptions()
    omegaFixOpts.SetFixMaxMatch(10) # allow multiple MCSS matches
    omegaFixOpts.SetFixDeleteH(True) # only use heavy atoms

    #atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_Aromaticity | oechem.OEExprOpts_RingMember
    #bondexpr = oechem.OEExprOpts_BondOrder | oechem.OEExprOpts_RingMember
    #atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_RingMember
    #bondexpr = oechem.OEExprOpts_BondOrder | oechem.OEExprOpts_RingMember
    #omegaFixOpts.SetAtomExpr(atomexpr)
    #omegaFixOpts.SetBondExpr(bondexpr)
    omegaFixOpts.SetFixMol(fixmol)
    omegaOpts.SetConfFixOptions(omegaFixOpts)
    omegaOpts.SetWarts(False) # expand molecule title
    omegaOpts.SetStrictStereo(True) # set strict stereochemistry
    omegaOpts.SetIncludeInput(False) # don't include input
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

    # Select the conformer with the fewest steric clashes
    bump_check = BumpCheck(receptor)
    poses = list()
    for pose in mol.GetConfs():
        poses.append( (pose, bump_check.count(pose)) )
    poses = sorted(poses, key=lambda x : x[1])
    mol.SetActive(poses[0][0])
    #print([pose[1] for pose in poses])

    #from openeye import oedocking
    #score = oedocking.OEScore(oedocking.OEScoreType_Chemgauss4)
    #score.Initialize(receptor)
    #scores = list()
    #for pose in mol.GetConfs():
    #    scores.append(score.ScoreLigand(pose))
    #print(scores)

    # Convert to single-conformer molecule
    mol = oechem.OEGraphMol(mol)

    return mol

if __name__ == '__main__':
    fragment = 'x2646' # TRY-UNI-714a760b-6 (the main amino pyridine core)
    fragment = 'x10789' # TRY-UNI-2eddb1ff-7 (beta-lactam an chloride)

    # Read receptor
    print('Reading receptor...')
    from openeye import oechem
    receptor = oechem.OEGraphMol()
    with oechem.oemolistream(f'receptors/monomer/Mpro-{fragment}_0_bound-receptor.oeb.gz') as infile:
        oechem.OEReadMolecule(infile, receptor)
        print(receptor.NumAtoms())

    # Read target molecules
    print('Reading target molecules...')
    from openeye import oechem
    prefix = 'primary_amine_enumeration_for_chodera_lab_FEP-permuted'
    prefix = 'boronic_ester_enumeration_for_chodera_lab_FEP-permuted'
    targets_filename = prefix + '.csv'
    target_molecules = list()
    with oechem.oemolistream(targets_filename) as ifs:
        for mol in ifs.GetOEGraphMols():
            target_molecules.append( oechem.OEGraphMol(mol) )
    print(f'There are {len(target_molecules)} target molecules')

    # Expand uncertain stereochemistry
    print('Expanding uncertain stereochemistry...')
    target_molecules = expand_stereochemistry(target_molecules)
    print(f'There are {len(target_molecules)} target molecules')

    # Read reference molecule with coordinates
    refmol_filename = f'receptors/monomer/Mpro-{fragment}_0_bound-ligand.mol2'
    with oechem.oemolistream(refmol_filename) as ifs:
        for mol in ifs.GetOEGraphMols():
            refmol = mol
            break
    print(f'Reference molecule as {refmol.NumAtoms()} atoms')

    # Get core fragment
    print('Identifying core fragment...')
    core_fragment = GetCoreFragment(refmol, target_molecules)
    oechem.OESuppressHydrogens(core_fragment)
    print(f'Core fragment has {core_fragment.NumAtoms()} atoms')

    # Write core fragment (without modifying it)
    with oechem.oemolostream(f'core-fragment-{fragment}.mol2') as ofs:
        oechem.OEWriteMolecule(ofs, oechem.OEGraphMol(core_fragment))

    # Expand conformers
    with oechem.oemolostream(prefix + f'-conformers-{fragment}.sdf') as ofs:
        # Write reference molecule copy
        oechem.OEWriteMolecule(ofs, oechem.OEGraphMol(refmol))

        from rich.progress import track
        for mol in track(target_molecules, 'Expanding conformers'):
            pose = generate_restricted_conformers(receptor, core_fragment, mol)
            if pose is not None:
                oechem.OEWriteMolecule(ofs, pose)

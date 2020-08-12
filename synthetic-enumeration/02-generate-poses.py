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
    #for frag in track(frags, description='Finding common fragments'):
    for frag in frags:
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

def generate_restricted_conformers(receptor, refmol, mol):
    """
    Generate and select a conformer of the specified molecule using the reference molecule

    Parameters
    ----------
    receptor : openeye.oechem.OEGraphMol
        Receptor (already prepped for docking) for identifying optimal pose
    refmol : openeye.oechem.OEGraphMol
        Reference molecule which shares some part in common with the proposed molecule
    mol : openeye.oechem.OEGraphMol
       The molecule to enumerate conformers for
    """
    from openeye import oechem, oeomega

    # Get core fragment
    #print('Identifying core fragment...')
    atomexpr = oechem.OEExprOpts_Hybridization
    bondexpr = oechem.OEExprOpts_Aromaticity

    core_fragment = GetCoreFragment(refmol, [mol], atomexpr=atomexpr, bondexpr=bondexpr)
    oechem.OESuppressHydrogens(core_fragment)
    fixmol = core_fragment
    #print(f'  Core fragment has {core_fragment.NumAtoms()} heavy atoms')

    MIN_CORE_ATOMS = 6
    if core_fragment.NumAtoms() < MIN_CORE_ATOMS:
        return None

    # Write core fragment (without modifying it)
    #with oechem.oemolostream(f'{prefix}-core-{fragment}.mol2') as ofs:
    #    oechem.OEWriteMolecule(ofs, oechem.OEGraphMol(core_fragment))

    # Create an Omega instance
    #omegaOpts = oeomega.OEOmegaOptions()
    omegaOpts = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Dense)

    # DEBUG
    atomexpr = oechem.OEExprOpts_DefaultAtoms
    bondexpr = oechem.OEExprOpts_DefaultBonds
    atomexpr = oechem.OEExprOpts_Hybridization
    bondexpr = oechem.OEExprOpts_Aromaticity

    # Set the fixed reference molecule
    omegaFixOpts = oeomega.OEConfFixOptions()
    omegaFixOpts.SetFixMaxMatch(10) # allow multiple MCSS matches
    omegaFixOpts.SetFixDeleteH(True) # only use heavy atoms
    omegaFixOpts.SetFixMol(fixmol)
    omegaFixOpts.SetAtomExpr(atomexpr) # DEBUG
    omegaFixOpts.SetBondExpr(bondexpr) # DEBUG
    omegaOpts.SetConfFixOptions(omegaFixOpts)

    molBuilderOpts = oeomega.OEMolBuilderOptions()
    molBuilderOpts.SetStrictAtomTypes(False) # don't give up if MMFF types are not found
    omegaOpts.SetMolBuilderOptions(molBuilderOpts)

    omegaOpts.SetWarts(False) # expand molecule title
    omegaOpts.SetStrictStereo(False) # set strict stereochemistry
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
def get_series(mol):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import Descriptors
    series_SMARTS_dict = {
        #"3-aminopyridine": "[R1][C,N;R0;!$(NC(=O)CN)]C(=O)[C,N;R0;!$(NC(=O)CN)][c]1cnccc1",
        "3-aminopyridine-like": "[R1]!@[C,N]C(=O)[C,N]!@[R1]",
        "3-aminopyridine-strict": "c1ccncc1NC(=O)!@[R1]",
        "Ugi": "[c,C:1][C](=[O])[N]([c,C,#1:2])[C]([c,C,#1:3])([c,C,#1:4])[C](=[O])[NH1][c,C:5]",
        "quinolones": "NC(=O)c1cc(=O)[nH]c2ccccc12",
        "piperazine-chloroacetamide": "O=C(CCl)N1CCNCC1",
    }

    smi = oechem.OECreateSmiString(mol)

    # Filter out covalent
    try:
        if oechem.OEGetSDData(mol,'acrylamide')=='True' or oechem.OEGetSDData(mol,'chloroacetamide')=='True':
            return None
    except Exception as e:
        print(e)

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

def generate_restricted_conformers_star(args):
    return generate_restricted_conformers(*args)

def generate_poses(receptor, refmol, target_molecules, output_filename):
    """
    Parameters
    ----------
    receptor : openeye.oechem.OEGraphMol
        Receptor (already prepped for docking) for identifying optimal pose
    refmol : openeye.oechem.OEGraphMol
        Reference molecule which shares some part in common with the proposed molecule
    target_molecules : list of OEMol
        List of molecules to build
    output_filename : str
        Output filename for generated conformers
    """
    # Expand uncertain stereochemistry
    print('Expanding uncertain stereochemistry...')
    target_molecules = expand_stereochemistry(target_molecules)
    print(f'  There are {len(target_molecules)} target molecules')

    # Identify optimal conformer for each molecule
    with oechem.oemolostream(output_filename) as ofs:
        # Write reference molecule copy
        refmol_copy = oechem.OEGraphMol(refmol)
        oechem.OESetSDData(refmol_copy, 'clash_score', '0.0')
        oechem.OEWriteMolecule(ofs, refmol_copy)

        from rich.progress import track
        #for mol in track(target_molecules, f'Generating poses for {len(target_molecules)} target molecules'):
        from multiprocessing import Pool
        from tqdm import tqdm

        pool = Pool()
        args = [ (receptor, refmol, mol) for mol in target_molecules ]
        for pose in track(pool.imap_unordered(generate_restricted_conformers_star, args), total=len(args), description='Enumerating conformers...'):
            if pose is not None:
                oechem.OEWriteMolecule(ofs, pose)
        pool.close()
        pool.join()

        #for mol in tqdm(target_molecules):
        #    pose = generate_restricted_conformers(receptor, core_fragment, mol)
        #    if pose is not None:
        #        oechem.OEWriteMolecule(ofs, pose)

if __name__ == '__main__':
    #fragment = 'x2646' # TRY-UNI-714a760b-6 (the main amino pyridine core)
    #fragment = 'x10789' # TRY-UNI-2eddb1ff-7 (beta-lactam an chloride)

    # TODO: Figure out why this single-threading workaround is needed to prevent issues
    from openeye import oechem
    #oechem.OESetMemPoolMode(oechem.OEMemPoolMode_SingleThreaded |
    #                        oechem.OEMemPoolMode_UnboundedCache)

    assay_data_filename = 'activity-data-2020-07-29.csv'
    fragments = {
        'x10789' : 'TRY-UNI-2eddb1ff-7',
        }

    # Load assay data if available
    assayed_molecules = list()
    with oechem.oemolistream(assay_data_filename) as ifs:
        for mol in ifs.GetOEGraphMols():
            assayed_molecules.append( oechem.OEGraphMol(mol) )
    print(f'  There are {len(assayed_molecules)} assayed molecules')

    # Load all fragments
    for prefix in [
                'activity-data-2020-08-11',
                #'aminopyridine-retrospective-jdc-2020-08-11',
                #'fastgrant-table1',
                #'aminopyridine_compounds_for_FEP_benchmarking',
                #'nucleophilic_displacement_enumeration_for_FEP-permuted',
                #'activity-data-2020-07-29',
                #'primary_amine_enumeration_for_chodera_lab_FEP-permuted',
                #'boronic_ester_enumeration_for_chodera_lab_FEP-permuted',
        ]:
        for fragment in fragments:

            # Read receptor
            print('Reading receptor...')
            from openeye import oechem
            receptor = oechem.OEGraphMol()
            receptor_filename = f'../receptors/monomer/Mpro-{fragment}_0_bound-receptor.oeb.gz'
            from openeye import oedocking
            oedocking.OEReadReceptorFile(receptor, receptor_filename)
            print(f'  Receptor has {receptor.NumAtoms()} atoms.')

            # Read reference fragment with coordinates
            refmol_filename = f'../receptors/monomer/Mpro-{fragment}_0_bound-ligand.mol2'
            refmol = None
            with oechem.oemolistream(refmol_filename) as ifs:
                for mol in ifs.GetOEGraphMols():
                    refmol = mol
                    break
            if refmol is None:
                raise Exception(f'Could not read {refmol_filename}')
            print(f'Reference molecule has {refmol.NumAtoms()} atoms')
            # Replace title
            refmol.SetTitle(fragments[fragment])
            # Copy data from assayed molecules (if present)
            for mol in assayed_molecules:
                if mol.GetTitle() == fragments[fragment]:
                    print(f'{refmol.GetTitle()} found in target_molecules; copying SDData')
                    oechem.OECopySDData(refmol, mol)
                    break

            # Read target molecules
            target_molecules_filename = prefix + f'.csv'
            print('Reading target molecules...')
            from openeye import oechem
            target_molecules = list()
            with oechem.oemolistream(target_molecules_filename) as ifs:
                for mol in ifs.GetOEGraphMols():
                    target_molecules.append( oechem.OEGraphMol(mol) )
            if len(target_molecules) == 0:
                raise Exception('No target molecules specified; check filename!')
            print(f'  There are {len(target_molecules)} target molecules')

            # Filter series and include only those that include the required scaffold
            filter_series = '3-aminopyridine-like'
            #filter_series = None
            if filter_series is not None:
                print(f'Filtering out series {filter_series}...')
                target_molecules = [ mol for mol in target_molecules if get_series(mol) == filter_series ]
                print(f'  There are {len(target_molecules)} target molecules')
                with oechem.oemolostream(f'filtered.mol2') as ofs:
                    for mol in target_molecules:
                        oechem.OEWriteMolecule(ofs, oechem.OEGraphMol(mol))

            # Filter series to include only those with IC50s
            filter_IC50 = True
            if filter_IC50:
                print(f'Retaining only molecules with IC50s...')
                target_molecules = [ mol for mol in target_molecules if len(oechem.OEGetSDData(mol, 'f_avg_pIC50'))>0 ]
                print(f'  There are {len(target_molecules)} target molecules')


            # Generate poses for all molecules
            output_filename = f'{prefix}-dockscores-{fragment}.sdf'
            generate_poses(receptor, refmol, target_molecules, output_filename)

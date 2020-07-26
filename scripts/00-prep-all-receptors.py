"""
Prepare all SARS-CoV-2 Mpro structures for docking and simulation in monomer and dimer forms

This should be run from the covid-moonshot/scripts directory

"""

# If structural data is not present, download and unpack it
import os

structures_path = '../structures'
output_basepath = '../receptors'

def download_url(url, save_path, chunk_size=128):
    import os
    base_path, filename = os.split(save_path)
    if not os.path.exists(base_path):
        os.makedirs(base_path)

    import requests
    r = requests.get(url, stream=True)
    with open(save_path, 'wb') as fd:
        for chunk in r.iter_content(chunk_size=chunk_size):
            fd.write(chunk)

def read_pdb_file(pdb_file):
    print(f'Reading receptor from {pdb_file}...')

    from openeye import oechem
    ifs = oechem.oemolistream()
    ifs.SetFlavor(oechem.OEFormat_PDB, oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)  # noqa

    if not ifs.open(pdb_file):
        oechem.OEThrow.Fatal("Unable to open %s for reading." % pdb_file)

    mol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ifs, mol):
        oechem.OEThrow.Fatal("Unable to read molecule from %s." % pdb_file)
    ifs.close()

    return (mol)

def prepare_receptor(complex_pdb_filename, output_basepath, dimer=False):
    """
    Parameters
    ----------
    complex_pdb_filename : str
        The complex PDB file to read in
    output_basepath : str
        Base path for output
    dimer : bool, optional, default=False
        If True, generate the dimer as the biological unit
    """
    import os
    basepath, filename = os.path.split(complex_pdb_filename)
    prefix, extension = os.path.splitext(filename)
    prefix = os.path.join(output_basepath, prefix)

    # Check if receptor already exists
    receptor_filename = f'{prefix}-receptor.oeb.gz'
    thiolate_receptor_filename = f'{prefix}-receptor-thiolate.oeb.gz'
    if os.path.exists(receptor_filename) and os.path.exists(thiolate_receptor_filename):
        return

    # Read in PDB file
    pdbfile_lines = [ line for line in open(complex_pdb_filename, 'r') if 'UNK' not in line ]

    # If monomer is specified, drop crystal symmetry lines
    if not dimer:
        pdbfile_lines = [ line for line in pdbfile_lines if 'REMARK 350' not in line ]

    # Reconstruct PDBFile contents
    pdbfile_contents = ''.join(pdbfile_lines)

    # Read the receptor and identify design units
    from openeye import oespruce, oechem
    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile(delete=False, mode='wt', suffix='.pdb') as pdbfile:
        pdbfile.write(pdbfile_contents)
        pdbfile.close()
        complex = read_pdb_file(pdbfile.name)
        # TODO: Clean up

    print('Identifying design units...')
    design_units = list(oespruce.OEMakeDesignUnits(complex))
    if len(design_units) == 1:
        design_unit = design_units[0]
    elif len(design_units) > 1:
        print('More than one design unit found---using first one')
        design_unit = design_units[0]
    elif len(design_units) == 0:
        raise Exception('No design units found')

    # Prepare the receptor
    print('Preparing receptor...')
    from openeye import oedocking
    protein = oechem.OEGraphMol()
    design_unit.GetProtein(protein)
    ligand = oechem.OEGraphMol()
    design_unit.GetLigand(ligand)

    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, protein, ligand)
    oedocking.OEWriteReceptorFile(receptor, receptor_filename)

    with oechem.oemolostream(f'{prefix}-protein.pdb') as ofs:
        oechem.OEWriteMolecule(ofs, protein)
    with oechem.oemolostream(f'{prefix}-ligand.mol2') as ofs:
        oechem.OEWriteMolecule(ofs, ligand)
    with oechem.oemolostream(f'{prefix}-ligand.pdb') as ofs:
        oechem.OEWriteMolecule(ofs, ligand)
    with oechem.oemolostream(f'{prefix}-ligand.sdf') as ofs:
        oechem.OEWriteMolecule(ofs, ligand)

    # Filter out UNK from PDB files (which have covalent adducts)
    pdbfile_lines = [ line for line in open(f'{prefix}-protein.pdb', 'r') if 'UNK' not in line ]
    with open(f'{prefix}-protein.pdb', 'wt') as outfile:
        outfile.write(''.join(pdbfile_lines))

    # Adjust protonation state of CYS145 to generate thiolate form
    print('Deprotonating CYS145...')
    pred = oechem.OEAtomMatchResidue(["CYS:145: :A"])
    for atom in protein.GetAtoms(pred):
        if oechem.OEGetPDBAtomIndex(atom) == oechem.OEPDBAtomName_SG:
            oechem.OESuppressHydrogens(atom)
            atom.SetFormalCharge(-1)
            atom.SetImplicitHCount(0)
    # Adjust protonation states
    print('Re-optimizing hydrogen positions...')
    place_hydrogens_opts = oechem.OEPlaceHydrogensOptions()
    place_hydrogens_opts.SetBypassPredicate(pred)
    protonate_opts = oespruce.OEProtonateDesignUnitOptions(place_hydrogens_opts)
    success = oespruce.OEProtonateDesignUnit(design_unit, protonate_opts)
    design_unit.GetProtein(protein)

    # Old hacky way to adjust protonation states
    #opts = oechem.OEPlaceHydrogensOptions()
    #opts.SetBypassPredicate(pred)
    #describe = oechem.OEPlaceHydrogensDetails()
    #success = oechem.OEPlaceHydrogens(protein, describe, opts)
    #if success:
    #    oechem.OEUpdateDesignUnit(design_unit, protein, oechem.OEDesignUnitComponents_Protein)

    # Write thiolate form of receptor
    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, protein, ligand)
    oedocking.OEWriteReceptorFile(receptor, thiolate_receptor_filename)

    with oechem.oemolostream(f'{prefix}-protein-thiolate.pdb') as ofs:
        oechem.OEWriteMolecule(ofs, protein)

    # Filter out UNK from PDB files (which have covalent adducts)
    pdbfile_lines = [ line for line in open(f'{prefix}-protein-thiolate.pdb', 'r') if 'UNK' not in line ]
    with open(f'{prefix}-protein-thiolate.pdb', 'wt') as outfile:
        outfile.write(''.join(pdbfile_lines))


if __name__ == '__main__':
    # Prep all receptors
    import glob, os

    if not os.path.exists(structures_path):
        # Download ZIP file
        url = 'https://fragalysis.diamond.ac.uk/media/targets/Mpro.zip'
        zip_path = os.path.join(structures_path, 'Mpro.zip')
        download_url(url, zip_path)
        # Unpack ZIP file
        from zipfile import ZipFile
        with ZipFile(, 'r') as zip_obj:
           zip_obj.extractall(structures_path)

    # Get list of all PDB files to prep
    source_pdb_files = glob.glob(os.path.join(structures_path, "aligned/Mpro-*_0/Mpro-*_0_bound.pdb"))

    # Create output directory
    os.makedirs(output_basepath, exist_ok=True)

    for dimer in [False, True]:
        if dimer:
            output_basepath = 'receptors/dimer'
        else:
            output_basepath = 'receptors/monomer'

        os.makedirs(output_basepath, exist_ok=True)

        def prepare_receptor_wrapper(complex_pdb_file):
            try:
                prepare_receptor(complex_pdb_file, output_basepath, dimer=dimer)
            except Exception as e:
                print(e)

        # Process all receptors in parallel
        from multiprocessing import Pool
        from tqdm import tqdm
        with Pool() as pool:
            max_ = len(source_pdb_files)
            with tqdm(total=max_) as pbar:
                for i, _ in enumerate(pool.imap_unordered(prepare_receptor_wrapper, source_pdb_files)):
                    pbar.update()

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

def upload_receptors(all_fragments, project=2509):
    receptor_datasets = dict()
    for fragment in all_fragments:
        dataset_name = f'Mpro-{fragment}-receptor'
        receptor_filename = os.path.join(f'../receptors/Mpro-{fragment}-receptor.oeb.gz')
        print(receptor_filename)

        # Check if exists
        filters = {"name": f"Mpro-{fragment}-receptor", "project": project}
        existing_datasets = [dataset for dataset in APISession.list_resources(Dataset, filters=filters)]
        if len(existing_datasets) == 0:
            print('  Uploading...')
            id = Dataset.upload(APISession, dataset_name, receptor_filename, project=project)
            receptor_datasets[fragment] = id
        else:
            print('  Already exists, skipping.')
            receptor_datasets[fragment] = existing_datasets[0].id

    print(receptor_datasets)
    import pickle
    with open('receptor-datasets.pkl', 'wb') as handle:
        pickle.dump(receptor_datasets, handle, protocol=pickle.HIGHEST_PROTOCOL)

def dock_to_receptors(all_fragments, conformers_fileid=16535, project=2509):
    for fragment in all_fragments:
        filters = {"name": f"Mpro-{fragment}-receptor", "project": project}
        receptor_datasets = [dataset for dataset in APISession.list_resources(Dataset, filters=filters)]
        receptor_fileid = receptor_datasets[0].id
        cmd = f'ocli --json jobs start 10634 "docking to {fragment}" --in={conformers_fileid} --receptor={receptor_fileid} --out "covid_submissions_03_26_2020 - docked to {fragment}"'
        print(cmd)

        import subprocess
        returned_output = subprocess.check_output(cmd, shell=True)
        print(returned_output)

if __name__ == '__main__':
    # Dock the ligands
    import os
    from openeye import oechem

    # Generate list of all X-ray fragments
    fragment_molecules = read_csv_molecules(os.path.join('../molecules', 'mpro_fragments_03_25_2020.csv'))
    all_fragments = [ oechem.OEGetSDData(molecule, "fragments") for molecule in fragment_molecules ]

    # Filter fragments to only retain those with prepped receptor structures
    all_fragments = [ fragment for fragment in all_fragments if os.path.exists(os.path.join(f'../receptors/Mpro-{fragment}-receptor.oeb.gz')) ]

    from orionclient.session import APISession
    from orionclient.types import Dataset

    # Upload receptor datasets
    upload_receptors(all_fragments)

    # Get all receptor datasets
    dock_to_receptors(all_fragments)

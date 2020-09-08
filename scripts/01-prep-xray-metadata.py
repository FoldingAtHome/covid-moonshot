"""
Compile metadata on X-ray structures for FAH

This script ingests `../structures/metadata.csv`:
```
,crystal_name,RealCrystalName,smiles,new_smiles,alternate_name,site_name,pdb_entry
1,Mpro-1q2w-2020-04-Bonanno_0,Mpro-1q2w-2020-04-Bonanno,C[C@H](O)CC(C)(C)O,NA,NA,Mpro-SARS1,1Q2W
3,Mpro-1wof-2020-04-Yang_0,Mpro-1wof-2020-04-Yang,CCOC(O)CC[C@H](C[C@@H]1CCNC1O)N[C@H](O)[C@H](CC(C)C)NC(O)[C@@H](NC(O)[C@H](C)NC(O)C1CC(C)ON1)C(C)C,NA,NA,Mpro-SARS1,1WOF
4,Mpro-2a5i-2020-04-Lee_0,Mpro-2a5i-2020-04-Lee,CCOC(O)[C@@H](O)C[C@@H](O)N(CCC(N)O)NC(O)[C@H](CC1CCCCC1)N[C@H](O)[C@H](CC(C)C)N[C@H](O)OCC1CCCCC1,NA,NA,Mpro-SARS1,2A5I
...
```
and adds new unique systems to `../xray-fah/metadata.csv`
```
RUN,crystal_name,RealCrystalName,smiles,new_smiles,alternate_name,site_name,pdb_entry
RUN0,Mpro-1q2w-2020-04-Bonanno_0,Mpro-1q2w-2020-04-Bonanno,C[C@H](O)CC(C)(C)O,NA,NA,Mpro-SARS1,1Q2W
RUN1,Mpro-1wof-2020-04-Yang_0,Mpro-1wof-2020-04-Yang,CCOC(O)CC[C@H](C[C@@H]1CCNC1O)N[C@H](O)[C@H](CC(C)C)NC(O)[C@@H](NC(O)[C@H](C)NC(O)C1CC(C)ON1)C(C)C,NA,NA,Mpro-SARS1,1WOF
RUN2,Mpro-2a5i-2020-04-Lee_0,Mpro-2a5i-2020-04-Lee,CCOC(O)[C@@H](O)C[C@@H](O)N(CCC(N)O)NC(O)[C@H](CC1CCCCC1)N[C@H](O)[C@H](CC(C)C)N[C@H](O)OCC1CCCCC1,NA,NA,Mpro-SARS1,2A5I
...
```
This MUST be run before `01-prep-xray-for-fah.py`.

"""

if __name__ == '__main__':
    # Parse arguments
    import argparse

    parser = argparse.ArgumentParser(description='Ingest new metadata')
    parser.add_argument('--structures', dest='structures_path', type=str, default='../structures',
                        help='directory of receptor conformations (default: ../structures)')
    parser.add_argument('--metadata', dest='metadata_filename', type=str, default='../fah-xray/fah-metadata.csv',
                        help='FAH metadata filename (default: ../fah-xray/fah-metadata.csv))')

    args = parser.parse_args()

    # Open FAH metadata file
    from collections import OrderedDict
    metadata = OrderedDict()
    crystal_names = set()
    import os, csv
    if os.path.exists(args.metadata_filename):
        print(f'Reading known structures from {args.metadata_filename}')
        with open(args.metadata_filename, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                run = row['RUN']
                metadata[run] = row
                crystal_names.add(row['crystal_name'])
        print(f'{len(metadata)} entries in metadata file')

    # Read DiamondMX/XChem structure medatadata
    print(f'Ingesting new structures from {args.structures_path}')
    structure_metadata_filename = os.path.join(args.structures_path, 'metadata.csv')
    with open(structure_metadata_filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['crystal_name'] not in crystal_names:
                # Add new entry
                run = f'RUN{len(metadata)}'
                row.pop('')
                row['RUN'] = run
                row.move_to_end('RUN', last=False)
                metadata[run] = row
        print(f'{len(metadata)} entries in metadata file')

    # Write metadata
    print(f'Writing out new {args.metadata_filename}')
    with open(args.metadata_filename, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=metadata['RUN0'].keys())
        writer.writeheader()
        for run in metadata.keys():
            writer.writerow(metadata[run])

"""
Rename compounds in virtual libraries and reformat CSV
"""

libraries = [
    'amidation',
    'mitsunobu',
    'red_ami',
    'red_amination_1',
    'williamson_ether'
    ]

for library in libraries:
    import csv
    csv_filename = f'{library}.csv'
    out_filename = f'{library}-renamed.csv'
    with open(csv_filename, 'r') as csvfile, open(out_filename, 'wt') as outfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        next(csvreader, None)  # skip the headers
        outfile.write(f'SMILES,TITLE\n') # write header
        for row in csvreader:
            title = row[0]
            smiles = row[1]
            title = f'{library}-{title}'
            # Write new entry
            outfile.write(f'{smiles},{title}\n')

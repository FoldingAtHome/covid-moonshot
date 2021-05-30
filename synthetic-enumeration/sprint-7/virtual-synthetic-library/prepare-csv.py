# Add naming to the file

def convert(input_filename, output_filename, name_prefix):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        index = 0
        import csv
        csvreader = csv.reader(infile, delimiter=',')
        for row in csvreader:
            smiles = row[0]
            if smiles == 'product_SMILES': # header
                outfile.write(f'SMILES,title\n')            
            else:
                outfile.write(f'{smiles},{name_prefix}-{index:07d}\n')
                index += 1

convert('source/reductive_aminations_for_computational_triage.csv', 'reductive-aminations.csv', 'amination')

            

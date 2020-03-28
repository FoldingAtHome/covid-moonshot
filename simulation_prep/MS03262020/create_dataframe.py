import os, glob, tqdm
import pandas as pd

infiles = ['covid_submissions_03_26_2020 - docked.sdf']
columns = ['fragments','OEChem','Hybrid2']
data = [[0]*len(columns)]

for file in infiles:
    with open(file) as f:
        lines = f.readlines()
    lines = [line.split('\n')[0] for line in lines]
    ligand_ndx = 0
    for j,line in tqdm.tqdm(enumerate(lines)):
        if '$$$$' in line:
            ligand_ndx += 1
            data.append([0]*len(columns))
            continue
        for field in columns:
            if field in line:
                if field == 'OEChem':
                    data[ligand_ndx][columns.index(field)] = lines[j-1]
                    break
                else:
                    data[ligand_ndx][columns.index(field)] = lines[j+1]
                    break
    df = pd.DataFrame(data, columns=columns)
    print(df)
    df.to_pickle('ligands.pkl')

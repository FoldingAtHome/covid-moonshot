#!/usr/bin/env python

import os, shutil, tqdm

runs_per_project = 500
project = 14364 # start RL from this project number
project_type = 'gromacs' # gromacs or openmm
system_type = 'L' # RL or L
dataset = 'MS03262020'
debug = False
total_projects = 1985


start = range(1,total_projects,runs_per_project)
projects = [project + x for x in range(len(start))]

print(start)
print(projects)

lig = 1 # ligands are 1-indexed

for project_ndx, start in tqdm.tqdm(enumerate(start)):
    project = projects[project_ndx]
    if not os.path.exists(f'projects/p{project}'):
        os.makedirs(f'projects/p{project}')

    run = 0 # run numbers are 0-indexed

    while run < runs_per_project * (project_ndx+1):

        if lig == total_projects: # if we get to the end, stop
            sys.exit()

        if not os.path.exists(f'{dataset}_{system_type}/LIG{lig}/conf.gro'): # if no topology generated, skip
            if debug:
                print(f'*** Found bad/incomplete ligand: {dataset}_{system_type}/LIG{lig}')
            lig += 1
            continue

        # skip if ligand was already processed
        if os.path.exists(f'projects/p{project}/RUN{run}/conf.gro'):
            lig += 1
            run += 1
            continue
    
        # otherwise copy the data over to the output directory
        if debug:
            print(f'copying {dataset}_{system_type}/LIG{lig} to projects/p{project}/RUN{run}')
        else:
            os.makedirs(f'projects/p{project}/RUN{run}')
            if project_type == 'gromacs':
                shutil.copy2(f'{dataset}_{system_type}/LIG{lig}/conf.gro',f'projects/p{project}/RUN{run}/conf.gro')
                shutil.copy2(f'{dataset}_{system_type}/LIG{lig}/topol.top',f'projects/p{project}/RUN{run}/topol.top')
            with open(f'projects/p{project}/RUN{run}/note','w') as f:
                f.write(f'{dataset}\n{system_type}\nLIG{lig}')
        
        # prepare for next iteration    
        run += 1
        lig += 1        

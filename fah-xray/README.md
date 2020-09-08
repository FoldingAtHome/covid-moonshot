# Prepare all Mpro X-ray structures from Fragalysis for simulation on Folding@home in multiple states

## Process

First, run `../scripts/01-prep-xray-metadata.py` to prepare metadata.

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

Next, prepare X-ray structures in parallel with `bsub < run-lsf.sh` after adjusting job ID span.

This script prepares the folowing projects:
``` 
13430 : apo Mpro monomer His41(0) Cys145(0)
13431 : apo Mpro monomer His41(+) Cys145(-)
13432 : holo Mpro monomer His41(0) Cys145(0)
13433 : holo Mpro monomer His41(+) Cys145(-)
13434 : apo Mpro dimer His41(0) Cys145(0)
13435 : apo Mpro dimer His41(+) Cys145(-)
13436 : holo Mpro dimer His41(0) Cys145(0)
13437 : holo Mpro dimer His41(+) Cys145(-)
```
Each RUN corresponds to a different fragment structure

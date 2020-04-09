# Docking and free energy calculations for the COVID Moonshot

This repository contains all the scripts and output for docking and free energy calculations of [COVID Moonshot](https://covid.postera.ai/covid) submissions run on [Folding@home](http://foldingathome.org).

## Author

John D. Chodera (MSKCC) `<john.chodera@choderalab.org>`

## Manifest
* `diamond-structures/` - source structures of SARS-CoV-2 main viral protease from [DiamondMX/XChem](https://www.diamond.ac.uk/covid-19/for-scientists/Main-protease-structure-and-XChem.html)
* `receptors/` - receptor structures prepared for docking and simulation by the [OpenEye Spruce Toolkit](https://docs.eyesopen.com/toolkits/python/sprucetk/index.html)
* `scripts/` - scripts for preparing receptors, docking ligands, and setting up free energy calculations
* `environment.yml` - [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) used for calculations
* `moonshot-submissions/` - input files and docking results for submitted batches of [COVID Moonshot](https://covid.postera.ai/covid) compounds
* `covalent-docking/` - input files and constrained docking results for subset of [COVID Moonshot](https://covid.postera.ai/covid) compounds with covalent warheads

## Outdated
* `redock-fragments/` - generation of ROC to assess recovery of hits from missess in screened fragments
* `2020-04-01-in-stock-and-one-step-synthesis/` - docking of high-priority compounds for initial triage round
* `melatonin/` - docking of melatonin and metabolites after observation that x0104 is extremely similar to melatonin
* `simulation_prep/` - historical scripts and input files used in preparing gromacs expanded ensemble free energy calculations for Folding@home

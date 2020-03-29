# Free energy calculations for the COVID Moonshot

John D. Chodera <john.chodera@choderalab.org>

## Manifest
* `diamond-structures/` - structures of SARS-CoV-2 main viral protease from [DiamondMX/XChem](https://www.diamond.ac.uk/covid-19/for-scientists/Main-protease-structure-and-XChem.html)
* `receptors/` - prepared receptor structures (in OpenEye `.oeb.gz` format as well as PDB) with prepared fragments
* `molecules/` - input molecule sets for docking
* `scripts/` - scripts for preparing receptors and docking and scoring ligands
* `docking/` - docked structures (and scores) for use in free energy calculations
* `simulation_prep/` - scripts used in preparing gromacs expanded ensemble free energy calculations for Folding@home

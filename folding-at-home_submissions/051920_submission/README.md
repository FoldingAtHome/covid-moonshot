#  Folding@Home absolute FEP from ligands docked to fragment-bound Mpro structures
## Release of ranked compounds 2020-05-12

### Files found in this release:
* COVID Moonshot 051220 Release.pdf - The announcement of this submitted batch of results.
* Plots/ - Contains all images found within the announcement.
* compound-set_FoldingAtHome_EE-FEP.sdf - SDF file containing all of the ligands submitted with tags for binding_free_energy, SMILES, associated receptor file, etc.
* master_results_WL0.12.pkl - A pickled dataframe containing more detailed results for each of the submitted ligands
* ../receptors.zip - A compressed directory containing all relevant protein structures.

##Overview

This submission contains results for the top 286 molecules ranked using absolute binding free energy perturbation (FEP) simulations performed on the massively parallel Folding@home distributed computing platform. A full synopsis of the methods and analysis used can be found in our preliminary submission post:
https://discuss.postera.ai/t/foldingathome-absolute-fep-from-ligands-docked-to-fragment-bound-mpro-structures/1352

The final score is a free energy quantity: the more negative, the better. A linear function can convert this quantity to a standard-concentration binding free energy.

Note that these 286 compounds represent a snapshot (from 05-19-2020) of data that meet our convergence criteria. These simulations will continue to run on Folding@home and will be used for future analyses by our lab and made publicly available at a later date.

These results were unable to be uploaded to fragalysis, but can be found hosted on github at the link below:
github/FoldingAtHome/covid-moonshot/master/folding-at-home_submissions/051920_submission

A full set of our results will also be uploaded to this repository soon, but the ones linked here are all non-covalent Moonshot submissions which match our criteria for convergence by this submission deadline.

The criteria is a Wang-Landau increment of less than 0.12 for the receptor-ligand systems (details can be found in the last forum post), and a negative computed free energy of binding.


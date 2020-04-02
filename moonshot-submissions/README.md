# Docked COVID Moonshot molecules and reference fragment sets

This directory contains docked COVID Moonshot user submissions and screened fragments for comparison.
For ensemble docking methodology, see the [scripts directory](../scripts).

![ensemble of docked molecules](https://github.com/foldingathome/covid-moonshot/raw/master/moonshot-submissions/docked-molecules.png "Ensemble of docked molecules")

## Manifest and current files
* `covid_submissions_03_31_2020-docked.{csv,sdf,pdb}` - COVID Moonshot molecules as of 18:36 PST March 31, 2020 docked into user-specified design fragment structures, with only the best scores (over all fragment structures) preserved
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID (CID)
  * `Hybrid2`: docking score (lower is better)
  * `fragments`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in `../receptors/Mpro-{fragment}-protein.pdb`
* `covid_submissions_03_31_2020.csv` - COVID Moonshot molecules as of 18:36 PST March 31, 2020
* `attic/` - older files
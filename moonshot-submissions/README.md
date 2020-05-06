# COVID Moonshot molecules docked and scored using OEDocking ensemble docking

This directory contains docked COVID Moonshot user submissions and screened fragments for comparison.
For ensemble docking methodology, see the [scripts directory](../scripts).

# All submissions as of 2020-05-03

All submissions as of [2020-05-03 commit](https://github.com/postera-ai/COVID_moonshot_submissions/blob/709a4e1d82e9dd6c12c2df811decdef41d33c559/covid_submissions_all_info.csv) docked to [Mpro monomer with CYS145(-) HIS41(+)](https://github.com/FoldingAtHome/covid-moonshot/tree/master/receptors/monomer):
* `covid_submissions_all_info_ensemble-oedocking.sdf` - formatted for [Fragalysis SDF spec](https://discuss.postera.ai/t/providing-computed-poses-for-others-to-look-at/1155/8)
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID (CID)
  * `Chemgauss4`: docking score (lower is better)
  * `docked_fragment`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in [`../receptors/monomers/Mpro-{fragment}-protein-thiolate.pdb`](https\
://github.com/FoldingAtHome/covid-moonshot/tree/master/receptors/monomer)
  * `covalent_distance_min`, `covalent_distance_mean`, `covalent_distance_stderr`: minimum, mean, and stderr of the mean distance (in A) between covalent warhead heavy atom and CYS145 SG
* `covid_submissions_all_info_ensemble-docked.csv` - same as above, but with all fields preserved (except Rationale, which contains problematic characters)
* `covid_submissions_all_info_ensemble-docked.{sdf,pdb}` - same as above, but with docked poses and all tags (in SDF)
* `covid_submissions_all_info_ensemble-docking-justscores.csv` - summary, sorted by hits
* `covid_submissions_all_info.csv` - all compound SMILES and CIDs as of 2020-05-03 from [repo](https://github.com/postera-ai/COVID_moonshot_submissions/blob/709a4e1d82e9dd6c12c2df811decdef41d33c559/covid_submissions_all_info.csv) with Rationale field deleted (since it contains prolematic characters) and the SMILES entry for `J__-UNK-57ae77f4-1` had to be edited from `CCCCCCCCCC(~C~[O-]C(C)C)C1=CC(P)=CC2CCNC12` to remove the problematic `~` characters in order to get OpenEye to read the whole file.


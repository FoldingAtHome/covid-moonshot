# Docked COVID Moonshot molecules and reference fragment sets

This directory contains docked COVID Moonshot user submissions and screened fragments for comparison.
For ensemble docking methodology, see the [scripts directory](../scripts).

![ensemble of docked molecules](https://github.com/foldingathome/covid-moonshot/raw/master/moonshot-submissions/docked-molecules.png "Ensemble of docked molecules")

# All submissions (covalent and noncovalent) as of 2020-04-23
* `covid_submissions_all_info_ensemble-oedocking.sdf` - formatted for [Fragalysis SDF spec](https://discuss.postera.ai/t/providing-computed-poses-for-others-to-look-at/1155)
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID (CID)
  * `Chemgauss4`: docking score (lower is better)
  * `docked_fragment`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in `../receptors/Mpro-{fragment}-protein.pdb`
  * `covalent_distance_min`, `covalent_distance_mean`, `covalent_distance_stderr`: minimum, mean, and stderr of the mean distance (in A) between covalent warhead atom and CYS145 SG
* `covid_submissions_all_info_ensemble-docked.csv` - same as above, but with all fields preserved (except Rationale, which contains problematic characters)
* `covid_submissions_all_info_ensemble-docked.{sdf,pdb}` - same as above, but with docked poses and all tags (in SDF)
* `covid_submissions_all_info_ensemble-docking-justscores.csv` - summary, sorted by hits
* `covid_submissions_all_info-2020-04-23.csv` - initial list from [repo](https://github.com/mc-robinson/COVID_moonshot_submissions/blob/master/covid_submissions_all_info.csv) with Rationale field deleted (since it contains prolematic characters)

## Final round 2 compounds 
* `COVID_MS_final_selection_round_2-docked-justscores.csv` - scored/docked compounds from covalent round 2
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID (CID)
  * `Chemgauss4`: docking score (lower is better)
  * `docked_fragment`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in `../receptors/Mpro-{fragment}-protein.pdb`
  * `covalent_distance_min`, `covalent_distance_mean`, `covalent_distance_stderr`: minimum, mean, and stderr of the mean distance (in A) between covalent warhead atom and CYS145 SG
* `2020_04_07_Nir_covalent_filtered_and_rejects_Holly_7_April-docked.csv` - same as above, but with all fields preserved (except Rationale, which contains problematic characters)
* `2020_04_07_Nir_covalent_filtered_and_rejects_Holly_7_April-docked.{sdf,pdb}` - same as above, but with docked poses and all tags (in SDF)
* `COVID_MS_final_selection_round_2.csv` - initial list

## Constrained covalent docking of all compounds with covalent warheads as of 2020_04_07 from Holly, using thiolate form of CYS145
* `2020_04_07_Nir_covalent_filtered_and_rejects_Holly_7_April-docked-justscores.csv` - noncovalent docking with just summary scores
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID (CID)
  * `Hybrid2`: docking score (lower is better)
  * `fragments`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in `../receptors/Mpro-{fragment}-protein.pdb`
* `2020_04_07_Nir_covalent_filtered_and_rejects_Holly_7_April-docked.csv` - same as above, but with all fields preserved (except Rationale, which contains problematic characters)
* `2020_04_07_Nir_covalent_filtered_and_rejects_Holly_7_April-docked.{sdf,pdb}` - same as above, but with docked poses and all tags (in SDF)
* `2020_04_07_Nir_covalent_filtered_and_rejects_Holly_7_April.csv` - initial list

## Noncovalent docking of all compounds as of 6 Apr 2020
* `covid_submissions_all_info-2020-04-06-docked-justscores.csv` - noncovalent docking with just summary scores
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID (CID)
  * `Hybrid2`: docking score (lower is better)
  * `fragments`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in `../receptors/Mpro-{fragment}-protein.pdb`
* `covid_submissions_all_info-2020-04-06-docked.csv` - same as above, but with all fields preserved (except Rationale, which contains problematic characters)
* `covid_submissions_all_info-2020-04-06-docked.{sdf,pdb}` - same as above, but with docked poses and all tags (in SDF)
* `covid_submissions_all_info-2020-04-06.csv` - All COVID Moonshot molecules updated 16:25 PST April 6, 2020 from [main github repo](https://discuss.postera.ai/t/updated-list-of-all-submissions/17/1)
One line had to be fixed:
```
CCCCCCCCCC(~C~[O-]C(C)C)C1=CC(P)=CC2CCNC12,J-UNK-57a-1,j,x0820,FALSE,FALSE,FALSE,https://covid.postera.ai/covid/submissions/57ae77f4-67f2-42f0-80f2-cf3ab5890fea,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE
```
to
```
CCCCCCCCCC(C[O-]C(C)C)C1=CC(P)=CC2CCNC12,J-UNK-57a-1,j,x0820,FALSE,FALSE,FALSE,https://covid.postera.ai/covid/submissions/57ae77f4-67f2-42f0-80f2-cf3ab5890fea,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE
```
in order to get OpenEye to read the whole file.

## Constrained docking of covalent inhibitors

Note that this approach does not always place the covalent warhead close to CYS145---I'm still debugging why this may be the case.

* `covid_submissions_with_warhead_info-docked-justscores.csv` - COVID Moonshot molecules as of Round 2 close on 2 Apr 2020, docked into user-specified design fragment structures, with only the best scores (over all fragment structures) preserved
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID (CID)
  * `Hybrid2`: docking score (lower is better)
  * `fragments`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in `../receptors/Mpro-{fragment}-protein.pdb`
* `covid_submissions_with_warhead_info-docked.csv` - same as above, but with all fields preserved (except Rationale, which contains problematic characters)
* `covid_submissions_with_warhead_info-docked.{sdf,pdb}` - same as above, but with docked poses and all tags (in SDF)
* `attic/` - older files

## Older files
* `attic/covid_submissions_03_31_2020-docked.{csv,sdf,pdb}` - COVID Moonshot molecules as of 18:36 PST March 31, 2020 docked into user-specified design fragment structures, with only the best scores (over all fragment structures) preserved
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID (CID)
  * `Hybrid2`: docking score (lower is better)
  * `fragments`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in `../receptors/Mpro-{fragment}-protein.pdb`
* `attic/covid_submissions_03_31_2020.csv` - COVID Moonshot molecules as of 18:36 PST March 31, 2020
* `attic/` - older files

# Synthetic enumeration for free energy calculations

Modeling of design conformations suggested by this [excellent blog post from Pat Walters](http://practicalcheminformatics.blogspot.com/2020/03/building-on-fragments-from-diamondxchem_30.html) and the [accompanying code](https://github.com/PatWalters/fragment_expansion/tree/master/oechem_eval).

## Manifest

### Compound sets
* `primary_amine_enumeration_for_chodera_lab_FEP.csv` - primary amine series (843 compounds)
* `boronic_ester_enumeration_for_chodera_lab_FEP.csv` - boronic ester series (122 compounds)
* `nucleophilic_displacement_enumeration_for_FEP.csv` - nucleophilic displacement series (15918)

### Docked conformers

* `primary_amine_enumeration_for_chodera_lab_FEP-permuted-conformers-x10789.sdf.gz`
* `boronic_ester_enumeration_for_chodera_lab_FEP-permuted-conformers-x10789.sdf.gz`
* `nucleophilic_displacement_enumeration_for_FEP-permuted-conformers-x10789.sdf.gz` - nucleophilic displacement series ()

### Scripts
* `01-fix-csv-files.sh` - permute columns of input files
* `02-generate-poses.py` - generate constrained input poses for a single fragment structure: `x10789` (`TRY-UNI-2eddb1ff-7`)

### Calculation metadata
* `2020-07-24.json`: `primary_amine_enumeration_for_chodera_lab_FEP.csv` forward only built from `x2646`
* `2020-07-27.json`: `primary_amine_enumeration_for_chodera_lab_FEP.csv` and `boronic_ester_enumeration_for_chodera_lab_FEP.csv` forward only built from `x10789`
* `2020-07-28.json`: `primary_amine_enumeration_for_chodera_lab_FEP.csv` and `boronic_ester_enumeration_for_chodera_lab_FEP.csv` backward only built from `x10789`
* `activity-data-2020-07-29.csv`: activity data downloaded from the [COVID Moonshot](https://covid.postera.ai/covid/activity_data.csv) on 2020-07-29
* `activity-data-2020-07-29-conformers-x10789.sdf.gz`: strictly filtered 3-aminopyridine related set for retrospective testing, with activity data preserved as SD tags (40 compounds)

## Procedure

Given a single reference fragment structure (already in `../receptors/` and prepared by `../scripts/00-prep-all-receptors.py`):
* expand uncertain stereochemistry for all target molecules
* identify the **common core** shared by the bound fragment structure and all target molecules
* identify the most likely protonation state in solution for each target molecule
* densely enumerate conformers with Omega, constraining the common core positions to the bound fragment structure
* pick the conformer with the least clashes with protein atoms

The SDF file written out contains the protonated fragment as molecule 0 followed by all docked fragments.

## TODO

* Don't use MMFF for omega:
```
Warning: OEMMFFParams::PrepMol() : unable to type atom 21 N
```

# Synthetic enumeration for free energy calculations

Modeling of design conformations suggested by this [excellent blog post from Pat Walters](http://practicalcheminformatics.blogspot.com/2020/03/building-on-fragments-from-diamondxchem_30.html) and the [accompanying code](https://github.com/PatWalters/fragment_expansion/tree/master/oechem_eval).

## Manifest

### Compound sets

* `primary_amine_enumeration_for_chodera_lab_FEP.csv` - primary amine series (843 compounds)
* `boronic_ester_enumeration_for_chodera_lab_FEP.csv` - boronic ester series (122 compounds)
* `nucleophilic_displacement_enumeration_for_FEP.csv` - nucleophilic displacement series (15918)
* `activity-data-2020-07-29.csv` - all compounds with activity data for retrospective benchmarking (888)
* `aminopyridine_compounds_for_FEP_benchmarking.csv` - 3-aminopyridine retrospective benchmarking compounds from Matthew Robinson (70)
* `fastgrant-table1.csv` - sentinel cases from Fast Grant application (Alpha Lee) (11)
* `RAL-THA-6b94ceba` - Ralph Robinson 3-aminopyridine P4 pocket exploration [`RAL-THA-6b94ceba`](https://postera.ai/covid/submissions/6b94ceba-f352-4275-ad8d-e766e56e6fa4)

### Docked conformers
* `primary_amine_enumeration_for_chodera_lab_FEP-permuted-conformers-x10789.sdf.gz`
* `boronic_ester_enumeration_for_chodera_lab_FEP-permuted-conformers-x10789.sdf.gz`
* `nucleophilic_displacement_enumeration_for_FEP-permuted-conformers-x10789.sdf.gz` - nucleophilic displacement series ()
* `aminopyridine_compounds_for_FEP_benchmarking-conformers-x10789.sdf` - 3-aminopyridine retrospective benchmarking compounds using single common fragment, prioritized by docked scores
* `aminopyridine_compounds_for_FEP_benchmarking-dockscores-x10789.sdf` - 3-aminopyridine retrospective benchmarking compounds using individual maximum common fragment enumeration, prioritized by docked scores

### Scripts
* `01-fix-csv-files.sh` - permute columns of input files
* `02-generate-poses.py` - generate constrained input poses for a single fragment structure: `x10789` (`TRY-UNI-2eddb1ff-7`)

### Calculation metadata
* `2020-07-24.json`: `primary_amine_enumeration_for_chodera_lab_FEP.csv` forward only built from `x2646`
* `2020-07-27.json`: `primary_amine_enumeration_for_chodera_lab_FEP.csv` and `boronic_ester_enumeration_for_chodera_lab_FEP.csv` forward only built from `x10789`
* `2020-07-28.json`: `primary_amine_enumeration_for_chodera_lab_FEP.csv` and `boronic_ester_enumeration_for_chodera_lab_FEP.csv` backward only built from `x10789`
* `2020-07-29-retrospective-aminopyridines.json`: retrospective 3-aminopyridine compounds that JDC fished out automatically, single common scaffold, sparse conformers ranked by steric clashes
* `2020-08-02-retrospective-aminopyridines-matt.json`: retrospective 3-aminopyridine compounds from Matthew Robinson, using a single common scaffold based on x10789, dense conformers ranked by dock scores (-Cl in wrong pocket)
* `2020-08-03-retrospective-aminopyridines-matt-dockscores.json`: retrospective 3-aminopyridine compounds from Matthew Robinson, using a individually selected common scaffolds based on x10789, dense conformers ranked by dock scores (-Cl in right pocket)

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

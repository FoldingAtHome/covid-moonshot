# Synthetic enumeration for free energy calculations

Modeling of design conformations suggested by this [excellent blog post from Pat Walters](http://practicalcheminformatics.blogspot.com/2020/03/building-on-fragments-from-diamondxchem_30.html) and the [accompanying code](https://github.com/PatWalters/fragment_expansion/tree/master/oechem_eval).

## Manifest
* `boronic_ester_enumeration_for_chodera_lab_FEP.csv` - boronic ester series
* `primary_amine_enumeration_for_chodera_lab_FEP.csv` - primary amine series
* `01-fix-csv-files.sh` - permute columns of input files
* `02-generate-poses.py` - generate constrained input poses for a single fragment structure: `x10789` (`TRY-UNI-2eddb1ff-7`)

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

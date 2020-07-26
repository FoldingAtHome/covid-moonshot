# Synthetic enumeration for free energy calculations

Modeling of design conformations suggested by this [excellent blog post from Pat Walters](http://practicalcheminformatics.blogspot.com/2020/03/building-on-fragments-from-diamondxchem_30.html) and the [accompanying code](https://github.com/PatWalters/fragment_expansion/tree/master/oechem_eval).

## Manifest
* `boronic_ester_enumeration_for_chodera_lab_FEP.csv` - boronic ester series
* `primary_amine_enumeration_for_chodera_lab_FEP.csv` - primary amine series
* `01-fix-csv-files.sh` - permute columns of input files
* `02-generate-poses.py` - generate constrained input poses for a single fragment structure: `x10789` (`TRY-UNI-2eddb1ff-7`)

## To generate models

* Download and unpack latest [Moonshot PDB files](https://fragalysis.diamond.ac.uk/media/targets/Mpro.zip)
* Prep receptor for x2646 (aminopyridine fragment) with `python 00-prep-all-receptors.py`

# Molecule sets for the COVID Moonshot

## Manifest
* `covid_submissions_03_26_2020.csv` - submissions to https://covid.postera.ai/covid downloaded from https://discuss.postera.ai/t/updated-list-of-all-submissions/17 as of 2:34 PST March 27, 2020, retaining only SMILES, CID, and list of fragments used as inspiration
* `covid_submissions_03_26_2020.sdf` - conversion of the above to SDF with `fragments` SD tag
* `melatonin-and-metabolites.csv` - melatonin and its metabolites
* `clinical-meleatonin-receptor-agonists.sdf` -
* `enamine-melatonin-realspace-analogues.sdf` -
* `Mpro-fragments.csv` - fragment SMILES and crystal names (`Mpro-x####`) extracted from `../diamond-structures/Mpro full XChem screen - experiment summary - ver-2020-03-25-annotated.xlsx`

## Notes
* Made a copy of [COVID Moonshot submissions sheet](https://discuss.postera.ai/t/updated-list-of-all-submissions/17), deleted all but first two columns (SMILES and TITLE), and exported to CSV
* Converted CSV to SDF files for upload to Orion:
```bash
convert.py clinical-melatonin-receptor-agonists.csv clinical-melatonin-receptor-agonists.sdf
convert.py melatonin-and-metabolites.csv melatonin-and-metabolites.sdf
convert.py enamine-melatonin-realspace-analogues.csv enamine-melatonin-realspace-analogues.sdf
convert.py covid_submissions_03_22_2020.csv covid_submissions_03_22_2020.sdf
```

# Broad Institute repurposing library ensemble docking screen

## Manifest

### Broad Drug Repurposing Hub library (2020-03-24)

Dataset from the [Broad Drug Repurposing Hub Library](https://clue.io/repurposing) as of 2020-03-24

* `broad-repurposing-docked.csv` - docked summary statistics
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID (CID)
  * `Hybrid2`: docking score (lower is better)
  * `fragments`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in `../receptors/Mpro-{fragment}-protein.pdb`
* `broad-repurposing-docked.sdf` - SDF files of top hits, with `fragments` SD tag indicating which structure they are docked to
* `broad-repurposing-library-20200324.{txt,csv}` - Broad Institute repurposing library downloaded from https://clue.io/repurposing#download-data

### Drugbank 5.1.5 (2020-01-03)

Dataset from [Drugbank](https://www.drugbank.ca/releases/latest#structures) 5.1.5 (2020-01-03)

* `drugbank-5.1.5-2020-01-03-docked.csv` - docked summary statistics
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID (CID)
  * `Hybrid2`: docking score (lower is better)
  * `fragments`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in `../receptors/Mpro-{fragment}-protein.pdb`
* `drugbank-5.1.5-2020-01-03-docked.sdf` - docked structures
* `drugbank-5.1.5-2020-01-03.csv`

## Method

See `../scripts/README.md`

```bash
# Iterate over all molecules
for JOBID in {1..10147}
do
  python ../scripts/02-dock-and-prep.py --molecules broad-repurposing-library-20200324.csv --receptors ../receptors --output repurposing-screen-docked --index $JOBID
done
# Aggregate data
python ../scripts/03-aggregate-docking-results.py --docked repurposing-screen-docked --output broad-repurposing-docked.csv --clean
python ../scripts/03-aggregate-docking-results.py --docked repurposing-screen-docked --output broad-repurposing-docked.sdf
```

# Broad Institute repurposing library ensemble docking screen

## Manifest
* `broad-repurposing-docked.csv` - top hits
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID (CID)
  * `Hybrid2`: docking score (lower is better)
  * `fragments`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in `../receptors/Mpro-{fragment}-protein.pdb`
* `broad-repurposing-docked.sdf` - SDF files of top hits, with `fragments` SD tag indicating which structure they are docked to
* `broad-repurposing-library-20200324.{txt,csv}` - Broad Institute repurposing library downloaded from https://clue.io/repurposing#download-data

## Method
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

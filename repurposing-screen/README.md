# Drug repurposing library ensemble docking screens

This directory contains the results of a quick docking screen to identify potentially useful molecules with FDA approval for some other indication, as well as other molecules of clinical interest, that may bind to the active site or dimer interface of SARS-CoV-2 main viral protease.

Currently, only the Broad repurposing library (below) has been fully docked.

## Manifest

### Broad Drug Repurposing Hub library (2020-03-24)

Dataset from the [Broad Drug Repurposing Hub Library](https://clue.io/repurposing) as of 2020-03-24

* `broad-repurposing-docked.csv` - docked summary statistics
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID (CID)
  * `Hybrid2`: docking score (lower is better)
  * `fragments`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in [`../receptors/monomer/Mpro-{fragment}-protein.pdb`](https://github.com/FoldingAtHome/covid-moonshot/tree/master/receptors/monomer)
* `broad-repurposing-docked.sdf` - SDF file containing docked poses of top hits, with `fragments` SD tag indicating which structure they are docked to
* `broad-repurposing-docked.{pdb,mol2}` - PDB and Tripos mol2 versions docked poses from `broad-repurposing-docked.sdf` without annotations
* `broad-repurposing-library-20200324.{txt,csv}` - undocked Broad Institute repurposing library downloaded from https://clue.io/repurposing#download-data

To inspect the top hits, please view `broad-repurposing-docked.csv`.

To view the docked poses, please first load [`../receptors/monomer/Mpro-x0336-protein.pdb`](`broad-repurposing-docked.csv`) (you may need to click the `Raw` button to download the PDB file to your machine) and load the `broad-repurposing-docked.{sdf,pdb,mol2}` file in your favorite format.

### Drugbank 5.1.5 (2020-01-03)

This docking screen is still in progress.

Dataset from [Drugbank](https://www.drugbank.ca/releases/latest#structures) 5.1.5 (2020-01-03)

* `drugbank-5.1.5-2020-01-03-docked.csv` - docked summary statistics
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID (CID)
  * `Hybrid2`: docking score (lower is better)
  * `fragments`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in `../receptors/Mpro-{fragment}-protein.pdb`
* `drugbank-5.1.5-2020-01-03-docked.sdf` - docked structures
* `drugbank-5.1.5-2020-01-03.csv`

## Method

See `../scripts/README.md` for complete details of methodology and an index of the scripts.

All structure preparation and docking utilized the [OpenEye Toolkits version 2019.10.2](https://docs.eyesopen.com/toolkits/python/releasenotes/releasenotes/index.html#release-highlights-2019-oct).


### Structure preparation

Each fragment-bound X-ray structure provided by DiamondMX was prepared using the [OpenEye Spruce Toolkit](https://docs.eyesopen.com/toolkits/python/sprucetk/index.html) to clean up the structure and extract completed protein and small molecule structures.
Protein structures with covalently bound molecules had the covalent adduct removed and their labeled CYS residues reverted to a standard cysteine.

Only the protease monomer was used in this docking screen.
All structures were prepared in `../receptors/monomer/`.


### Docking

For each molecule to be docked, a set of reasonable protonation and tautomeric states were enumerated using [`OEGetReasonableProtomers`](https://docs.eyesopen.com/toolkits/python/quacpactk/OEProtonFunctions/OEGetReasonableProtomers.html#OEProton::OEGetReasonableProtomers) from the [OpenEye QuacPac Toolkit](https://docs.eyesopen.com/toolkits/python/quacpactk/index.html).
Each protonation/tautomeric state was then expanded into a dense set of conformers using the [OpenEye Omega Toolkit](https://docs.eyesopen.com/toolkits/python/omegatk/index.html).
For each active site (non-covalent or covalent) and dimer interface fragment structure, one docked pose was generated using the binding site defined by the corresponding fragment using the corresponding X-ray structure.
Each multiconformer molecule was docked to the specified fragment site using the [OpenEye OEDocking Toolkit](https://docs.eyesopen.com/toolkits/python/dockingtk/index.html) using the [`Hybrid2` search method](https://docs.eyesopen.com/toolkits/python/dockingtk/docking.html#hybrid-method) (which exploits shape overlays with the corresponding bound fragment) and the [`High`](https://docs.eyesopen.com/toolkits/python/dockingtk/OEDockingConstants/OESearchResolution.html#OEDocking::OESearchResolution::High) search resolution.
The highest dock score across all docked structures was recorded, and the corresponding fragment structure noted.
The resulting hit list reports the best (most negative) docking score for each molecule, sorted from best to worst.

### Commands used

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

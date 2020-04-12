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

#### Primary conclusions so far

##### Vonoprozan

The top docking hit in the Broad repurposing screen is the proton pump inhibitor (PPI) [vonoprzzan](https://www.drugbank.ca/drugs/DB11739), developed by Takeda, used in the treatment of acid reflux diseases like GERD.
While most PPIs are prodrugs that undergo acid-catalyzed activation and then irreversibly inhibit the gastric proton pump enzyme H+,K+-ATPase by covalently attaching to specific cysteine residues, vonoprazan is the first member of a new class of potassium-competitive acid blocker (P-CAB) that instead reversibly competes for K+ ion binding [[DOI](https://doi.org/10.1007/s40265-015-0368-z)].

Recently, an [experimental repurposing screen](https://www.biorxiv.org/content/10.1101/2020.04.03.023846v1) of the [Prestwick Chemical Library of FDA-approved drugs](http://www.prestwickchemical.com/libraries-screening-lib-pcl.html) identified vonoprazan as having more potency in restroing cell viability 3 days following SARS-CoV-2 infection than the arbidol (a broad spectrum antiviral that blocks viral entry of many enveloped viruses) when applied at 10 uM concentration to VeroE6 cell lines[[DOI](https://www.biorxiv.org/content/10.1101/2020.04.03.023846v1)].
However, since three PPIs (omeprazole, chloroquine diphosphate, and hydroxychloroquine sulfate) were also identified as having antiviral activity, leading to the suggested hypothesis that these compounds may act by increasing the pH of the endosomal/golgian pathway and thereby limit processing of the viral spike protein by endosomal proteases, thereby blocking viral entry, rather than through direct binding of viral targets.
While vonoprozan does have measurable biodistribution to the lung, it appears to be to a much lesser degree than the target organ (stomach) [[DOI](https://doi.org/10.1038/s41401-019-0353-2)].

#### Fragment overlap

All the fragments were used to define a set of Gaussians "colored" by interaction types using the [OpenEye Shape Toolkit](https://docs.eyesopen.com/toolkits/python/shapetk/index.html).
For each docked compound in the repurposing set, the shape/color combination overlay between the docked compound and the merged fragments---which corrects for the overlap among fragment atoms to avoid overcounting---was used to measure how well the docked molecule fills the space spanned by the fragments.
The [Tversky measure](https://docs.eyesopen.com/toolkits/python/shapetk/shape_theory.html#molecular-shape) of the merged fragment query was used to determine how much of the space spanned by the fragments was filled by the docked molecule, and is reported as the `overlap_score`.
The results are sorted by this overlap score, which goes from ~1 (fills much of the fragment space) to ~0 (fills none of the fragment space).
For each docked molecule, we also identified all fragments that overlap significantly with the docked molecule using the Tversky measure for the fragment to see what fraction of the fragment overlapped with the docked molecule.
If this was >0.6, the fragment was added to the overlapping_fragments list.

See the following files:

* `broad-repurposing-docked-overlap/` - directory containing each docked repurposing compound with its highly-overlapping fragments, sorted by the fraction of the fragment space spanned (best to worst)
* `broad-repurposing-docked-overlap.csv` - docked compounds, sorted from best to worst fragment overlap
* `broad-repurposing-docked-overlap.sdf` - docked structures with annotations as SD tags, sorted from best to worst overlap score
* `broad-repurposing-docked-overlap.pdb` - docked structures, sorted from best to worst overlap score
* `broad-repurposing-docked-overlap.mol2` - docked structures, sorted from best to worst overlap score

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

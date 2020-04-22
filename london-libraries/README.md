# pyridine urea noncovalent series from Nir London

## Manifest
* `pyridine_urea-docked-justscore.csv` - simplified ranking of hybrid docked compounds to x0434
  * `SMILES`: SMILES for compound
  * `TITLE`: the compound ID
  * `Chemgauss4`: docking score (lower is better)
  * `docked_fragment`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in `../receptors/monomer/Mpro-{fragment}-protein.pdb`
* `pyridine_urea-docked.csv` - same as above, but with all fields preserved
* `pyridine_urea-docked.{sdf,pdb}` - same as above, but with docked poses and all tags (in SDF)
* `pyridine_urea.csv` - JDC conversion to CSV, addition of x0434 inspiration fragment
* `pyridine_urea.smi` - original file from Nir London

## Provenance (Nir London)

15,279 non-covalent pyridine-ureas based on the fragment from x0434.

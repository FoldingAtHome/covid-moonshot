### Sprint 2 synthetic purchasing

Results files accessed from folding at home server 10am BST 25-08-2020.  
3363 molecules analysed (4018 expanded stereoisomers) from second covid-moonshot sprint. 


### Identifying compounds for purchasing
Compounds filtered to purchase molecules with predicted affinity improvements of 1 kcal/mol, relative to molecule `TRY-UNI-2eddb1ff-7`

* Molecule `EN300-1940535` removed, as furan less preferable than comparible (in both structure and predicted affinity `EN300-330332`

Following molecules were removed due to alkyl-bromide residues
* `EN300-1274653`
* `EN300-70105` 
* `EN300-1915724`
* `EN300-1909681`
* `EN300-1274638`
* `EN300-1910987`
* `EN300-701252`
* `EN300-1274639`
* `EN300-1914939`

Molecules flagged for consideration are those containing ester functional group

| Title | Smiles | Notes |
|-------|--------|-------|
|EN300-1886131 | Cc1ccncc1NC(=O)Cc2cc(cc(c2)Cl)OCc3c(nc(s3)NC(=O)OC(C)(C)C)C |       |
|EN300-843777 |CCOC(=O)C(Cn1cc(nn1)COc2cc(cc(c2)Cl)CC(=O)Nc3cnccc3C)(F)F | Neighbouring F's will activate ester - don't make |
|EN300-60436 |CCOC(=O)c1c2ccccc2oc1COc3cc(cc(c3)Cl)CC(=O)Nc4cnccc4C | |
|EN300-301686 |Cc1ccncc1NC(=O)Cc2cc(cc(c2)Cl)OCc3c4ccccc4n(n3)C(=O)OC(C)(C)C| |

For the other esters - Melissa said that they would likely be hydrolysed, Ed Griffen didn't seem to think there was much cause for concern, and Matt Robinson was on the fence, so they are being submitted for purchasing (EN300-1886131, EN300-60436, EN300-301686), but if they show good affinity we can follow up.

### Manifest

* `analysis.json` output generated for sprint 2, using the automated folding@home pipeline `https://github.com/choderalab/covid-moonshot-infra`. Contains per-run data in both 'details' and 'analysis'
* `ligands.csv` Details of predicted high affinity ligands from sprint
* `purchasing.csv` an abridged version of `ligands.csv` containing only ligands proposed for purchasing
* `ligands.pdf` overview of ligands predicted affinity and 2D structure, ordered from predicted highest affinity first
* `base.pse` pymol session for viewing lowest-reverse-work structures from the simulations for highest-affinity compounds, alongside the crystallographic reference
* `ligands.mol2`, `ligands.sdf` ligands in different file formats
* `proteins.pdb` snapshot of protein from lowest-reverse-work structures. The order of these aligns with those in the `ligands.sdf` files
* `reference.mol2`, `reference.sdf` conformation of the 'reference' `TRY-UNI-2eddb1ff-7` molecule from the folding at home simulation, again the order of these aligns with `ligands.sdf`

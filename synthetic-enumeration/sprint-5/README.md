# Sprint 5 : P1' substituents optimization for isoquinoline / benzopyran series

## Objectives

This sprint evaluates potential P1' substituents building on the isoquinoline series (formerly 3-the aminopyridine series) to attempt to further optimize lead compound [MAT-POS-b3e365b9-1](http://postera.ai/covid/submissions/MAT-POS-b3e365b9-1).

### Relevant compounds

Current best designs in the isoquinoline series with activity data posted to https://covid.postera.ai/covid/activity_data
* [MAT-POS-b3e365b9-1](https://postera.ai/covid/submissions/MAT-POS-b3e365b9-1) - P2 chromane substituent
* [EDJ-MED-e4b030d8-11](https://postera.ai/covid/submissions/EDJ-MED-e4b030d8-11) - P2 chromane
* [MAT-POS-0c8fa4a7-1](https://postera.ai/covid/submissions/MAT-POS-0c8fa4a7-1) - P2 pocket dihydronaphtalene
* [ALP-POS-477dc5b7-2](https://postera.ai/covid/submissions/MAT-POS-0c8fa4a7-1) - P2 pocket tetrahydroquinoline
* [ALP-POS-3b848b35-2](https://postera.ai/covid/submissions/ALP-POS-3b848b35-2) - P4 beta-lactam (suspected liability)

Molecules with data in CDD that has not yet been posted to https://postera.ai/covid/activity_data
* [ALP-POS-477dc5b7-3](https://postera.ai/covid/submissions/ALP-POS-477dc5b7-3) - chromane with possible P1 substituent
* [EDJ-MED-e4b030d8-13](https://postera.ai/covid/submissions/EDJ-MED-e4b030d8-13) - chromane with methyl group pointing at P1

### Relevant X-ray structures

These structures appear to contain the common isoquinoline-based scaffold:

* [MAT-POS-b3e365b9-1](https://postera.ai/covid/submissions/MAT-POS-b3e365b9-1), the chirally-separated enantiomer of [VLA-UCB-1dbca3b4-15](https://postera.ai/covid/submissions/VLA-UCB-1dbca3b4-15) responsible for the dominant observe activity - [x11498](https://fragalysis.diamond.ac.uk/viewer/react/preview/direct/target/Mpro/mols/x11498/L/P/C)
* [MAT-POS-8a69d52e-7](https://covid.postera.ai/covid/submissions/MAT-POS-8a69d52e-7) - [x12073](https://fragalysis.diamond.ac.uk/viewer/react/preview/direct/target/Mpro/mols/x12073/L/P/C)

A superposition of these structures can be viewed on Fragalysis [here](https://fragalysis.diamond.ac.uk/viewer/react/projects/295/235).

Observations:
* His163: His163 ND deprotonated to hydrogen bond with Tyr161 OH; either His163 NE protonated and isoquinoline scaffold N deprotonated, or vice-versa to form hydrogen bond
* His164: His164 ND likely protonated to orient chromane Cl; NE likely protonated to donate hydrogen bond to Thr175 OG
* Gln189 NE-H2 hydrogen bond with chromane O
* linker C=O hydrogen bonds with protein backbone Glu166 NH
* Protonation states of Cys145 indeterminate: Cys145 SG- could interact with linker N and His41 NE protonated; His41 ND not interacting with anything but protonated could help stabilize ligand -Cl substituent

Required scaffold interactions:
* isoquinoline N - His163 NE
* linker C=O - Glu166 backbone N

Other scaffold interactions:
* isoquinoline `c1cn[c:1]c2ccccc12` - Phe140 O
* isoquinoline `c1cn[c:1]c2ccccc12` - Phe140 O

Other interactions of note:
* chromane O - Gln189 NE
* chromane Cl substituent - His164 ND
* chromane Cl substituent - His164 ND

Other Fragalysis structures of interest:

#### Benzopyrans
* [JAG-UCB-119787ef-1](https://postera.ai/covid/submissions/JAG-UCB-119787ef-1) - [x10898](https://fragalysis.diamond.ac.uk/viewer/react/preview/direct/target/Mpro/mols/x10898/L/P/C)
* [MAT-POS-968e8d9c-1](https://postera.ai/covid/submissions/MAT-POS-968e8d9c-1) - [x10906](https://fragalysis.diamond.ac.uk/viewer/react/preview/direct/target/Mpro/mols/x10906/L/P/C)
* [LON-WEI-0a73fcb8-7](https://postera.ai/covid/submissions/LON-WEI-0a73fcb8-7) - [x10942](https://fragalysis.diamond.ac.uk/viewer/react/preview/direct/target/Mpro/mols/x10942/L/P/C)
* [EDJ-MED-0e996074-1](https://postera.ai/covid/submissions/EDJ-MED-0e996074-1) - [x11159](https://fragalysis.diamond.ac.uk/viewer/react/preview/direct/target/Mpro/mols/x11159/L/P/C)
* [BRU-THA-a358fbdd-2](https://postera.ai/covid/submissions/RU-THA-a358fbdd-2) - [x11233](https://fragalysis.diamond.ac.uk/viewer/react/preview/direct/target/Mpro/mols/x11233/L/P/C)


## Plan

Our initial goal is to determine which of the 10 distinct intermediates we should make in bulk:
* 2 of these are commercial
* the other 8 are derived from 2 other intermediates

We aim to cluster related compounds to reduce the variance in our estimates.

## Experiment design

The design of this experiment is intended to help us quickly identify which intermediates are worth pursuing at scale.
We therefore sort the molecules we could generate by number of atoms and interleave molecules from every series,
working from small to large.

* Transformations between X-ray structures
* Transformations from X-ray structures to all compounds that have been assayed or queued for synthesis
* Transformations from X-ray structures to all proposed intermediates
* Transformations from all proposed intermediates to subsets (small transformations, good docking scores) of each intermediate, interleaved by intermediate
* Transformations from X-ray structures to all submitted designs not queued for synthesis
* Transformations from all proposed intermediates to all possible designs [CLUSTERED?]

## Notes

We just got data for these, and should be sure to include them:

Synthetic intermediates:
![Sprint 5 synthetic routes](sprint-5-synthetic-routes.jpg)

## Manifest

Scripts
* `01-aggregate-compounds.py` - generate the initial sorted list to dock
* `02-generate-poses.py` - generate poses
* `03-sort-poses-key-interactions.py` - sort poses by key interactions (used for Sprint 5)
* `03-sort-poses-docking-score.py` - sort poses by docking score
* `04-create-json.py` - create JSON for setting up perses calculations

Docked compounds
* `sprint-5-compounds.sdf` - compounds metadata with annotations
* `sprint-5-microstates-x11498.sdf` - docked microstates to frag x11498 structure
* `sprint-5-microstates-x12073.sdf` - docked microstates to frag x12073 structure
* `sprint-5-microstates-x11498-inters.sdf` - annotations from Tim Dudgeon
* `sprint-5-microstates-x12073-inters.sdf` - annotations from Tim Dudgeon
* `sprint-5-microstates-x11498-sorted.sdf` - sorted by number of fragment interactions and clashes
* `sprint-5-microstates-x12073-sorted.sdf` - sorted by number of fragment interactions and clashes

SMARTS strings to label compounds by source intermediates
* `intermediates/benzopyran_sprint_5_intermediates.csv`

Filtered for liabilities using Eli Lilly rules:
* `filtered/filtered_alkyl_halide_bb_ether_synthesis_for_FEP.smi`
* `filtered/filtered_amine_bb_amide_couplings_for_FEP.smi`
* `filtered/filtered_cooh_bb_amide_couplings_for_FEP.smi`

Submitted designs
* `submissions/submissions-2020-11-02.csv` - all submitted designs as of 2020-11-02

Sorted, compiled, annotated molecules
* `sorted/sprint-5.csv` - molecules to be docked

## Misc notes on filtering

Here are the filtered lists with a unique ID for every compound. Since there were so many molecules, I used Lilly's rules, which are easier to run in high-throughput (since I couldn't do all these by hand). You'll notice that the number of "demerits": "D()" is still included in the final files, which you can use if you want to shorten the list further. Of those that were not filtered out, the main demerit is simply having many atoms. Not sure, but could be a decent predictor of FEP challenges.

A group of commands like this should give you what you need in Pandas:
```
df = pd.read_csv('filtered_alkyl_halide_bb_ether_synthesis_for_FEP.smi', sep=' ', header=None)
df = df.rename(columns={0: 'SMILES', 1: 'ID', 3:'demerits', 4:'reason_for_demerits'})
df['demerits'] = df['demerits'].apply(lambda x: int(x.split('(')[-1].split(')')[0]))
df = df.drop(columns=2)
```

Larger D() score is indeed worse. According to Lilly's rules, everything >= D(100) is filtered out.

## Methodology

**Compounds evaluated:** Enumeration of new compound designs based was performed based on synthetic routes proposed by the med chem design team. In addition, all compound designs manually submitted in the COVID Moonshot website as of 2020-11-02 containing the amide linker with isoquinoline scaffold (SMARTS match `C(=O)Nc1cncc2ccccc12`) were included.

**Microstate enumeration:** Any atoms with undefined stereochemistry were enumerated with `OEOmega.Flipper` from the OpenEye Toolkit 2020.1.0, creating multiple enantiomers that are associated with each racemic design. A single reasonable tautomer was selected for each molecule using `OEGetReasonableProtomer`.

**Initial pose generation:** `OEOmega` was used with dense sampling to generate a dense set of conformations with fixed scaffold (either the benzopyran-linker-isoquinoline, `C1(CCOc2ccccc12)C(=O)Nc1cncc2ccccc12`, or if that fails, the liker-isoquinoline, `C(=O)Nc1cncc2ccccc12`) from the following X-ray structures where the protein and reference ligand were prepared with OpenEye SpruceTK in the monomeric form of Mpro using neutral His41 and Cys145. This structure had His163 was protonated to hydrogen bond with the isoquinoline scaffold N.
* [MAT-POS-b3e365b9-1](https://postera.ai/covid/submissions/MAT-POS-b3e365b9-1), the chirally-separated enantiomer of [VLA-UCB-1dbca3b4-15](https://postera.ai/covid/submissions/VLA-UCB-1dbca3b4-15) responsible for the dominant observe activity - [x11498](https://fragalysis.diamond.ac.uk/viewer/react/preview/direct/target/Mpro/mols/x11498/L/P/C)
* [MAT-POS-8a69d52e-7](https://covid.postera.ai/covid/submissions/MAT-POS-8a69d52e-7) - [x12073](https://fragalysis.diamond.ac.uk/viewer/react/preview/direct/target/Mpro/mols/x12073/L/P/C)
Relative alchemical free energy network: A star map was used to create a map from the reference X-ray structure to all proposed designs.

**Alchemical free energy calculations:** Nonequilibrium relative alchemical free energy calculations were prepared with [`perses`](http://github.com/choderalab/perses) 0.8.1 (based on OpenMM) was used with an addition of an RMSD restraint to keep the “core” ligand heavy atoms in the atom map and protein alpha carbons within 6.5A of ligand heavy atoms that turns on after the RMSD exceeds 2 A (with force constant K_RMSD ~ 1 kT/A^2) (see [this PR](https://github.com/choderalab/perses/pull/748)). The AMBER14SB force field was used the protein, TIP3P with AMBER-recommended counterion parameters for solvent, and the Open Force Field Initiative `openff-1.2.0` (“Parsley”) small molecule force field. Na+ and Cl- ions were added to achieve an effective ionic strength of 70 mM to match experimental conditions. No bond constraints were used for non-water molecules (due to a performance issue with constraints in certain alchemical transformations in [OpenMM](http://openmm.org) 7.4.2), and hydrogen mass repartitioning was used (hydrogen mass 4 amu) with timestep of 2 fs. A Monte Carlo barostat was employed with temperature of 300 K and pressure of 1 atm, and the BAOAB Langevin integrator was used, with an initial equilibration of 500 ps. Nonequilibrium free energy calculations were run on Folding@home, executing cycles of 500 ps equilibration in ligand A, 500 ps nonequilibrium switching from ligand A to ligand B, 500 ps equilibration in B, and 500 ps switching from B to A. The nonequilibrium protocol work is accumulated for forward (A->B) and backward (B->A) nonequilibrium switches. The Bennett acceptance ratio (BAR) is used to estimate the free energy difference between molecules. 40 independent trajectories of 2 nonequilibrium cycles/trajectory were used to collect work values for the ligand transformations in complex, and 25 trajectories of 2 cycles/trajectory for the ligand in solvent.

**Known issues:** For some `compounds` with defined stereochemistry, it appears that the single generated `microstate` sometimes has inverted stereochemistry for reasons we don't understand. See, for example, `ALP-POS-d3acb8cc-1` in `json/sprint-5-x12073-monomer-neutral.json`.

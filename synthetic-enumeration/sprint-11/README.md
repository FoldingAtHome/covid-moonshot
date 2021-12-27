# Sprint 11

Retrospective and prospective assessment of designs for optimizing sprio compounds with 5-membered fused rings, focusing on

## Reference structures used in this sprint

* `P1800_0A` : crystallized with [VLA-UCB-50c39ae8-2](https://covid.postera.ai/covid/submissions/VLA-UCB-50c39ae8-2) (IC50 246 nM) - base for chromane-5spiro-isoquinoline : `O=C1C[C@]2(CCOc3ccc(Cl)cc32)C(=O)N1c1cncc2ccccc12` - conformation 1, site A; reference ligand will be [VLA-UCB-50c39ae8-2](https://covid.postera.ai/covid/submissions/VLA-UCB-50c39ae8-2)

* `P2113_0B` : crystallized with [JOH-MSK-1f2dff76-2](https://covid.postera.ai/covid/submission/JOH-MSK-1f2dff76-2) - base for THIQ-5spiro-isoquinoline : enantiopure [MAT-POS-1bed62cf-3](https://covid.postera.ai/covid/submissions/MAT-POS-1bed62cf-3) - conformation 1, site B; reference ligand will be [JOH-MSK-1f2dff76-4](https://covid.postera.ai/covid/submission/JOH-MSK-1f2dff76-4) to avoid having to truncate the P3/P4 substituent in every transformation

* `P2222_0A` : [`MAT-POS-c7726e07-5`](https://covid.postera.ai/covid/submissions/MAT-POS-c7726e07-5) (IC50 74 nM) - base for THIQ-5spiro-isqoquinolines : enantiopure active form of [`EDJ-MED-8bb691af-4`](https://covid.postera.ai/covid/submissions/EDJ-MED-8bb691af-4) - conformation 2, site A; reference ligand will be [JOH-MSK-1f2dff76-4](https://covid.postera.ai/covid/submission/JOH-MSK-1f2dff76-4) to avoid having to truncate the P3/P4 substituent in every transformation

## Notes

This sprint was run on Folding@home using [perses](http://github.com/choderalab/perses) commit [`7d073a9`](https://github.com/choderalab/perses/commit/7d073a9dab1dd9f857c8f2e4b3eaf996ebd17a53) and periodic nonequilibrium cycling with 4 fs timesteps, 1 ns for each phase of 4 ns cycles, following 5 ns of equilibration.
Atom mappings were derived from docked geometries using the use_given_geometries=True options, which uses the `AtomMapper.generate_mapping_from_positions()` facility.

Both restrained and non-restrained calculations were run.

## Manifest
* `receptors/` - generated with [`fah_prep`](https://github.com/choderalab/fah_prep) commit [`94b25d5`](https://github.com/choderalab/fah_prep/commit/94b25d53303e5b9e924d2e18dc162406042ac6ef).
* `submissions/` - submissions as of 2021-12-26
* `sorted/` - filtered submissions

## Steps
```
for fragid in P1800_0A P2113_0B P2222_0A; do
prepare receptors -i ../../structures/Mpro/aligned/Mpro-${fragid}/ -f "Mpro-${fragid}*.pdb" -o ./receptors
done
```

## Notes

Make sure we get predictions for [EDJ-MED-fd8ed875](https://covid.postera.ai/covid/submissions/fd8ed875-7ab8-4dfb-a695-3f1dd5feea36) and [EDJ-MED-a12e3a20](https://covid.postera.ai/covid/submissions/a12e3a20-6154-48de-9476-1a7403650f51).


Spirocycles:
P2222 A B - tetrahydroisoquinoline (P2)
P2385 A B - identical to P2222?

P2113 A
P2113 B : displaced from A, alternative conformers that decouple His41 - Cys145  ( [MAT-POS-1bed62cf-3](https://covid.postera.ai/covid/submissions/1bed62cf-6fb2-4954-b5e3-72f9cb131639/3)  - 222 nM racemate) - matches bound pose of conformation of P1800 (quinoline-chromane)

P1800 A B - chromane (P2)


P2222 and P2385 curiously have no density


P1507 A B - benzofuran

P0143 A B - chromane-5spiro-isoquinoline

P0207 A
P2385 A B - THIQ-5spiro-quinoline slightly offset, similar to P2222
P2402 A B
P0207 A B
P0765 A B
P0811 A - chromane-6spiro-quinoline (six-membered spiro) in same pose as P2222
P0950 A B - THIQ-6spiro-quinoline (six-membered spiro) in same pose as P2222

P0207 A B - chromane-5spiro-iqoquinoline; A in same pose as P1800, B intermediate between P1800 and P2222

P2222 A B (quinoline-5spiro-THIQ) is significantly displaced toward P2 from P1800 (quinoline-chromane)

Start with P2222
fix quinoline ring
build in rest of spiro
keep protonation states consistent

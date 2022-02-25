# Sprint 12

Retrospective and prospective assessment of designs for optimizing sprio compounds with 5-membered fused rings, including all designs up to 9 Feb 2022

**NOTE:** This used the conda-forge perses 0.9.3, which was missing key FAH generator improvements. These improvements will be included in 0.9.4.

## Reference structures used in this sprint

* `P2385_0A` : crystallized with [MAT-POS-c7726e07-5](https://covid.postera.ai/covid/submissions/MAT-POS-c7726e07-5) (IC50 95 [81, 111] nM)

## Notes

### Protonation and tautomer expansion

Expansion of protonation states and computation of state penalties uses Schrodinger Epik:
```
# Convert single-conformer SDF to MAE
$SCHRODINGER/utilities/structconvert input.sdf input.mae
# Expand protonation states with state penalty threshold of ~6 kT
$SCHRODINGER/epik -WAIT -ms 16 -ph 7.3 -p 0.002 -imae input.mae -omae output.mae
# Convert expanded protonation states back to SDF to extract molecule and state penalty
$SCHRODINGER/epik -WAIT -ms 16 -ph 7.3 -p 0.002 -imae output.mae -omae output.sdf
```
State penalty is stored in `r_epik_State_Penalty` field (in kcal/mol).

This sprint was run on Folding@home using [perses](http://github.com/choderalab/perses) commit [`7d073a9`](https://github.com/choderalab/perses/commit/7d073a9dab1dd9f857c8f2e4b3eaf996ebd17a53) and periodic nonequilibrium cycling with 4 fs timesteps, 1 ns for each phase of 4 ns cycles, following 5 ns of equilibration.
Atom mappings were derived from docked geometries using the use_given_geometries=True options, which uses the `AtomMapper.generate_mapping_from_positions()` facility.

Both restrained and non-restrained calculations were run.

## Manifest
* `receptors/` - generated with [`fah_prep`](https://github.com/choderalab/fah_prep) commit [`94b25d5`](https://github.com/choderalab/fah_prep/commit/94b25d53303e5b9e924d2e18dc162406042ac6ef).
* `submissions/` - submissions
* `sorted/` - filtered submissions

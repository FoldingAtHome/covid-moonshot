# Sprint 11

Retrospective and prospective assessment of designs for optimizing sprio compounds with 5-membered fused rings, focusing on

## Reference structures used in this sprint

* `P1800_0A` : crystallized with [VLA-UCB-50c39ae8-2](https://covid.postera.ai/covid/submissions/VLA-UCB-50c39ae8-2) (IC50 246 nM) - base for chromane-5spiro-isoquinoline : `O=C1C[C@]2(CCOc3ccc(Cl)cc32)C(=O)N1c1cncc2ccccc12` - conformation 1, site A; reference ligand will be [VLA-UCB-50c39ae8-2](https://covid.postera.ai/covid/submissions/VLA-UCB-50c39ae8-2)


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

## Issues

Had to drop this structure because Spruce did not model residues 45-51 in chain B in a simulatable fashion

* `P2113_0B` : crystallized with [JOH-MSK-1f2dff76-2](https://covid.postera.ai/covid/submission/JOH-MSK-1f2dff76-2) - base for THIQ-5spiro-isoquinoline : enantiopure [MAT-POS-1bed62cf-3](https://covid.postera.ai/covid/submissions/MAT-POS-1bed62cf-3) - conformation 1, site B; reference ligand will be [JOH-MSK-1f2dff76-4](https://covid.postera.ai/covid/submission/JOH-MSK-1f2dff76-4) to avoid having to truncate the P3/P4 substituent in every transformation

It produced:
```
ATOM   2687  N   CYS B  44      13.574  -9.223  27.323  1.00 57.76           N
ATOM   2688  CA  CYS B  44      13.544  -8.135  28.304  1.00 60.17           C
ATOM   2689  C   CYS B  44      12.181  -7.962  28.981  1.00 61.28           C
ATOM   2690  O   CYS B  44      11.609  -6.862  28.952  0.00 99.99           O
ATOM   2691  CB  CYS B  44      14.003  -6.831  27.661  1.00 62.48           C
ATOM   2692  SG  CYS B  44      15.718  -6.853  27.097  1.00 70.92           S
ATOM   2693  OXT CYS B  44      11.598  -9.015  29.614  0.00 99.99           O
ATOM   7296  H   CYS B  44      13.281  -9.015  26.390  1.00 20.00           H
ATOM   7297  HA  CYS B  44      14.268  -8.377  29.096  1.00 20.00           H
ATOM   7298  HB2 CYS B  44      13.357  -6.626  26.795  1.00 20.00           H
ATOM   7299  HB3 CYS B  44      13.892  -6.024  28.400  1.00 20.00           H
ATOM   7300  HG  CYS B  44      16.197  -7.883  27.729  1.00 20.00           H
ATOM   2694  N   PRO B  52      19.928  -2.973  29.545  1.00 79.97           N
ATOM   2695  CA  PRO B  52      19.833  -4.438  29.554  1.00 79.99           C
ATOM   2696  C   PRO B  52      21.116  -5.111  29.091  1.00 79.95           C
ATOM   2697  O   PRO B  52      21.722  -4.686  28.109  1.00 80.34           O
ATOM   2698  CB  PRO B  52      18.673  -4.732  28.591  1.00 80.65           C
ATOM   2699  CG  PRO B  52      18.573  -3.530  27.718  1.00 81.05           C
ATOM   2700  CD  PRO B  52      19.094  -2.364  28.491  1.00 79.50           C
ATOM   7301  H2  PRO B  52      19.628  -2.626  30.433  1.00 20.00           H
ATOM   7302  H3  PRO B  52      20.881  -2.714  29.391  1.00 20.00           H
ATOM   7303  HA  PRO B  52      19.571  -4.797  30.560  1.00 20.00           H
ATOM   7304  HB2 PRO B  52      18.889  -5.628  27.990  1.00 20.00           H
ATOM   7305  HB3 PRO B  52      17.737  -4.882  29.149  1.00 20.00           H
ATOM   7306  HG2 PRO B  52      19.175  -3.676  26.809  1.00 20.00           H
ATOM   7307  HG3 PRO B  52      17.523  -3.356  27.439  1.00 20.00           H
ATOM   7308  HD2 PRO B  52      19.697  -1.707  27.847  1.00 20.00           H
ATOM   7309  HD3 PRO B  52      18.268  -1.787  28.932  1.00 20.00           H
```

* `P2222_0A` : [`MAT-POS-c7726e07-5`](https://covid.postera.ai/covid/submissions/MAT-POS-c7726e07-5) (IC50 74 nM) - base for THIQ-5spiro-isqoquinolines : enantiopure active form of [`EDJ-MED-8bb691af-4`](https://covid.postera.ai/covid/submissions/EDJ-MED-8bb691af-4) - conformation 2, site A; reference ligand will be [JOH-MSK-1f2dff76-4](https://covid.postera.ai/covid/submission/JOH-MSK-1f2dff76-4) to avoid having to truncate the P3/P4 substituent in every transformation

`receptors/dimer/Mpro-P2222_0A_bound-protein.pdb` cannot be simulated because it's missing 44-52 of chain B:
```
ATOM   2679  N   ILE B  43      13.465 -11.100  25.079  1.00 49.09           N
ATOM   2680  CA  ILE B  43      13.936 -11.647  26.351  1.00 51.80           C
ATOM   2681  C   ILE B  43      13.922 -10.588  27.451  1.00 52.80           C
ATOM   2682  O   ILE B  43      14.952 -10.354  28.100  0.00 99.99           O
ATOM   2683  CB  ILE B  43      15.322 -12.309  26.200  1.00 53.42           C
ATOM   2684  CG1 ILE B  43      16.455 -11.250  26.074  1.00 55.26           C
ATOM   2685  CG2 ILE B  43      15.318 -13.303  25.035  1.00 53.83           C
ATOM   2686  CD1 ILE B  43      17.858 -11.806  25.806  1.00 57.18           C
ATOM   2687  OXT ILE B  43      12.778  -9.903  27.720  0.00 99.99           O
ATOM   7272  H   ILE B  43      14.132 -10.818  24.390  1.00 20.00           H
ATOM   7273  HA  ILE B  43      13.236 -12.438  26.656  1.00 20.00           H
ATOM   7274  HB  ILE B  43      15.513 -12.881  27.121  1.00 20.00           H
ATOM   7275 HG12 ILE B  43      16.195 -10.574  25.246  1.00 20.00           H
ATOM   7276 HG13 ILE B  43      16.491 -10.680  27.015  1.00 20.00           H
ATOM   7277 HG21 ILE B  43      16.312 -13.765  24.943  1.00 20.00           H
ATOM   7278 HG22 ILE B  43      14.566 -14.085  25.221  1.00 20.00           H
ATOM   7279 HG23 ILE B  43      15.072 -12.774  24.102  1.00 20.00           H
ATOM   7280 HD11 ILE B  43      18.574 -10.974  25.737  1.00 20.00           H
ATOM   7281 HD12 ILE B  43      18.151 -12.474  26.629  1.00 20.00           H
ATOM   7282 HD13 ILE B  43      17.855 -12.368  24.861  1.00 20.00           H
ATOM   2688  N   ASN B  53      21.977  -6.061  29.521  1.00 65.29           N
ATOM   2689  CA  ASN B  53      23.020  -7.024  29.191  1.00 65.89           C
ATOM   2690  C   ASN B  53      22.340  -8.235  28.602  1.00 65.87           C
ATOM   2691  O   ASN B  53      21.910  -9.120  29.341  1.00 65.84           O
ATOM   2692  CB  ASN B  53      23.803  -7.402  30.445  1.00 67.81           C
ATOM   2693  CG  ASN B  53      25.075  -8.166  30.162  1.00 71.36           C
ATOM   2694  OD1 ASN B  53      25.132  -9.055  29.296  1.00 72.30           O
ATOM   2695  ND2 ASN B  53      26.125  -7.830  30.898  1.00 72.09           N
ATOM   7283  H1  ASN B  53      21.206  -6.177  28.895  1.00 20.00           H
ATOM   7284  H2  ASN B  53      21.669  -6.211  30.460  1.00 20.00           H
ATOM   7285  H3  ASN B  53      22.340  -5.133  29.434  1.00 20.00           H
ATOM   7286  HA  ASN B  53      23.710  -6.595  28.450  1.00 20.00           H
ATOM   7287  HB2 ASN B  53      24.065  -6.478  30.981  1.00 20.00           H
ATOM   7288  HB3 ASN B  53      23.158  -8.025  31.081  1.00 20.00           H
```

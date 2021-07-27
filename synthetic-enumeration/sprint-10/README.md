# Sprint 10

Retrospective and prospective assessment of designs replacing the P1 substituent based on 'x10959' : 'ADA-UCB-6c2cb422-1' including all designs and data prior to 2021-07-26.

This sprint was run on Folding@home using [perses](http://github.com/choderalab/perses) commit [`7d073a9`](https://github.com/choderalab/perses/commit/7d073a9dab1dd9f857c8f2e4b3eaf996ebd17a53) and periodic nonequilibrium cycling with 4 fs timesteps, 1 ns for each phase of 4 ns cycles, following 5 ns of equilibration.
Atom mappings were derived from docked geometries using the use_given_geometries=True options, which uses the `AtomMapper.generate_mapping_from_positions()` facility.

Both restrained and non-restrained calculations were run.

## Manifest
* `receptors/` - generated with [`fah_prep`](https://github.com/choderalab/fah_prep) commit [`94b25d5`](https://github.com/choderalab/fah_prep/commit/94b25d53303e5b9e924d2e18dc162406042ac6ef).

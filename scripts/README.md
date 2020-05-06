# Hybrid docking

This script uses the [OpenEye OEDocking toolkit](https://docs.eyesopen.com/toolkits/python/dockingtk/index.html) to perform hybrid docking of small molecules.

## Manifest
* `00-prep-all-receptors.py` - prepare `OEReceptor` files from all DiamondMX structures using [Spruce](https://docs.eyesopen.com/toolkits/python/sprucetk/index.html); note that generated PDB files for the protein contain `UNK` atoms in the case of covalent inhibitors---these must be filtered out
* `01-check-receptors.py` - check that `Mpro-*-protein.pdb` files can be assigned AMBER14SB parameters if `UNK` atoms are filtered out
* `02-dock-and-prep.py` - perform ensemble docking of a single compound and (optionally) prepare it for simulation
* `03-aggregate-docking-results.py` - compile docking results
* `attic/` - older scripts

## Ensemble docking protocol

Each fragment-bound X-ray structure provided by DiamondMX was prepared using the [OpenEye Spruce Toolkit](https://docs.eyesopen.com/toolkits/python/sprucetk/index.html) to clean up the structure and extract protein and small molecule structures with missing atoms added and missing loops either reconstructed or truncated.
CYS145 was deprotonated and the hydrogens reoptimized, which produced a doubly-protonated HIS41.

For each COVID Moonshot designed molecule, a set of reasonable protonation and tautomeric states were enumerated using [`OEGetReasonableProtomers`](https://docs.eyesopen.com/toolkits/python/quacpactk/OEProtonFunctions/OEGetReasonableProtomers.html#OEProton::OEGetReasonableProtomers) from the [OpenEye QuacPac Toolkit](https://docs.eyesopen.com/toolkits/python/quacpactk/index.html).
Each protomer was then expanded into a dense set of conformers using the [OpenEye Omega Toolkit](https://docs.eyesopen.com/toolkits/python/omegatk/index.html).
For each active site (non-covalent or covalent) and dimer interface fragment structure, one docked pose was generated using the binding site defined by the corresponding fragment using the corresponding X-ray structure.

Each multiconformer molecule was docked to the specified fragment site using the [OpenEye OEDocking Toolkit](https://docs.eyesopen.com/toolkits/python/dockingtk/index.html) using the [`Hybrid2` search method](https://docs.eyesopen.com/toolkits/python/dockingtk/docking.html#hybrid-method) (which exploits shape overlays with the corresponding bound fragment) and the [`High`](https://docs.eyesopen.com/toolkits/python/dockingtk/OEDockingConstants/OESearchResolution.html#OEDocking::OESearchResolution::High) search resolution.
The highest dock score across all docked structures was recorded, and the corresponding fragment structure noted.

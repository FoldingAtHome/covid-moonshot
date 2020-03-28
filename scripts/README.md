# Hybrid docking

This script uses the [OpenEye OEDocking toolkit](https://docs.eyesopen.com/toolkits/python/dockingtk/index.html) to perform hybrid docking of small molecules.

## Manifest
* `00-prep-all-receptors.py` - prepare `OEReceptor` files from all DiamondMX structures using [Spruce](https://docs.eyesopen.com/toolkits/python/sprucetk/index.html); note that generated PDB files for the protein contain `UNK` atoms in the case of covalent inhibitors---these must be filtered out
* `01-check-receptors.py` - check that `Mpro-*-protein.pdb` files can be assigned AMBER14SB parameters if `UNK` atoms are filtered out
* `02-dock-ligands-to-corresponding-receptors.py` - dock small molecule ligands to their corresponding receptor structures
* `02-dock-ligands-to-corresponding-receptors-multiprocessing.py` - dock small molecule ligands to their corresponding receptor structures
* `03-dock-ligands-to-one-receptor.py` - dock small molecule ligands to a pre-specified receptor (optional utility)

## Docking protocol

Each fragment-bound X-ray structure provided by DiamondMX was prepared using the [OpenEye Spruce Toolkit](https://docs.eyesopen.com/toolkits/python/sprucetk/index.html) to clean up the structure and extract completed protein and small molecule structures.

For each COVID Moonshot designed molecule, a set of reasonable protonation and tautomeric states were enumerated using [`OEGetReasonableProtomers`](https://docs.eyesopen.com/toolkits/python/quacpactk/OEProtonFunctions/OEGetReasonableProtomers.html#OEProton::OEGetReasonableProtomers) from the [OpenEye QuacPac Toolkit](https://docs.eyesopen.com/toolkits/python/quacpactk/index.html).
Each protomer was then expanded into conformers using the [OpenEye Omega Toolkit](https://docs.eyesopen.com/toolkits/python/omegatk/index.html).
For each fragment structure the Moonshot compound design listed as inspiration, one docked pose was generated using the binding site defined by the specified inspiration fragment.
Each multiconformer molecule was docked to the specified fragment site using the [OpenEye OEDocking Toolkit](https://docs.eyesopen.com/toolkits/python/dockingtk/index.html) using the [`Hybrid2` search method](https://docs.eyesopen.com/toolkits/python/dockingtk/docking.html#hybrid-method) (which exploits shape overlays with the corresponding bound fragment) and the [`High`](https://docs.eyesopen.com/toolkits/python/dockingtk/OEDockingConstants/OESearchResolution.html#OEDocking::OESearchResolution::High) search resolution.

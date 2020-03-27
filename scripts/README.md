# Hybrid docking

This script uses the [OpenEye OEDocking toolkit](https://docs.eyesopen.com/toolkits/python/dockingtk/index.html) to perform hybrid docking of small molecules.

## Manifest
* `01-prep-all-receptors.py` - prepare `OEReceptor` files from all DiamondMX structures using [Spruce](https://docs.eyesopen.com/toolkits/python/sprucetk/index.html)
* `02-dock-ligands-to-corresponding-receptors.py` - dock small molecule ligands to their corresponding receptor structures
* `03-dock-ligands-to-one-receptor.py` - dock small molecule ligands to a pre-specified receptor

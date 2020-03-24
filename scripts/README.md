# Hybrid docking

This script uses the [OpenEye OEDocking toolkit](https://docs.eyesopen.com/toolkits/python/dockingtk/index.html) to perform hybrid docking of small molecules.

## Manifest
* `01-prep-receptor.py` - prepare `OEReceptor` files from DiamondMX structures using [Spruce](https://docs.eyesopen.com/toolkits/python/sprucetk/index.html)
* `02-dock-ligands.py` - dock small molecule ligands

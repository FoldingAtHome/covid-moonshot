# Older unused docking scripts

## Manifest
* `02-dock-ligands-to-corresponding-receptors.py` - dock small molecule ligands to their corresponding receptor structures
* `02-dock-ligands-to-corresponding-receptors-multiprocessing.py` - dock small molecule ligands to their corresponding receptor structures
* `02-dock-ligands-to-all-receptors-orion.py` - upload receptors to Orion for parallel docking IN THE CLOUD!
* `03-dock-ligands-to-one-receptor.py` - dock small molecule ligands to a pre-specified receptor (optional utility)
* `04-parameterize-ligands.py` - pre-generate JSON parameter caches for [`openmmforcefields.generators.SystemGenerator`](https://github.com/openmm/openmmforcefields#automating-force-field-management-with-systemgenerator) in `json-files/`
* `05-combine-all-docking-results.py` - combine all Orion datarecord results into CSV and SDF files

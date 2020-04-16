# SARS-CoV known inhibitors

See https://discuss.postera.ai/t/a-brief-exploration-of-past-sars-small-molecule-inhibitors/895?u=johnchodera

## Manifest

* `SARS_2003_inhibitors-docked-justscore.csv` - docked SARS-CoV inhibitors compiled by Matt Robinson
  * `SMILES`: SMILES for compound
  * `TITLE`: the paper_id
  * `Hybrid2`: docking score (lower is better)
  * `docked_fragment`: fragment ID for corresponding fragment X-ray structure for best docked pose; corresponding structures are in `../receptors/monomer/Mpro-{fragment}-protein.pdb`
* `SARS_2003_inhibitors-docked.csv` - docked SARS-CoV inhibitors compiled by Matt Robinson, with all annotation preserved
* `SARS_2003_inhibitors-docked.sdf` - docked SARS-CoV inhibitors compiled by Matt Robinson, with all annotation preserved
* `SARS_2003_inhibitors.csv` - SARS-CoV inhibitors compiled by Matt Robinson (columns permuted so SMILES is first)

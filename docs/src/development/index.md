# Development notes

The detailed record of how the pipeline was tested and changed: what each
real-data test found, how it was diagnosed, and why the code is the way it
is. The rest of the documentation gives the conclusions; these pages give the
evidence, for anyone checking a result or changing the code it concerns.

| Notes | Covers |
|:--|:--|
| [ZTF investigation log](ztf-investigation.md) | the first real-data tests: five bugs, a clouded frame, the quality gate, plate solving, photometry, variable-star false positives |
| [Pan-STARRS1 working notes](pan-starrs1-notes.md) | every failure on real Pan-STARRS1 image sets, the configuration now in `pipeline_config`, epochs, completeness, the blind test |
| [Design decisions](design-decisions.md) | choices revisited: the PSF fallback, ZOGY's astrometric noise term, performance of reference building, technologies tried and rejected |

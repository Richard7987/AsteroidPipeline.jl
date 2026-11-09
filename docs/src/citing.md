# How to cite

If you use AsteroidPipeline.jl in scientific work, please cite it.
[`CITATION.cff`](https://github.com/Richard7987/AsteroidPipeline.jl/blob/main/CITATION.cff)
has the reference, which GitHub's "Cite this repository" button gives as APA
or BibTeX:

```bibtex
@software{asteroidpipeline,
  author  = {Bernal Estrada, Alejandro},
  title   = {{AsteroidPipeline.jl}},
  version = {1.5.0},
  year    = {2026},
  doi     = {10.5281/zenodo.23871668},
  url     = {https://github.com/Richard7987/AsteroidPipeline.jl},
}
```

That DOI always leads to the latest version; each version has its own, listed on
[Zenodo](https://doi.org/10.5281/zenodo.23871668), for citing the exact one you used.

## What it builds on

Please also cite what your results depend on.

| Used for | Reference |
|:--|:--|
| astrometric and photometric calibration: Gaia DR3 | [Gaia Collaboration 2016](https://doi.org/10.1051/0004-6361/201629272), [2023](https://doi.org/10.1051/0004-6361/202243940) |
| Gaia through VizieR, and VSX | [Ochsenbein et al. 2000](https://doi.org/10.1051/aas:2000169) |
| Gaia through NOIRLab's Astro Data Lab | [Nikutta et al. 2020](https://doi.org/10.1016/j.ascom.2020.100411) |
| known objects: SkyBoT | [Berthier et al. 2016](https://doi.org/10.1093/mnras/stw492) |
| cross-matching: SIMBAD | [Wenger et al. 2000](https://doi.org/10.1051/aas:2000332) |
| cross-matching: VSX | Watson, Henden & Price 2006, Society for Astronomical Sciences Annual Symposium 25, 47 |
| asteroid magnitudes: the H-G model | Bowell et al. 1989, in *Asteroids II*, University of Arizona Press, 524–556 |
| analytic PSF: Moffat profile | Moffat 1969, Astronomy & Astrophysics 3, 455 |
| difference imaging: ZOGY | [Zackay, Ofek & Gal-Yam 2016](https://doi.org/10.3847/0004-637X/830/1/27) |
| plate solving: astrometry.net | [Lang et al. 2010](https://doi.org/10.1088/0004-6256/139/5/1782) |
| NEO score: digest2 | [Keys et al. 2019](https://doi.org/10.1088/1538-3873/ab1157) |
| periodogram: Lomb–Scargle | [Lomb 1976](https://doi.org/10.1007/BF00648343), [Scargle 1982](https://doi.org/10.1086/160554) |
| rotation periods: Fourier series of the phase | [Harris et al. 1989](https://doi.org/10.1016/0019-1035(89)90015-8) |
| orbits: Gauss's method, Lambert's problem, universal-variable propagation | Curtis 2014, *Orbital Mechanics for Engineering Students*, 3rd ed., Butterworth-Heinemann |
| orbits: systematic ranging | [Virtanen, Muinonen & Bowell 2001](https://doi.org/10.1006/icar.2001.6592) |
| orbits: rejecting outlying observations | [Carpino, Milani & Chesley 2003](https://doi.org/10.1016/S0019-1035(03)00051-4) |
| world coordinates in FITS | [Greisen & Calabretta 2002](https://doi.org/10.1051/0004-6361:20021326) |
| the language | [Bezanson et al. 2017](https://doi.org/10.1137/141000671) |

If your images come from a survey, cite it as well. The validation here used
Pan-STARRS1 ([Chambers et al. 2016](https://arxiv.org/abs/1612.05560)) and ZTF
([Bellm et al. 2019](https://doi.org/10.1088/1538-3873/aaecbe),
[Masci et al. 2019](https://doi.org/10.1088/1538-3873/aae8ac)). Surveys and
catalogues also ask for an acknowledgement in their own words: check their
terms of use.

The Minor Planet Center's orbits and Confirmation Page, and JPL's Horizons and
Small-Body Identification services, are online services: name the ones you
used, with the date you queried them.

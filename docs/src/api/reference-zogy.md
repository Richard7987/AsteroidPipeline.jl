# Reference & ZOGY Differencing

Building a deep static-sky reference stack (`build_reference`,
`load_frame`), estimating a frame's point-spread function empirically or
via an analytic Moffat fallback (`estimate_psf`, `fit_moffat_psf`), and
ZOGY proper image subtraction (`zogy_subtract`).

```@autodocs
Modules = [AsteroidPipeline]
Pages = ["reference.jl", "psf.jl", "zogy.jl"]
```

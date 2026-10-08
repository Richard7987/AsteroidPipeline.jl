# Photometry & Completeness

Gaia-calibrated magnitudes for candidates (`aperture_flux`,
`photometric_zeropoint`, `candidate_magnitudes` — what lets
`mpc80_report` fill in a magnitude the way IASC's Astrometrica reports do),
and synthetic-asteroid injection to measure how complete a configuration
really is (`inject_movers`, `injection_recovery`; see
`examples/injection_test.jl`): a handful of real objects can't trace a
completeness curve, thousands of injected ones can.

```@autodocs
Modules = [AsteroidPipeline]
Pages = ["calibrated_photometry.jl", "injection.jl"]
```

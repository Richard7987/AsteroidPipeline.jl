"""
    crossmatch_catalog(candidates, catalog::Symbol; radius::Real)

Cross-match candidate tracklets against a known-object catalog
(`:skybot`, `:vsx`, or `:simbad`) within `radius` arcseconds, to separate
known objects from candidates warranting human verification.
"""
function crossmatch_catalog(candidates, catalog::Symbol; radius::Real)
    error("not yet implemented")
end

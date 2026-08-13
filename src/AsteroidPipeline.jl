module AsteroidPipeline

using FITSIO
using Photometry
using LombScargle
using HTTP
using CSV
using TypedTables: Table

include("detection.jl")
include("linking.jl")
include("crossmatch.jl")

export detect_sources, link_candidates, crossmatch_catalog

end # module AsteroidPipeline

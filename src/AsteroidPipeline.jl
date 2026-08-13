module AsteroidPipeline

using FITSIO
using Photometry
using LombScargle

include("detection.jl")
include("linking.jl")
include("crossmatch.jl")

export detect_sources, link_candidates, crossmatch_catalog

end # module AsteroidPipeline

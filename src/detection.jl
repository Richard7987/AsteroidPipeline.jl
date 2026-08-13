"""
    detect_sources(image::AbstractMatrix{<:Real}; threshold::Real)

Detect point sources in a background-subtracted FITS frame.

Returns a table of candidate source positions and fluxes, to be consumed
by [`link_candidates`](@ref) for inter-frame motion matching.
"""
function detect_sources(image::AbstractMatrix{<:Real}; threshold::Real)
    error("not yet implemented")
end

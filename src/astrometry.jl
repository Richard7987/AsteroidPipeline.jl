"""
    load_wcs(header::AbstractString) -> WCSTransform

Parse the astrometric solution (WCS) from a raw FITS header string — e.g.
`FITSIO.read_header(hdu, String)` — returning the primary `WCS.WCSTransform`.

Throws an error if the header defines no WCS solution. If it defines more
than one (e.g. alternate WCS descriptions with an `a`/`b`/... suffix), the
first one is returned.
"""
function load_wcs(header::AbstractString)
    solutions = WCS.from_header(String(header))
    isempty(solutions) && error("no WCS solution found in header")
    return first(solutions)
end

"""
    pix_to_sky(wcs::WCSTransform, x::Real, y::Real) -> (ra, dec)

Convert a single 1-based pixel position `(x, y)` — matching both Julia's
array indexing and the FITS convention — to sky coordinates (right
ascension, declination, in degrees) using the astrometric solution `wcs`.
"""
function pix_to_sky(wcs::WCSTransform, x::Real, y::Real)
    world = pix_to_world(wcs, Float64[x, y])
    return (ra=world[1], dec=world[2])
end

"""
    astrometric_calibrate(tracklets, wcs_per_frame, timestamps)

Convert the pixel positions in `tracklets` (as returned by
[`link_candidates`](@ref)) to sky coordinates, producing candidates ready
for [`crossmatch_catalog`](@ref).

`wcs_per_frame[k]` is the `WCS.WCSTransform` describing frame `k`'s
astrometric solution (see [`load_wcs`](@ref)); `timestamps[k]` is that
frame's observation epoch (Julian Date), matching the `timestamps` used
by `link_candidates`.

Returns a flat table with one row per tracklet point, with columns `id`
(the tracklet's index in `tracklets`), `frame`, `x`, `y`, `ra`, `dec`
(degrees), and `epoch` (Julian Date). Since it already carries `id`,
`ra`, `dec`, and `epoch`, this table's rows can be passed directly to
[`crossmatch_catalog`](@ref).
"""
function astrometric_calibrate(tracklets, wcs_per_frame, timestamps)
    id = Int[]
    frame = Int[]
    x = Float64[]
    y = Float64[]
    ra = Float64[]
    dec = Float64[]
    epoch = Float64[]

    for (tracklet_id, tracklet) in enumerate(tracklets), p in tracklet
        sky = pix_to_sky(wcs_per_frame[p.frame], p.x, p.y)
        push!(id, tracklet_id)
        push!(frame, p.frame)
        push!(x, p.x)
        push!(y, p.y)
        push!(ra, sky.ra)
        push!(dec, sky.dec)
        push!(epoch, timestamps[p.frame])
    end

    return Table(; id, frame, x, y, ra, dec, epoch)
end

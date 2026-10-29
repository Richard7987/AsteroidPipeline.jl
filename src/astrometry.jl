"""
    load_wcs(header::AbstractString) -> WCSTransform

Parse the astrometric solution (WCS) from a raw FITS header string — e.g.
`FITSIO.read_header(hdu, String)` — returning the primary `WCS.WCSTransform`.

Throws an error if the header defines no celestial (RA/Dec) WCS solution:
IRAF writes a pixel-to-pixel one (`WCSDIM`, `LTM1_1`…) into frames it
processed, which would otherwise pass for a solution. If it defines more
than one (e.g. alternate WCS descriptions with an `a`/`b`/... suffix), the
first celestial one is returned.

Pan-STARRS1 headers carry `CNPIX1`/`CNPIX2`, a legacy IRAF keyword pair
that wcslib takes for the start of a second, implicit WCS; it builds that
from uninitialised fields and fails with an error that depends on heap
state ("singular", "Invalid parameter value"). On any wcslib error from a
header carrying those keywords, the parse is retried with only those two
removed.

For a frame with no WCS at all, see [`solve_wcs`](@ref), which needs only
the pointing, or [`plate_solve`](@ref) (nova.astrometry.net).
"""
function load_wcs(header::AbstractString)
    solutions = try
        WCS.from_header(String(header))
    catch e
        (e isa ErrorException && _has_fits_card(header, ("CNPIX1", "CNPIX2"))) || rethrow()
        WCS.from_header(_strip_fits_cards(String(header), ("CNPIX1", "CNPIX2")))
    end
    celestial = filter(w -> startswith(w.ctype[1], "RA") && startswith(w.ctype[2], "DEC"), solutions)
    isempty(celestial) && error("no celestial WCS solution found in header")
    return first(celestial)
end

"""
    _strip_fits_cards(header, keywords) -> String

`header` with any 80-character FITS card whose keyword starts with one of
`keywords` removed. Used by [`load_wcs`](@ref) to drop the specific
legacy keywords that trigger a real wcslib parsing bug — see its
docstring.
"""
function _strip_fits_cards(header::AbstractString, keywords)
    ncards = length(header) ÷ 80
    cards = [header[80(i-1)+1:80i] for i in 1:ncards]
    keep = filter(c -> !any(k -> startswith(c, k), keywords), cards)
    return join(keep)
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

# Whether `header` has an 80-character FITS card whose keyword starts with
# one of `keywords` (see `_strip_fits_cards`).
function _has_fits_card(header::AbstractString, keywords)
    ncards = length(header) ÷ 80
    return any(i -> any(k -> startswith(header[80(i-1)+1:80i], k), keywords), 1:ncards)
end

# the keywords frame_epoch looks for, in order, and the first a header has
const _EPOCH_KEYS = ("DATE-OBS", "MJD-OBS", "OBSMJD", "MJD_OBS", "JD", "OBSJD", "JD-OBS")
_first_key(header, keys) = (i = findfirst(k -> haskey(header, k), keys); i === nothing ? nothing : keys[i])

"""
    frame_epoch(header::FITSHeader; timestamp_key=nothing, exptime_key="EXPTIME") -> Float64

A frame's observation epoch as a UTC Julian Date at **mid-exposure** — what
MPC reports require, and what Astrometrica writes into IASC reports.

`timestamp_key` names the header keyword holding the time: a number is read
as a Modified Julian Date (as a Julian Date if the keyword names one, or if
it is above 2,400,000), a string as an ISO-8601 date-time (fractional
seconds kept). Left `nothing`, the first of the usual keywords the header
has is used: `DATE-OBS` (with `TIME-OBS`, `UT` or `UTC` when it holds only
the date), `MJD-OBS`, `OBSMJD`, `MJD_OBS`, `JD`, `OBSJD`, `JD-OBS` — so
images from any telescope that follows one of these conventions need no
setting. `DATE-OBS` first matters for Pan-STARRS1: its `MJD-OBS` is the shutter
time in **TAI**, 37 s later than its `DATE-OBS` (UTC), despite
`TIMESYS = 'UTC'`. ZTF has no `DATE-OBS`; its `OBSMJD` (UTC) is found.

If the header has `exptime_key` (seconds), half of it is added: both surveys
stamp the *start* of the exposure. `exptime_key=nothing` skips that, for a
timestamp that already marks mid-exposure.
"""
function frame_epoch(header::FITSHeader; timestamp_key::Union{Nothing,AbstractString}=nothing,
                     exptime_key::Union{Nothing,AbstractString}="EXPTIME")
    key = timestamp_key === nothing ? _first_key(header, _EPOCH_KEYS) : timestamp_key
    key === nothing && throw(ArgumentError("no observation time in the header (looked for " *
                                           join(_EPOCH_KEYS, ", ") * "); pass timestamp_key"))
    value = header[key]
    if value isa AbstractString && occursin(r"^\s*\d{4}-\d\d-\d\d\s*$", value)
        # a date alone: the time is in a keyword of its own
        tkey = _first_key(header, ("TIME-OBS", "UT", "UTC", "UT-START"))
        tkey === nothing || (value = strip(value) * "T" * strip(string(header[tkey])))
    end
    jd = if value isa AbstractString
        m = match(r"^(\d{4}-\d\d-\d\d(?:T\d\d:\d\d:\d\d)?)(\.\d+)?", strip(value))
        m === nothing && throw(ArgumentError("$key = $(repr(value)) is not an ISO-8601 date-time"))
        datetime2julian(DateTime(m[1])) + (m[2] === nothing ? 0.0 : parse(Float64, m[2])) / 86400
    elseif occursin("JD", key) && !occursin("MJD", key) || value > 2.4e6
        Float64(value)
    else
        Float64(value) + 2400000.5
    end
    if exptime_key !== nothing && haskey(header, exptime_key)
        jd += header[exptime_key] / 2 / 86400
    end
    return jd
end

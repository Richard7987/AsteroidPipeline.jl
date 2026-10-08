"""
    load_wcs(header::AbstractString) -> WCSTransform

Parse the astrometric solution (WCS) from a raw FITS header string — e.g.
`FITSIO.read_header(hdu, String)` — returning the primary `WCS.WCSTransform`.

Throws an error if the header defines no WCS solution. If it defines more
than one (e.g. alternate WCS descriptions with an `a`/`b`/... suffix), the
first one is returned.

Real Pan-STARRS1 (PS1) headers — used by real IASC practice campaigns —
carry `CNPIX1`/`CNPIX2` (a legacy IRAF/DSS-plate-astrometry keyword pair,
unrelated to the header's actual CTYPE/CRVAL/CRPIX/CDELT WCS, present
alongside it) that `WCS.jl`/wcslib mistakes for the start of a *separate*,
implicit DSS-style WCS description; with none of that convention's other
required keywords present, wcslib builds a degenerate, all-zero linear
transform for it and raises "Linear transformation matrix is singular" —
even though the header's real, complete WCS parses fine on its own.
Confirmed directly (bisecting a real PS1 header down to the single
offending keyword): removing just `CNPIX1`/`CNPIX2` is sufficient: the
same header then parses cleanly, and the first (real) solution is
unaffected. Handled here as a targeted retry — only when the header
actually carries those keywords, only stripping those two — rather than
a general parsing workaround.

The retry used to trigger only on wcslib's "singular" message, which
turned out to be just one face of the problem: the implicit solution's
undefined fields are whatever memory wcslib got, so the error it raises
depends on the process's heap state. Measured on a real 2025 PS1 header:
300 parses in a fresh process all failed "singular", but inside a longer
`run_pipeline` call the same header failed "Invalid parameter value"
instead — intermittently, crashing the run. Stripped of `CNPIX1`/`CNPIX2`
the same header parsed cleanly 300 times out of 300, so any wcslib error
on a header carrying them now gets the retry.

For a frame with no WCS at all, see [`plate_solve`](@ref) — `run_pipeline`
uses it as a fallback when given `plate_solve_api_key`.
"""
function load_wcs(header::AbstractString)
    solutions = try
        WCS.from_header(String(header))
    catch e
        (e isa ErrorException && _has_fits_card(header, ("CNPIX1", "CNPIX2"))) || rethrow()
        WCS.from_header(_strip_fits_cards(String(header), ("CNPIX1", "CNPIX2")))
    end
    isempty(solutions) && error("no WCS solution found in header")
    return first(solutions)
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

"""
    frame_epoch(header::FITSHeader; timestamp_key=nothing, exptime_key="EXPTIME") -> Float64

A frame's observation epoch as a UTC Julian Date at **mid-exposure** — what
MPC reports require, and what Astrometrica writes into IASC reports.

`timestamp_key` names the header keyword holding the time: a number is read
as a Modified Julian Date, a string as an ISO-8601 date-time (fractional
seconds kept). Left `nothing`, it is `DATE-OBS` when the header has one, else
`MJD-OBS`. That default matters for Pan-STARRS1: its `MJD-OBS` is the shutter
time in **TAI**, not UTC — 37 s later than `DATE-OBS` ("UTC start of
exposure") in every real IASC header checked, 2019 and 2025 alike, despite
the same header's `TIMESYS = 'UTC'` (its `SHUTOPEN` card, labelled TAI, is
the same instant as `MJD-OBS`). Reading `MJD-OBS` put every epoch 37 s late —
caught by comparing this pipeline's MPC report for a real object with
Astrometrica's. ZTF's `OBSMJD` is UTC; pass it explicitly as
`timestamp_key="OBSMJD"` there.

If the header has `exptime_key` (seconds), half of it is added: both surveys
stamp the *start* of the exposure. `exptime_key=nothing` skips that, for a
timestamp that already marks mid-exposure.
"""
function frame_epoch(header::FITSHeader; timestamp_key::Union{Nothing,AbstractString}=nothing,
                     exptime_key::Union{Nothing,AbstractString}="EXPTIME")
    key = timestamp_key === nothing ? (haskey(header, "DATE-OBS") ? "DATE-OBS" : "MJD-OBS") : timestamp_key
    value = header[key]
    jd = if value isa AbstractString
        m = match(r"^(\d{4}-\d\d-\d\d(?:T\d\d:\d\d:\d\d)?)(\.\d+)?", strip(value))
        m === nothing && throw(ArgumentError("$key = $(repr(value)) is not an ISO-8601 date-time"))
        datetime2julian(DateTime(m[1])) + (m[2] === nothing ? 0.0 : parse(Float64, m[2])) / 86400
    else
        Float64(value) + 2400000.5
    end
    if exptime_key !== nothing && haskey(header, exptime_key)
        jd += header[exptime_key] / 2 / 86400
    end
    return jd
end

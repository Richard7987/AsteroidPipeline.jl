"""
    gaia_reference_stars(ra, dec, radius_deg; epoch_jd::Real, mag_range=(12.0, 21.0)) -> Table

Gaia DR3 stars within `radius_deg` of `(ra, dec)` (degrees), with G
magnitude inside `mag_range` — a live network request, to the ESA Gaia
archive's TAP service (`gaiadr3.gaia_source`), falling back to VizieR's
copy (`I/355/gaiadr3`) if the archive fails with a network error, a
timeout or a 5xx status. Two services because one did go down mid-run:
VizieR's TAP answered 503 to every request for a sustained stretch while
this was being validated, which would otherwise have stopped every
`refine_astrometry` run cold.

Positions are propagated by proper motion from Gaia DR3's own epoch
(J2016.0) to `epoch_jd` (a Julian Date, normally the frame's own); stars
without a proper motion keep their catalogued position. A decade of proper motion is easily a
significant fraction of an arcsecond for nearby stars, which is the
whole error budget [`refine_wcs`](@ref) is trying to beat.

Returns a table with columns `ra`, `dec` (degrees, at `epoch_jd`) and
`gmag`, ready to pass to [`refine_wcs`](@ref).
"""
function gaia_reference_stars(ra::Real, dec::Real, radius_deg::Real;
                               epoch_jd::Real, mag_range=(12.0, 21.0))
    0 < radius_deg <= 5 || throw(ArgumentError("radius_deg must be in (0, 5]"))
    lo, hi = mag_range
    # Both queries alias their columns to the same names, so one parser
    # serves either service.
    cone(racol, deccol) = "CONTAINS(POINT('ICRS', $racol, $deccol), " *
                          "CIRCLE('ICRS', $ra, $dec, $radius_deg)) = 1"
    esa = "SELECT ra AS ra_deg, dec AS dec_deg, pmra AS pm_ra, pmdec AS pm_dec, " *
          "phot_g_mean_mag AS g FROM gaiadr3.gaia_source " *
          "WHERE $(cone("ra", "dec")) AND phot_g_mean_mag BETWEEN $lo AND $hi"
    vizier = "SELECT RA_ICRS AS ra_deg, DE_ICRS AS dec_deg, pmRA AS pm_ra, pmDE AS pm_dec, " *
             "Gmag AS g FROM \"I/355/gaiadr3\" " *
             "WHERE $(cone("RA_ICRS", "DE_ICRS")) AND Gmag BETWEEN $lo AND $hi"
    rows = try
        _tap_query(_GAIA_ARCHIVE_TAP_URL, esa)
    catch e
        _is_service_failure(e) || rethrow()
        @warn "Gaia archive unavailable; falling back to VizieR" exception=e
        _tap_query(_VIZIER_TAP_URL, vizier)
    end

    years = (epoch_jd - _GAIA_DR3_EPOCH_JD) / 365.25
    out_ra = Float64[]; out_dec = Float64[]; gmag = Float64[]
    for r in rows
        (ismissing(r.ra_deg) || ismissing(r.dec_deg)) && continue
        pmra = ismissing(r.pm_ra) ? 0.0 : Float64(r.pm_ra)   # mas/yr, already * cos(dec)
        pmde = ismissing(r.pm_dec) ? 0.0 : Float64(r.pm_dec)
        d = Float64(r.dec_deg) + pmde * years / 3.6e6
        push!(out_ra, Float64(r.ra_deg) + pmra * years / 3.6e6 / cosd(d))
        push!(out_dec, d)
        push!(gmag, ismissing(r.g) ? NaN : Float64(r.g))
    end
    return Table(ra=out_ra, dec=out_dec, gmag=gmag)
end

const _GAIA_DR3_EPOCH_JD = 2457389.0   # J2016.0
const _GAIA_ARCHIVE_TAP_URL = "https://gea.esac.esa.int/tap-server/tap/sync"

"""
    refine_wcs(image, wcs, stars; threshold=10.0, mask=nothing, max_offset=60.0,
               match_radius=2.0, min_matches=6, clip_sigma=3.0)
        -> (wcs, refined, n_matches, rms_arcsec, offset_px)

Re-derive a frame's astrometric solution from reference stars actually
detected in it, rather than trusting the header's — the same job
Astrometrica's "Data Reduction" does before any measurement is reported.

`image` is in the `(y, x)` order [`detect_sources`](@ref) takes; `wcs`
is the header's approximate solution; `stars` is a table of reference
stars with `ra`, `dec` columns (degrees, at the frame's epoch), e.g. from
[`gaia_reference_stars`](@ref). Stars are detected at `threshold` sigma
(`mask` is forwarded to `detect_sources`), and every catalog star is
projected through `wcs`. The header's error is found first as one
global pixel offset — the peak of a vote over every detection/star pair
within `max_offset` pixels — then each star is matched to the nearest
detection within `match_radius` pixels of its offset-corrected position,
and a linear tangent-plane (TAN) solution — 6 parameters: shift, scale,
rotation, shear — is fitted by least squares to those matches, rejecting
matches beyond `clip_sigma` times the residual RMS and refitting until
none are rejected.

Why this is needed rather than optional for real IASC/Pan-STARRS1 data:
on a real 2025 practice set (field `XY54_p10`), matched against 38 Gaia
stars, the header WCS was off by a near-constant ~26.6 px (6.8") — far
beyond the <1" an MPC report needs, though PS1's own `CERROR` keyword
claims 0.06". Refit this way, the residual RMS fell to ~0.17", and an
unknown object's positions agreed with Astrometrica's own
Gaia-calibrated measurement of it. The header's `PCA*` polynomial
distortion terms were checked as the explanation first and ruled out:
evaluated at the chip's corners they move a position by <0.5 px.

Returns a named tuple: the refined `wcs` (a plain TAN `WCSTransform`,
tangent point at `wcs`'s own `crval`), `refined` (`false` if fewer than
`min_matches` stars matched — `wcs` is then the input, unchanged — so a
caller can warn instead of failing outright), `n_matches`, `rms_arcsec`
(the fit's residual RMS; `NaN` if not refined), and `offset_px` (the
global `(dx, dy)` offset found, detected minus predicted).
"""
function refine_wcs(image::AbstractMatrix{<:Real}, wcs::WCSTransform, stars;
                    threshold::Real=10.0, mask::Union{Nothing,AbstractMatrix{Bool}}=nothing,
                    max_offset::Real=60.0, match_radius::Real=2.0,
                    min_matches::Integer=6, clip_sigma::Real=3.0)
    unrefined(n, offset) = (wcs=wcs, refined=false, n_matches=n, rms_arcsec=NaN, offset_px=offset)

    detections = detect_sources(image; threshold, mask)
    det_xy = [(d.xcen, d.ycen) for d in detections]
    n1, n2 = size(image)

    predicted = Tuple{Float64,Float64,Float64,Float64}[]   # (x, y, ra, dec)
    for s in stars
        x, y = world_to_pix(wcs, Float64[s.ra, s.dec])
        (-max_offset < x < n2 + max_offset && -max_offset < y < n1 + max_offset) || continue
        push!(predicted, (x, y, Float64(s.ra), Float64(s.dec)))
    end
    (isempty(det_xy) || isempty(predicted)) && return unrefined(0, (NaN, NaN))

    offset = _vote_offset(det_xy, predicted, max_offset)
    matches = _match_stars(det_xy, predicted, offset, match_radius)
    length(matches) < min_matches && return unrefined(length(matches), offset)

    crval = (wcs.crval[1], wcs.crval[2])
    fit = _fit_tan(matches, crval, clip_sigma)
    fit.n < min_matches && return unrefined(fit.n, offset)
    return (wcs=fit.wcs, refined=true, n_matches=fit.n, rms_arcsec=fit.rms_arcsec, offset_px=offset)
end

"""
    _vote_offset(det_xy, predicted, max_offset; bin=2.0) -> (dx, dy)

The global pixel offset (detected minus predicted) most detection/star
pairs agree on: every pair within `max_offset` votes for its own offset
in `bin`-pixel cells, and the winning cell's pairs are averaged. Random
pairings spread their votes thinly over the whole `max_offset` square;
true pairings all land in one cell, so this finds the header's error
without needing it to already be small.
"""
function _vote_offset(det_xy, predicted, max_offset::Real; bin::Real=2.0)
    votes = Dict{Tuple{Int,Int},Vector{Tuple{Float64,Float64}}}()
    for (px, py, _, _) in predicted, (dx_, dy_) in det_xy
        ox, oy = dx_ - px, dy_ - py
        (abs(ox) <= max_offset && abs(oy) <= max_offset) || continue
        push!(get!(votes, (floor(Int, ox / bin), floor(Int, oy / bin)), Tuple{Float64,Float64}[]), (ox, oy))
    end
    isempty(votes) && return (0.0, 0.0)
    # Count each cell together with its 8 neighbours, so a true offset
    # sitting on a cell boundary isn't split across two cells.
    neighbourhood(c) = (get(votes, (c[1] + a, c[2] + b), ()) for a in -1:1, b in -1:1)
    best = argmax(c -> sum(length, neighbourhood(c)), collect(keys(votes)))
    pairs = collect(Iterators.flatten(neighbourhood(best)))
    return (median(first.(pairs)), median(last.(pairs)))
end

function _match_stars(det_xy, predicted, offset, match_radius::Real)
    matches = Tuple{Float64,Float64,Float64,Float64}[]   # (x, y, ra, dec)
    used = Set{Int}()
    for (px, py, ra, dec) in predicted
        tx, ty = px + offset[1], py + offset[2]
        best, best_dist = 0, Float64(match_radius)
        for (i, (x, y)) in enumerate(det_xy)
            d = hypot(x - tx, y - ty)
            d <= best_dist && ((best, best_dist) = (i, d))
        end
        (best == 0 || best in used) && continue
        push!(used, best)
        push!(matches, (det_xy[best]..., ra, dec))
    end
    return matches
end

"""
    _gnomonic(ra, dec, crval) -> (xi, eta)

Standard (tangent-plane) coordinates of `(ra, dec)` about tangent point
`crval`, in degrees — FITS's intermediate world coordinates for a TAN
projection, so a linear fit of pixels to these *is* a TAN WCS solution.
"""
function _gnomonic(ra::Real, dec::Real, crval)
    a0, d0 = crval
    cosc = sind(d0) * sind(dec) + cosd(d0) * cosd(dec) * cosd(ra - a0)
    xi = cosd(dec) * sind(ra - a0) / cosc
    eta = (cosd(d0) * sind(dec) - sind(d0) * cosd(dec) * cosd(ra - a0)) / cosc
    return rad2deg(xi), rad2deg(eta)
end

function _fit_tan(matches, crval, clip_sigma::Real)
    keep = trues(length(matches))
    local a, b, rms
    while true
        m = matches[keep]
        X = hcat(ones(length(m)), [p[1] for p in m], [p[2] for p in m])
        std_coords = [_gnomonic(p[3], p[4], crval) for p in m]
        a = X \ first.(std_coords)
        b = X \ last.(std_coords)
        all_X = hcat(ones(length(matches)), [p[1] for p in matches], [p[2] for p in matches])
        all_std = [_gnomonic(p[3], p[4], crval) for p in matches]
        resid = hypot.(all_X * a .- first.(all_std), all_X * b .- last.(all_std)) .* 3600
        rms = sqrt(mean(resid[keep] .^ 2))
        new_keep = resid .<= max(clip_sigma * rms, 1e-3)
        (new_keep == keep || count(new_keep) < 3) && break
        keep = new_keep
    end
    cd = [a[2] a[3]; b[2] b[3]]
    crpix = -(cd \ [a[1], b[1]])
    wcs = WCSTransform(2; ctype=["RA---TAN", "DEC--TAN"], crval=[crval[1], crval[2]],
                       crpix=crpix, cdelt=[1.0, 1.0], pc=cd)
    return (wcs=wcs, n=count(keep), rms_arcsec=rms)
end

"""
    _field_reference_stars(image, wcs, epoch_jd) -> Table

[`gaia_reference_stars`](@ref) covering all of `image`'s footprint
(through `wcs`), with a 10% margin for the header WCS's own error.
"""
function _field_reference_stars(image::AbstractMatrix, wcs::WCSTransform, epoch_jd::Real)
    n1, n2 = size(image)
    centre = pix_to_world(wcs, Float64[n2 / 2, n1 / 2])
    corners = (pix_to_world(wcs, Float64[x, y]) for (x, y) in ((1, 1), (n2, 1), (1, n1), (n2, n1)))
    radius = 1.1 * maximum(_angular_distance_arcsec(centre..., c...) for c in corners) / 3600
    return gaia_reference_stars(centre[1], centre[2], radius; epoch_jd)
end

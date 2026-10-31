"""
    gaia_reference_stars(ra, dec, radius_deg; epoch_jd::Real, mag_range=(12.0, 21.0),
                         cache=gaia_cache_dir()) -> Table

Gaia DR3 stars within `radius_deg` of `(ra, dec)` (degrees), with G
magnitude inside `mag_range`, from the first of three copies of the
catalogue that answers: the ESA Gaia archive (`gaiadr3.gaia_source`),
VizieR (`I/355/gaiadr3`, Strasbourg) and NOIRLab's Astro Data Lab
(`gaia_dr3.gaia_source`, Tucson). One is tried after another on a network
error, a timeout or a 5xx status: the first two, both in Europe, were down
together for days in October 2026.

Gaia DR3 doesn't change, so each answer is kept in `cache` (a directory;
`nothing` to keep none) and the same query is answered from it after,
without the network: a field analyzed again, or observed again, costs no
query.

Positions are propagated by proper motion from Gaia DR3's own epoch
(J2016.0) to `epoch_jd` (a Julian Date, normally the frame's own); stars
without a proper motion keep their catalogued position. A decade of proper motion is easily a
significant fraction of an arcsecond for nearby stars, which is the
whole error budget [`refine_wcs`](@ref) is trying to beat.

Returns a table with columns `ra`, `dec` (degrees, at `epoch_jd`) and
`gmag`, ready to pass to [`refine_wcs`](@ref).
"""
function gaia_reference_stars(ra::Real, dec::Real, radius_deg::Real;
                               epoch_jd::Real, mag_range=(12.0, 21.0),
                               cache::Union{Nothing,AbstractString}=gaia_cache_dir())
    0 < radius_deg <= 5 || throw(ArgumentError("radius_deg must be in (0, 5]"))
    lo, hi = mag_range
    # Every query aliases its columns to the same names, so one parser serves
    # any service.
    cone(racol, deccol) = "CONTAINS(POINT('ICRS', $racol, $deccol), " *
                          "CIRCLE('ICRS', $ra, $dec, $radius_deg)) = 1"
    esa = "SELECT ra AS ra_deg, dec AS dec_deg, pmra AS pm_ra, pmdec AS pm_dec, " *
          "phot_g_mean_mag AS g FROM gaiadr3.gaia_source " *
          "WHERE $(cone("ra", "dec")) AND phot_g_mean_mag BETWEEN $lo AND $hi"
    vizier = "SELECT RA_ICRS AS ra_deg, DE_ICRS AS dec_deg, pmRA AS pm_ra, pmDE AS pm_dec, " *
             "Gmag AS g FROM \"I/355/gaiadr3\" " *
             "WHERE $(cone("RA_ICRS", "DE_ICRS")) AND Gmag BETWEEN $lo AND $hi"
    # Data Lab's ADQL has no CONTAINS on this table; q3c is its own cone search
    noirlab = "SELECT ra AS ra_deg, dec AS dec_deg, pmra AS pm_ra, pmdec AS pm_dec, " *
              "phot_g_mean_mag AS g FROM gaia_dr3.gaia_source " *
              "WHERE q3c_radial_query(ra, dec, $ra, $dec, $radius_deg) = 't' " *
              "AND phot_g_mean_mag BETWEEN $lo AND $hi"
    services = ((_GAIA_ARCHIVE_TAP_URL, esa, "the Gaia archive"), (_VIZIER_TAP_URL, vizier, "VizieR"),
                (_NOIRLAB_TAP_URL, noirlab, "NOIRLab's Data Lab"))

    cached = cache === nothing ? nothing : joinpath(cache, string(hash(esa); base=16) * ".csv")
    text = cached !== nothing && isfile(cached) ? read(cached, String) : nothing
    if text === nothing
        for (k, (url, query, name)) in enumerate(services)
            try
                text = _tap_csv(url, query)
                break
            catch e
                (_is_service_failure(e) && k < length(services)) || rethrow()
                @warn "$name unavailable; trying $(services[k+1][3])" exception = e
            end
        end
        if cached !== nothing
            mkpath(cache)
            tmp = tempname(cache)
            write(tmp, text)
            mv(tmp, cached; force=true)
        end
    end
    rows = CSV.File(IOBuffer(text))

    years = (epoch_jd - _GAIA_DR3_EPOCH_JD) / 365.25
    out_ra = Float64[]; out_dec = Float64[]; gmag = Float64[]
    for r in rows
        (ismissing(r.ra_deg) || ismissing(r.dec_deg)) && continue
        # missing, or NaN as some services write it: no proper motion
        pm(v) = ismissing(v) || !isfinite(v) ? 0.0 : Float64(v)
        pmra = pm(r.pm_ra)   # mas/yr, already * cos(dec)
        pmde = pm(r.pm_dec)
        d = Float64(r.dec_deg) + pmde * years / 3.6e6
        push!(out_ra, Float64(r.ra_deg) + pmra * years / 3.6e6 / cosd(d))
        push!(out_dec, d)
        push!(gmag, ismissing(r.g) ? NaN : Float64(r.g))
    end
    return Table(ra=out_ra, dec=out_dec, gmag=gmag)
end

"""
    gaia_cache_dir() -> String

Where [`gaia_reference_stars`](@ref) keeps the answers it got:
`\$XDG_CACHE_HOME` (default `~/.cache`) `/AsteroidPipeline/gaia`. Deleting
it only costs the queries again.
"""
gaia_cache_dir() = joinpath(get(ENV, "XDG_CACHE_HOME", joinpath(homedir(), ".cache")), "AsteroidPipeline", "gaia")

const _GAIA_DR3_EPOCH_JD = 2457389.0   # J2016.0
const _GAIA_ARCHIVE_TAP_URL = "https://gea.esac.esa.int/tap-server/tap/sync"
const _NOIRLAB_TAP_URL = "https://datalab.noirlab.edu/tap/sync"

"""
    refine_wcs(image, wcs, stars; threshold=10.0, mask=nothing, max_offset=60.0,
               match_radius=2.0, min_matches=6, clip_sigma=3.0, min_significance=8.0)
        -> (wcs, refined, n_matches, rms_arcsec, offset_px, n_detections, matches)

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

If the offset's vote stands less than `min_significance` robust standard
deviations above the other offsets', the frame is left unrefined: a
damaged frame (smeared, or full of detector defects) gives thousands of
detections, and random pairs then make a false offset win, with hundreds
of matches whose RMS is that of chance, `match_radius / √2`.

Pan-STARRS1 headers can be off by several arcseconds, although their
`CERROR` claims 0.06", far beyond the <1" an MPC report needs; refitted,
the residual RMS is ~0.2".

Returns a named tuple: the refined `wcs` (a plain TAN `WCSTransform`,
tangent point at `wcs`'s own `crval`), `refined` (`false` if fewer than
`min_matches` stars matched — `wcs` is then the input, unchanged — so a
caller can warn instead of failing outright), `n_matches`, `rms_arcsec`
(the fit's residual RMS; `NaN` if not refined), and `offset_px` (the
global `(dx, dy)` offset found, detected minus predicted), and `n_detections`
(sources detected at `threshold` — compared with the catalog's star count it
flags frames full of detector defects; see [`run_pipeline`](@ref)), and
`matches` (the matched stars, `(x, y, ra, dec)`, before the fit's clipping —
what [`refine_wcs_set`](@ref) pools across frames).
"""
function refine_wcs(image::AbstractMatrix{<:Real}, wcs::WCSTransform, stars;
                    threshold::Real=10.0, mask::Union{Nothing,AbstractMatrix{Bool}}=nothing,
                    max_offset::Real=60.0, match_radius::Real=2.0,
                    min_matches::Integer=6, clip_sigma::Real=3.0, min_significance::Real=8.0)
    detections = detect_sources(image; threshold, mask)
    unrefined(n, offset, matches=_NO_MATCHES) = (wcs=wcs, refined=false, n_matches=n, rms_arcsec=NaN,
                                                 offset_px=offset, n_detections=length(detections), matches)
    det_xy = [(d.xcen, d.ycen) for d in detections]
    n1, n2 = size(image)

    predicted = Tuple{Float64,Float64,Float64,Float64}[]   # (x, y, ra, dec)
    for s in stars
        x, y = world_to_pix(wcs, Float64[s.ra, s.dec])
        (-max_offset < x < n2 + max_offset && -max_offset < y < n1 + max_offset) || continue
        push!(predicted, (x, y, Float64(s.ra), Float64(s.dec)))
    end
    (isempty(det_xy) || isempty(predicted)) && return unrefined(0, (NaN, NaN))

    offset, significance = _vote_offset(det_xy, predicted, max_offset)
    significance < min_significance && return unrefined(0, offset)
    matches = _match_stars(det_xy, predicted, offset, match_radius)
    length(matches) < min_matches && return unrefined(length(matches), offset, matches)

    crval = (wcs.crval[1], wcs.crval[2])
    fit = _fit_tan(matches, crval, clip_sigma)
    fit.n < min_matches && return unrefined(fit.n, offset, matches)
    return (wcs=fit.wcs, refined=true, n_matches=fit.n, rms_arcsec=fit.rms_arcsec, offset_px=offset,
            n_detections=length(detections), matches)
end

const _NO_MATCHES = Tuple{Float64,Float64,Float64,Float64}[]

"""
    refine_wcs_set(images, wcss, stars; masks=nothing, wide_offset_arcsec=77.0,
                   max_rms=1.0, iterations=5, rematch_arcsec=2.0, kwargs...) -> Vector

[`refine_wcs`](@ref) for the frames of one field taken in sequence, fitted
together: one result per frame, as `refine_wcs`'s, plus `joint` (whether
the frame's solution came from the joint fit) and `shift_only` (whether only
the header's offset was corrected — see below).

Fitted one frame at a time, a solution rests on the few reference stars
matched in that frame (6 to 16 on IASC's Pan-STARRS1 chips), and frames
can disagree by several pixels, which makes every star seem to move. The
headers' *relative* alignment is reliable even when their absolute
position is 30″ off, so each frame's matches are carried to every other
frame's pixels through the headers, and each frame is fitted to all of
them. The joint solutions are kept if the headers bear that out: fitted
again with an offset of their own, no frame's stars may need one larger
than 3 standard errors plus 0.25 px, and the RMS may not exceed `max_rms`
arcseconds.

Each frame is first matched as [`refine_wcs`](@ref) does (`kwargs` are
passed on); one that fails is retried searching `wide_offset_arcsec`, as
some headers are further off than `refine_wcs`'s default search. Then the
stars are matched again from the fitted solution (within `rematch_arcsec`),
up to `iterations` times while more of them match: the header's error
isn't one offset across a whole chip, so the first match finds only part
of the stars. If the joint fit is rejected, a frame whose own fit failed
but with 4 or more stars agreeing on the header's offset gets that offset
corrected (`shift_only`), which keeps it aligned with the others.
"""
function refine_wcs_set(images::AbstractVector, wcss::AbstractVector, stars;
                        masks=nothing, wide_offset_arcsec::Real=77.0, max_rms::Real=1.0,
                        iterations::Integer=5, rematch_arcsec::Real=2.0,
                        clip_sigma::Real=3.0, min_matches::Integer=6, kwargs...)
    n = length(images)
    length(wcss) == n || throw(ArgumentError("images and wcss must have the same length"))
    masks === nothing && (masks = fill(nothing, n))
    # tolerances are in arcseconds; each frame's search is in its own pixels
    scale = [_pixel_scale_arcsec(w, size(img)) for (w, img) in zip(wcss, images)]
    wide_px = wide_offset_arcsec ./ scale
    rematch_px = rematch_arcsec ./ scale
    # Matched again from each fit until no more stars match: one offset for the
    # whole chip matches only the part where the header's error is that offset.
    function _rematched(r, k)
        for _ in 1:iterations
            r.refined || break
            again = refine_wcs(images[k], r.wcs, stars; mask=masks[k], clip_sigma, min_matches, kwargs...,
                               max_offset=rematch_px[k], min_significance=0)
            again.refined && again.n_matches > r.n_matches || break
            # the offset stays the header's
            r = merge(again, (offset_px=r.offset_px,))
        end
        return r
    end
    tag(r; joint=false, shift_only=false) = merge(r, (; joint, shift_only))
    single = map(1:n) do k
        r = refine_wcs(images[k], wcss[k], stars; mask=masks[k], clip_sigma, min_matches, kwargs...)
        if !r.refined
            wide = refine_wcs(images[k], wcss[k], stars; mask=masks[k], clip_sigma, min_matches, kwargs...,
                              max_offset=wide_px[k])
            (wide.refined || wide.n_matches > r.n_matches) && (r = wide)
        end
        _rematched(r, k)
    end
    # A frame left with far fewer stars than the best one starts again from the
    # best one's solution, carried over through the headers.
    best = argmax(k -> single[k].refined ? single[k].n_matches : -1, 1:n)
    if single[best].refined
        for k in 1:n
            (k == best || single[k].n_matches >= 0.6 * single[best].n_matches) && continue
            guess = _carried_wcs(wcss[k], wcss[best], single[best].wcs, size(images[k]))
            r = refine_wcs(images[k], guess, stars; mask=masks[k], clip_sigma, min_matches, kwargs...,
                           max_offset=rematch_px[k], min_significance=0)
            r.refined || continue
            r = _rematched(merge(r, (offset_px=single[k].offset_px,)), k)
            r.n_matches > single[k].n_matches && (single[k] = r)
        end
    end

    if n >= 2
        joint = map(1:n) do j
            pooled = Tuple{Float64,Float64,Float64,Float64}[]
            source = Int[]
            for k in 1:n, (x, y, ra, dec) in single[k].matches
                p = k == j ? (x, y) : Tuple(world_to_pix(wcss[j], pix_to_world(wcss[k], Float64[x, y])))
                push!(pooled, (p[1], p[2], ra, dec)); push!(source, k)
            end
            length(pooled) < min_matches && return nothing
            f = _fit_tan(pooled, (wcss[j].crval[1], wcss[j].crval[2]), clip_sigma)
            (f.n >= min_matches && f.rms_arcsec <= max_rms && _frames_agree(pooled, source, (wcss[j].crval[1], wcss[j].crval[2]))) ? f : nothing
        end
        if all(!isnothing, joint)
            return [tag((wcs=f.wcs, refined=true, n_matches=f.n, rms_arcsec=f.rms_arcsec, offset_px=r.offset_px,
                         n_detections=r.n_detections, matches=r.matches); joint=true)
                    for (f, r) in zip(joint, single)]
        end
    end
    return map(1:n) do k
        r = single[k]
        (r.refined || r.n_matches < 4) && return tag(r)
        tag(merge(r, (wcs=_shifted_wcs(wcss[k], r.offset_px, size(images[k])),)); shift_only=true)
    end
end

# Whether the frames' headers agree with each other: the pooled stars are
# fitted again with an offset of its own for each frame's stars (relative to
# the first frame's), and no offset may exceed 3 standard errors plus 0.25 px.
# Mean residuals per frame would miss it: the joint fit absorbs much of an offset.
function _frames_agree(pooled, source, crval)
    frames = unique(source)
    length(frames) < 2 && return true
    std_coords = [_gnomonic(p[3], p[4], crval) .* 3600 for p in pooled]       # arcsec
    X = hcat(ones(length(pooled)), [p[1] for p in pooled], [p[2] for p in pooled],
             [Float64(s == k) for s in source, k in frames[2:end]])
    size(X, 1) <= size(X, 2) + 1 && return true
    a = X \ first.(std_coords)
    b = X \ last.(std_coords)
    resid2 = (X * a .- first.(std_coords)) .^ 2 .+ (X * b .- last.(std_coords)) .^ 2
    sigma = sqrt(sum(resid2) / (2 * (size(X, 1) - size(X, 2))))              # per coordinate
    pixel = sqrt(abs(a[2] * b[3] - a[3] * b[2]))                             # arcsec per pixel
    n1 = count(==(frames[1]), source)
    return all(enumerate(frames[2:end])) do (i, k)
        se = sigma * sqrt(1 / count(==(k), source) + 1 / n1)
        hypot(a[3+i], b[3+i]) <= 3 * sqrt(2) * se + 0.25 * pixel
    end
end

# A solution for the frame of header `wcs`, from another frame's refitted
# solution `solved` (header `solved_header`): a grid of this frame's pixels,
# carried to the other frame's through the two headers, and fitted as TAN.
function _carried_wcs(wcs::WCSTransform, solved_header::WCSTransform, solved::WCSTransform, image_size)
    n1, n2 = image_size
    points = map(vec([(x, y) for x in range(1, n2; length=7), y in range(1, n1; length=7)])) do (x, y)
        ra, dec = pix_to_world(solved, world_to_pix(solved_header, pix_to_world(wcs, Float64[x, y])))
        (x, y, ra, dec)
    end
    return _fit_tan(points, (wcs.crval[1], wcs.crval[2]), Inf).wcs
end

# `wcs` with every position moved by `offset` pixels, as a TAN solution:
# fitted to a grid of the header's own positions, shifted.
function _shifted_wcs(wcs::WCSTransform, offset, image_size)
    n1, n2 = image_size
    grid = [(x, y) for x in range(1, n2; length=7), y in range(1, n1; length=7)]
    points = map(vec(grid)) do (x, y)
        ra, dec = pix_to_world(wcs, Float64[x, y])
        (x + offset[1], y + offset[2], ra, dec)
    end
    return _fit_tan(points, (wcs.crval[1], wcs.crval[2]), Inf).wcs
end

"""
    _vote_offset(det_xy, predicted, max_offset; bin=2.0) -> ((dx, dy), significance)

The global pixel offset (detected minus predicted) most detection/star
pairs agree on: every pair within `max_offset` votes for its own offset
in `bin`-pixel cells, and the winning cell's pairs are averaged. Random
pairings spread their votes thinly over the whole `max_offset` square;
true pairings all land in one cell, so this finds the header's error
without needing it to already be small.

Also returns the winning cell's significance: its excess over the median
cell, in robust standard deviations of the cells (at least the median's
Poisson noise). A frame with thousands of spurious detections fills every
cell, and its winner is no higher than chance makes the highest of a few
thousand cells, about 4. `Inf` when most cells are empty: random pairs
can't fill a winner then, however few stars a frame holds.
"""
function _vote_offset(det_xy, predicted, max_offset::Real; bin::Real=2.0)
    votes = Dict{Tuple{Int,Int},Vector{Tuple{Float64,Float64}}}()
    for (px, py, _, _) in predicted, (dx_, dy_) in det_xy
        ox, oy = dx_ - px, dy_ - py
        (abs(ox) <= max_offset && abs(oy) <= max_offset) || continue
        push!(get!(votes, (floor(Int, ox / bin), floor(Int, oy / bin)), Tuple{Float64,Float64}[]), (ox, oy))
    end
    isempty(votes) && return (0.0, 0.0), 0.0
    # Count each cell together with its 8 neighbours, so a true offset
    # sitting on a cell boundary isn't split across two cells.
    neighbourhood(c) = (get(votes, (c[1] + a, c[2] + b), ()) for a in -1:1, b in -1:1)
    count_at(c) = sum(length, neighbourhood(c))
    best = argmax(count_at, collect(keys(votes)))
    pairs = collect(Iterators.flatten(neighbourhood(best)))
    offset = (median(first.(pairs)), median(last.(pairs)))
    lo, hi = floor(Int, -max_offset / bin), floor(Int, max_offset / bin)
    counts = [count_at((a, b)) for a in lo:hi, b in lo:hi]
    m = median(counts)
    σ = max(1.4826 * median(abs.(counts .- m)), sqrt(m))
    return offset, σ == 0 ? Inf : (count_at(best) - m) / σ
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
        # Monotone: a rejected match never comes back, so each pass drops at
        # least one or stops. Re-admitting them could oscillate forever on
        # near-noiseless data, where residuals sit right at the threshold.
        new_keep = keep .& (resid .<= max(clip_sigma * rms, 1e-3))
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

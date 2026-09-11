"""
    link_across_nights(night_movers::AbstractVector;
                        match_radius_arcsec::Real=5.0,
                        max_accel_arcsec_per_day2::Real=Inf,
                        min_nights::Integer=2)

Group tracklets from *different* observing nights that are consistent
with one real moving object, the same problem [`link_candidates`](@ref)
solves within a single night, one level up: on sky coordinates rather
than pixels (different nights can have entirely different pointings/WCS
solutions), and over gaps of days rather than minutes, where the target's
own real orbital curvature — not just its instantaneous rate — starts to
matter.

`night_movers[k]` is night `k`'s own [`astrometric_calibrate`](@ref) (or
[`search_field`](@ref)) output table (`id`, `frame`, `x`, `y`, `ra`,
`dec`, `epoch`; `x`/`y` are unused here). Each night's tracklets (per
`id`) are independently fit to a linear rate in a local tangent-plane
projection around their own mean position — via [`_linfit`](@ref),
reused as-is from [`link_candidates`](@ref)'s own within-night fit, not
reimplemented — then that rate is used to extrapolate forward to every
other night's tracklets' own mean epoch. A cross-night pair is accepted
whenever the extrapolated position lands within `match_radius_arcsec` of
the other tracklet's actual mean position, loosened (if
`max_accel_arcsec_per_day2` is finite) by an added
`max_accel_arcsec_per_day2 * Δt_days^2` allowance for real curvature
over the longer gap. Accepted pairs are grouped transitively (one
tracklet can plausibly pair with several other nights' tracklets, which
must then all collapse into a single group) via union-find; only groups
spanning at least `min_nights` distinct nights are kept. Tracklets with
fewer than 2 points (no rate to fit) are skipped — they carry no motion
information to extrapolate from and so cannot be matched here.

`max_accel_arcsec_per_day2` defaults to `Inf` (off), mirroring
[`link_candidates`](@ref)'s own `max_speed::Real=Inf` default exactly:
a single generous `match_radius_arcsec` does the real matching work by
default. **Neither default is derived from this project's own measured
data** — unlike e.g. `variability_chi2`'s `systematic_error_fraction`
(calibrated against three real confirmed variables), no real multi-night
same-object case was available to calibrate against when this was
written; treat `match_radius_arcsec=5.0` as a starting point to tune
against your own campaign's real cadence and target population, not a
validated constant.

Returns a `Vector` of groups, each a `Vector` of `(night, id)`
`NamedTuple`s — the tracklets, across nights, judged to be the same
object — mirroring [`link_candidates`](@ref)'s own "vector of vectors"
return shape rather than a flat table.
"""
function link_across_nights(night_movers::AbstractVector;
                             match_radius_arcsec::Real=5.0,
                             max_accel_arcsec_per_day2::Real=Inf,
                             min_nights::Integer=2)
    Node = NamedTuple{(:night, :id),Tuple{Int,Int}}
    nodes = Node[]
    fits = NamedTuple[]  # (ra0, dec0, epoch0, rate_ra, rate_dec) per node, same order as nodes

    for (night, candidates) in enumerate(night_movers)
        for id in unique(candidates.id)
            rows = collect(filter(r -> r.id == id, candidates))
            length(rows) < 2 && continue

            ra0, dec0 = mean(r.ra for r in rows), mean(r.dec for r in rows)
            epoch0 = mean(r.epoch for r in rows)
            cos_dec0 = cosd(dec0)
            t = [r.epoch - epoch0 for r in rows]
            dra = [(r.ra - ra0) * cos_dec0 * 3600 for r in rows]
            ddec = [(r.dec - dec0) * 3600 for r in rows]
            ra_int, rate_ra = _linfit(t, dra)
            dec_int, rate_dec = _linfit(t, ddec)

            push!(nodes, (night=night, id=id))
            push!(fits, (ra0=ra0 + ra_int / cos_dec0 / 3600, dec0=dec0 + dec_int / 3600,
                          epoch0=epoch0, rate_ra=rate_ra, rate_dec=rate_dec, cos_dec0=cos_dec0))
        end
    end

    n = length(nodes)
    parent = collect(1:n)
    find(i) = (while parent[i] != i; i = parent[i]; end; i)
    function union!(i, j)
        ri, rj = find(i), find(j)
        ri != rj && (parent[ri] = rj)
    end

    for i in 1:n, j in 1:n
        i == j && continue
        nodes[i].night == nodes[j].night && continue
        a, b = fits[i], fits[j]
        dt = b.epoch0 - a.epoch0
        dt <= 0 && continue  # only extrapolate forward in time; the (j,i) pair covers the reverse

        pred_ra = a.ra0 + (a.rate_ra * dt) / a.cos_dec0 / 3600
        pred_dec = a.dec0 + a.rate_dec * dt / 3600
        sep_arcsec = hypot((pred_ra - b.ra0) * a.cos_dec0, pred_dec - b.dec0) * 3600

        tolerance = match_radius_arcsec
        isfinite(max_accel_arcsec_per_day2) && (tolerance += max_accel_arcsec_per_day2 * dt^2)

        sep_arcsec <= tolerance && union!(i, j)
    end

    groups = Dict{Int,Vector{Node}}()
    for i in 1:n
        push!(get!(groups, find(i), Node[]), nodes[i])
    end

    return [g for g in values(groups) if length(unique(node.night for node in g)) >= min_nights]
end

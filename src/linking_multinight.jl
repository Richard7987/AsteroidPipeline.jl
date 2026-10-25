"""
    link_across_nights(night_movers::AbstractVector;
                        match_radius_arcsec::Real=5.0,
                        max_accel_arcsec_per_day2::Real=Inf,
                        min_nights::Integer=2)

Group tracklets from *different* observing nights that are consistent
with one moving object: [`link_candidates`](@ref)'s problem one level up,
on sky coordinates rather than pixels (each night has its own pointing),
and over gaps of days, where orbital curvature starts to matter.

`night_movers[k]` is night `k`'s [`astrometric_calibrate`](@ref) (or
[`search_field`](@ref)) table (`id`, `frame`, `x`, `y`, `ra`, `dec`,
`epoch`). Each night's tracklets are fit to a linear rate in a tangent
plane around their mean position ([`_linfit`](@ref)) and extrapolated to
every other night's tracklets' mean epoch. A pair is accepted when the
extrapolation lands within `match_radius_arcsec`, plus
`max_accel_arcsec_per_day2 * Δt_days^2` for curvature if that is finite.
Accepted pairs are grouped transitively (union-find), and groups spanning
at least `min_nights` nights are kept. Tracklets with fewer than 2 points
have no rate and are skipped.

The defaults are not calibrated on real multi-night data: tune
`match_radius_arcsec` to your cadence and targets.

Returns groups, each a `Vector` of `(night, id)` named tuples.
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

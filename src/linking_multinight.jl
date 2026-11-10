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

# Earth's mean motion, rad/day.
const _EARTH_MEAN_MOTION = 2π / 365.256

"""
    link_nights(night_movers, observatory; sigma=0.5, max_days=15, outlier_sigma=3.0, min_nights=3)

Link tracklets seen on different nights into one object each, by orbit:
[`link_across_nights`](@ref)'s problem, solved with [`fit_orbit`](@ref)
instead of straight lines.

`night_movers[k]` is night `k`'s table of tracklet points (`id`, `ra`,
`dec` in degrees, `epoch` as JD UTC), as [`astrometric_calibrate`](@ref)
gives them; `observatory` is the MPC site code, one for all nights or one
per night. Tracklets with fewer than 2 points are skipped. The observers'
positions come from JPL Horizons (a live request).

Two tracklets `max_days` apart or less are first compared by their motion:
the mean angular velocity between them must equal the mean of their own
velocities, as it does for any motion of constant acceleration. The
allowance is 3σ of the measured velocities (each point `sigma` arcsec),
plus the next term of the motion, bounded by `6|v|ω²Δt²/12` with ω the
tracklet's own angular rate plus the Earth's mean motion. The diurnal
parallax is taken out first, for a grid of distances from 0.02 to 50 AU;
a pair passes at any one of them. A pair that passes is then fitted
with an orbit, started by ranging over those distances only; it is linked
when every observation is within `outlier_sigma` times `sigma` of it.
Linked pairs are grown, best fitting first, by the tracklet from another
night whose addition fits best, refitted from the orbit already found,
while the orbit still fits every observation. Each tracklet ends in at
most one linkage: those spanning more nights first, then those with
smaller residuals.

Returns the linkages spanning at least `min_nights` nights. Two nights are
not enough: two tracklets moving like asteroids often fit some orbit, and
false tracklets were linked on two nights by the dozen, never on three
([Linking across nights](../validation/linking.md)). Each linkage is a
named tuple: `tracklets` (`(night, id)` named tuples, in time order) and
`orbit` (the [`fit_orbit`](@ref) result).
"""
function link_nights(night_movers::AbstractVector, observatory; sigma::Real=0.5, max_days::Real=15,
                     outlier_sigma::Real=3.0, min_nights::Integer=3)
    codes = observatory isa AbstractString ? fill(String(observatory), length(night_movers)) : String.(collect(observatory))
    length(codes) == length(night_movers) || throw(ArgumentError("one observatory code, or one per night"))
    tracklets = NamedTuple[]
    for (night, rows) in enumerate(night_movers), id in unique(rows.id)
        pts = sort(collect(filter(r -> r.id == id, rows)); by=r -> r.epoch)
        length(pts) < 2 && continue
        push!(tracklets, (night=night, id=id, code=codes[night], ra=Float64[p.ra for p in pts],
                          dec=Float64[p.dec for p in pts], t=[p.epoch + _tt_minus_utc(p.epoch) for p in pts]))
    end
    R = Dict{Tuple{String,Float64},Vector{Float64}}()
    for code in unique(codes), center in (code, nothing)
        t = unique(reduce(vcat, [k.t for k in tracklets if k.code == code]; init=Float64[]))
        isempty(t) && continue
        for (ti, p) in zip(t, _observer_heliocentric_tt(t, center))
            R[(center === nothing ? "500" : code, ti)] = collect(p)
        end
    end
    full = [merge(k, (R=[R[(k.code, ti)] for ti in k.t], G=[R[("500", ti)] for ti in k.t])) for k in tracklets]
    return _link_nights(full, Float64(sigma), Float64(max_days), Float64(outlier_sigma), min_nights)
end

# Each tracklet's geocentric mean direction and angular velocity (unit
# vectors, rad/day) at its mean epoch, the observer's offset from the
# geocentre taken out for an object `Δ` AU away, and their 1σ.
function _geocentric_motion(k, Δ, σ_rad)
    u = [normalize(k.R[i] .+ Δ .* _line_of_sight(k.ra[i], k.dec[i]) .- k.G[i]) for i in eachindex(k.t)]
    t̄ = mean(k.t)
    dt = k.t .- t̄
    S = sum(abs2, dt)
    u0 = mean(u)
    v = S == 0 ? zeros(3) : sum(dt .* u) ./ S
    return (u=u0, v=v, t=t̄, σu=σ_rad / sqrt(length(dt)), σv=S == 0 ? Inf : σ_rad / sqrt(S))
end

function _motion_consistent(a, b)
    Δt = b.t - a.t
    chord = (b.u .- a.u) ./ Δt
    v̄ = (a.v .+ b.v) ./ 2
    ω = norm(v̄) + _EARTH_MEAN_MOTION
    allowed = 3 * sqrt((a.σu^2 + b.σu^2) / Δt^2 + (a.σv^2 + b.σv^2) / 4) + 6 * norm(v̄) * ω^2 * Δt^2 / 12
    return norm(chord .- v̄) <= allowed
end

function _linkage_orbit(group, σ, outlier_sigma; kwargs...)
    ra = reduce(vcat, [k.ra for k in group]); dec = reduce(vcat, [k.dec for k in group])
    t = reduce(vcat, [k.t for k in group]); R = reduce(vcat, [k.R for k in group])
    o = sortperm(t)
    f = try
        _fit_orbit(ra[o], dec[o], t[o], R[o], fill(σ, length(t)), Inf; kwargs...)
    catch
        return nothing
    end
    all(i -> hypot(f.residuals.ra[i], f.residuals.dec[i]) <= outlier_sigma * σ, eachindex(f.residuals.ra)) || return nothing
    return f
end

function _link_nights(tracklets, σ, max_days, outlier_sigma, min_nights)
    n = length(tracklets)
    distances = exp10.(range(log10(0.02), log10(50.0); length=24))
    σ_rad = deg2rad(σ / 3600)
    motion = [_geocentric_motion(k, Δ, σ_rad) for k in tracklets, Δ in distances]
    order = sortperm([mean(k.t) for k in tracklets])
    pairs = Dict{Tuple{Int,Int},Any}()
    for (x, i) in enumerate(order)
        for j in order[x+1:end]
            tracklets[i].night == tracklets[j].night && continue
            mean(tracklets[j].t) - mean(tracklets[i].t) > max_days && break
            ok = [_motion_consistent(motion[i, d], motion[j, d]) for d in eachindex(distances)]
            any(ok) || continue
            # the orbit starts from the distances the motion allows, and those next to them
            near = [any(ok[max(1, d - 1):min(end, d + 1)]) for d in eachindex(distances)]
            f = _linkage_orbit(tracklets[[i, j]], σ, outlier_sigma; distances=distances[near])
            f === nothing || (pairs[(i, j)] = f)
        end
    end

    neighbours = [Int[] for _ in 1:n]
    for (i, j) in keys(pairs)
        push!(neighbours[i], j); push!(neighbours[j], i)
    end
    linkages = Dict{Vector{Int},Any}()
    grown = Set{Tuple{Int,Int}}()
    for ((i, j), f) in sort(collect(pairs); by=p -> p[2].rms)
        (i, j) in grown && continue
        members = [i, j]
        while true
            linkages[sort(members)] = f
            nights = Set(tracklets[m].night for m in members)
            best = nothing
            for c in unique(reduce(vcat, neighbours[members]))
                (c in members || tracklets[c].night in nights) && continue
                g = _linkage_orbit(tracklets[vcat(members, c)], σ, outlier_sigma; from=(f.state, f.epoch))
                g === nothing || (best !== nothing && g.rms >= best[2].rms) || (best = (c, g))
            end
            best === nothing && break
            members, f = vcat(members, best[1]), best[2]
        end
        union!(grown, ((a, b) for a in members, b in members))
    end

    taken = falses(n)
    out = NamedTuple[]
    for (members, f) in sort(collect(linkages); by=p -> (-length(p[1]), p[2].rms))
        any(taken[members]) && continue
        length(members) < min_nights && continue
        taken[members] .= true
        sorted = sort(members; by=m -> mean(tracklets[m].t))
        push!(out, (tracklets=[(night=tracklets[m].night, id=tracklets[m].id) for m in sorted], orbit=f))
    end
    return out
end

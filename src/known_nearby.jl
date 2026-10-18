#=
Is a candidate really new? Identification ([`analyze_set`](@ref)) only
accepts a known object within 5" (SkyBoT) or 15" (MPC orbits) of a
candidate in three frames. An object whose catalogued orbit is poor can be
much further off — on a real 2026-10-04 set SkyBoT placed 2019 UA185 far from
where it was — and is then reported as a new candidate. What tells it apart
is its motion: a known object predicted near the candidate, moving at the
same rate in the same direction, is almost certainly that object. This is
the check Astrometrica's known-object overlay and the MPC's MPChecker
support by eye, done for every candidate.
=#

"""
    tracklet_motion(points) -> (rate_arcsec_h, pa_deg, rms_arcsec)

Straight-line motion through a tracklet's `points` (rows with `ra`, `dec`
in degrees and `epoch` in days): its rate in arcsec/hour, its position angle
in degrees east of north, and the RMS distance of the points from the fitted
line in arcsec — how well they line up, as a real object's do (IASC's
"true signature": a straight line, a constant rate).
"""
function tracklet_motion(points)
    pts = sort(collect(points); by=p -> p.epoch)
    ra0, dec0, t0 = pts[1].ra, pts[1].dec, pts[1].epoch
    ξ = [mod(p.ra - ra0 + 180, 360) - 180 for p in pts] .* cosd(dec0) .* 3600
    η = [(p.dec - dec0) * 3600 for p in pts]
    t = [(p.epoch - t0) * 24 for p in pts]
    length(pts) >= 2 && t[end] > 0 || return (rate_arcsec_h=0.0, pa_deg=0.0, rms_arcsec=0.0)
    X = hcat(ones(length(t)), t)
    a, b = X \ ξ, X \ η
    rms = sqrt(mean((X * a .- ξ) .^ 2 .+ (X * b .- η) .^ 2))
    return (rate_arcsec_h=hypot(a[2], b[2]), pa_deg=mod(atand(a[2], b[2]), 360), rms_arcsec=rms)
end

# Whether a known object predicted `offset` (arcsec east, north: measured
# minus predicted) away from a tracklet is that tracklet with a poor orbit:
# the same rate (10%, or 2"/h for slow objects) and direction (5°), and the
# offset mostly *along* the track — an orbit's error is mostly in time, the
# object running early or late on its own path. Looser tests matched
# unrelated main-belt asteroids: a real candidate (XY54_p10's NHU0001)
# had two within 8' moving within 25% and 20° of it, offset across the track.
function _same_motion(m, pred, offset)
    abs(m.rate_arcsec_h - pred.rate_arcsec_h) <= max(0.10 * m.rate_arcsec_h, 2.0) || return false
    abs(mod(m.pa_deg - pred.pa_deg + 180, 360) - 180) <= 5 || return false
    across = abs(offset[1] * cosd(m.pa_deg) - offset[2] * sind(m.pa_deg))
    return across <= max(30.0, 0.15 * hypot(offset...))
end

"""
    nearby_known_objects(candidates; observatory, mpcorb=nothing, radius_arcmin=10,
                         match_arcsec=15, skybot=true) -> (checks, warnings)

For every tracklet in `candidates` (an [`astrometric_calibrate`](@ref)
table, or its rows), the known minor planets predicted within
`radius_arcmin` of it at its own epochs — from `mpcorb` (as
[`read_mpcorb`](@ref) returns; observer positions from JPL Horizons) and
SkyBoT — compared with it.

`checks` maps each tracklet id to `(verdict, motion, objects)`: `motion` is
its own [`tracklet_motion`](@ref); `objects` are named tuples `name`,
`source`, `separation_arcsec` (median over the frames), `mag` (predicted
V), the object's predicted `rate_arcsec_h` and `pa_deg`, and `agreement`:
`:position` (within `match_arcsec`), `:motion` (further off, but at the same
rate and in the same direction, displaced along its track — the signature of
a poorly known orbit) or `:none`. `verdict` is
`:known` if any object agrees in position, `:possibly_known` if one agrees in
motion, else `:new`. A service that can't be reached is left out, with a
line in `warnings`.
"""
function nearby_known_objects(candidates; observatory::AbstractString, mpcorb=nothing, radius_arcmin::Real=10,
                              match_arcsec::Real=15, skybot::Bool=true)
    rows = collect(candidates)
    warnings = String[]
    ObjT = NamedTuple{(:name, :source, :separation_arcsec, :mag, :rate_arcsec_h, :pa_deg, :agreement),
                      Tuple{String,String,Float64,Float64,Float64,Float64,Symbol}}
    found = Dict{Int,Vector{ObjT}}()
    tracklets = Dict(id => sort(filter(r -> r.id == id, rows); by=r -> r.epoch) for id in unique(r.id for r in rows))
    isempty(tracklets) && return (checks=Dict{Int,Any}(), warnings=warnings)
    motions = Dict(id => tracklet_motion(pts) for (id, pts) in tracklets)
    radius = radius_arcmin * 60

    function consider!(id, name, source, pred)   # pred: (ra, dec, epoch) at the tracklet's epochs
        pts = tracklets[id]
        seps = [_angular_distance_arcsec(p.ra, p.dec, q.ra, q.dec) for (p, q) in zip(pts, pred)]
        sep = median(seps)
        sep <= radius || return
        m = tracklet_motion(pred)
        offset = ((mod(pts[1].ra - pred[1].ra + 180, 360) - 180) * cosd(pts[1].dec) * 3600,
                  (pts[1].dec - pred[1].dec) * 3600)
        agreement = sep <= match_arcsec ? :position : _same_motion(motions[id], m, offset) ? :motion : :none
        push!(get!(found, id, ObjT[]), (name=String(name), source=source, separation_arcsec=sep,
                                        mag=pred[1].mag, rate_arcsec_h=m.rate_arcsec_h, pa_deg=m.pa_deg,
                                        agreement=agreement))
    end

    if mpcorb !== nothing
        try
            orbits = mpcorb isa AbstractString ? read_mpcorb(mpcorb) : mpcorb
            epochs = sort(unique(r.epoch for r in rows))
            observers = Dict(zip(epochs, _observer_heliocentric(epochs, observatory)))
            t0 = epochs[1]
            # objects move far less than a degree during one set (unless very
            # fast, and those are rarely missed by the 15" identification)
            screen = radius / 3600 + 1.0
            near = filter(orbits) do o
                p = _mpcorb_radec(o, t0, observers[t0])
                any(pts -> _angular_distance_arcsec(p.ra, p.dec, pts[1].ra, pts[1].dec) <= screen * 3600,
                    values(tracklets))
            end
            for (id, pts) in tracklets, o in near
                pred = map(pts) do p
                    obs = observers[p.epoch]
                    q = _mpcorb_radec(o, p.epoch, obs)
                    (ra=q.ra, dec=q.dec, epoch=p.epoch, mag=_hg_magnitude(o.H, o.G, q.r, q.Δ, hypot(obs...)))
                end
                consider!(id, o.name, "MPC orbit", pred)
            end
        catch e
            _is_service_failure(e) || e isa HTTP.Exceptions.HTTPError || e isa Base.IOError || rethrow()
            push!(warnings, "MPC orbits unavailable for the nearby-object check: $(sprint(showerror, e))")
        end
    end

    if skybot
        try
            for (id, pts) in tracklets
                have = Set(_bare_name(o.name) for o in get(found, id, ObjT[]))
                first_, last_ = pts[1], pts[end]
                a = Dict(m.name => m for m in _skybot_matches(first_, radius / 3600; observatory))
                b = last_.epoch == first_.epoch ? a :
                    Dict(m.name => m for m in _skybot_matches(last_, radius / 3600; observatory))
                for name in keys(a)
                    (haskey(b, name) && !(_bare_name(name) in have)) || continue
                    # SkyBoT gives the first and last epochs; in between, linear
                    f(p) = (p.epoch - first_.epoch) / max(last_.epoch - first_.epoch, eps())
                    pred = [(ra=a[name].ra + f(p) * (b[name].ra - a[name].ra),
                             dec=a[name].dec + f(p) * (b[name].dec - a[name].dec), epoch=p.epoch, mag=a[name].mv)
                            for p in pts]
                    consider!(id, name, "SkyBoT", pred)
                end
            end
        catch e
            (e isa HTTP.Exceptions.HTTPError || e isa SkyBoTServiceError || e isa Base.IOError) || rethrow()
            push!(warnings, "SkyBoT unavailable for the nearby-object check: $(sprint(showerror, e))")
        end
    end

    checks = Dict{Int,Any}()
    for id in keys(tracklets)
        objects = sort(get(found, id, ObjT[]); by=o -> o.separation_arcsec)
        verdict = any(o -> o.agreement === :position, objects) ? :known :
                  any(o -> o.agreement === :motion, objects) ? :possibly_known : :new
        checks[id] = (verdict=verdict, motion=motions[id], objects=objects)
    end
    return (checks=checks, warnings=warnings)
end

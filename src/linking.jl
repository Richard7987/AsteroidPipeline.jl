"""
    link_candidates(detections_per_frame, timestamps;
                     max_speed::Real=Inf, match_radius::Real=2.0,
                     min_frames::Integer=length(detections_per_frame))

Match source detections across a sequence of frames by consistent linear
motion, producing asteroid candidate tracklets.

`detections_per_frame` is a vector of per-frame detection tables (as
returned by [`detect_sources`](@ref), each with `x`, `y` columns — or,
when present, the sub-pixel `xcen`, `ycen` columns, used instead);
`timestamps` gives the observation time of each frame, in the same time
unit as `max_speed` (pixels per unit time).

Every pair of detections in the first two frames — or, when `min_frames`
is less than the number of frames, in *any* two frames — is treated as a
trial linear motion vector (rejected up front if its implied speed exceeds
`max_speed`), and matched — closest detection within `match_radius`
pixels of the linear extrapolation — against every other frame; this seed
step alone is exact whenever frames 1 and 2's positions already pin down
the true velocity. If that seed tracklet reaches at least 3 frames, its
velocity is then **refit** by least squares over all of its matched
points, and every frame is re-matched against the refined line — this is
what makes the tracklet robust to the seed pair's own position error
(from centroiding noise, say): a frame the noisy 2-point seed missed can
be picked up by the refit, and conversely a frame it only matched by
chance can drop out. If the refit's own velocity exceeds `max_speed`, the
refit is discarded and the tracklet keeps its original seed-built points
rather than being dropped outright (the seed itself already passed the
`max_speed` check). `min_frames` is applied once, after any refit.
Distinct seed pairs whose refit converges to the same final point set
(more than one can lock onto the same real object) are deduplicated.

Seeding from frames 1 and 2 only used to apply regardless of
`min_frames`, so an object missing from either — crossing a detector gap,
masked over a star, or leaving the field — could never be linked, even
when `min_frames` said one missing frame was fine. Measured with
synthetic objects injected into a real IASC set ([`inject_movers`](@ref)):
of fast movers detected cleanly in three of four frames, those missing
frame 1 or 2 were never recovered.

Returns a vector of tracklets; each tracklet is a vector of
`(frame, x, y)` named tuples for detections consistent with the same
linear motion, sorted by frame.
"""
function link_candidates(detections_per_frame, timestamps;
                          max_speed::Real=Inf, match_radius::Real=2.0,
                          min_frames::Integer=length(detections_per_frame))
    nframes = length(detections_per_frame)
    nframes == length(timestamps) ||
        throw(ArgumentError("detections_per_frame and timestamps must have the same length"))

    Point = NamedTuple{(:frame, :x, :y),Tuple{Int,Float64,Float64}}
    tracklets = Vector{Point}[]
    nframes < 2 && return tracklets

    seen = Set{Vector{Point}}()
    # Seeds from frames 1 and 2 alone can only find objects present in both;
    # once min_frames allows a missing frame, any pair has to be able to seed.
    seed_pairs = min_frames >= nframes ? [(1, 2)] : [(a, b) for a in 1:nframes-1 for b in a+1:nframes]

    for (a, b) in seed_pairs
        dt_ab = timestamps[b] - timestamps[a]
        dt_ab == 0 && throw(ArgumentError("timestamps[$a] and timestamps[$b] must differ"))

        for d1 in detections_per_frame[a], d2 in detections_per_frame[b]
            (x1, y1), (x2, y2) = _position(d1), _position(d2)
            dx, dy = x2 - x1, y2 - y1
            speed = hypot(dx, dy) / abs(dt_ab)
            speed > max_speed && continue

            vx, vy = dx / dt_ab, dy / dt_ab
            tracklet = Point[]
            for k in 1:nframes
                if k == a
                    push!(tracklet, (frame=a, x=Float64(x1), y=Float64(y1)))
                elseif k == b
                    push!(tracklet, (frame=b, x=Float64(x2), y=Float64(y2)))
                else
                    dt = timestamps[k] - timestamps[a]
                    best = _closest_detection(detections_per_frame[k], x1 + vx * dt, y1 + vy * dt, match_radius)
                    best === nothing || push!(tracklet, _point(k, best))
                end
            end

            if length(tracklet) >= 3
                refit = _refit_tracklet(tracklet, timestamps, detections_per_frame, match_radius, max_speed)
                refit !== nothing && (tracklet = refit)
            end

            length(tracklet) >= min_frames && (tracklet in seen || (push!(seen, tracklet); push!(tracklets, tracklet)))
        end
    end

    return tracklets
end

"""
    _position(d) -> (x, y)

A detection's position: its sub-pixel centroid (`xcen`, `ycen`) when the
row has one — as [`detect_sources`](@ref)'s do — else its `x`, `y`.
"""
_position(d) = hasproperty(d, :xcen) ? (d.xcen, d.ycen) : (d.x, d.y)

_point(frame::Integer, d) = let (x, y) = _position(d)
    (frame=Int(frame), x=Float64(x), y=Float64(y))
end

"""
    _closest_detection(detections, pred_x, pred_y, match_radius)

The detection in `detections` (a table with `x`, `y` columns) closest to
`(pred_x, pred_y)`, if within `match_radius` pixels; `nothing` otherwise.
"""
function _closest_detection(detections, pred_x::Real, pred_y::Real, match_radius::Real)
    best, best_dist = nothing, match_radius
    for d in detections
        x, y = _position(d)
        dist = hypot(x - pred_x, y - pred_y)
        if dist <= best_dist
            best, best_dist = d, dist
        end
    end
    return best
end

"""
    _linfit(t, y) -> (intercept, slope)

Ordinary least-squares fit of `y = intercept + slope * t`.
"""
function _linfit(t::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(t)
    tbar, ybar = sum(t) / n, sum(y) / n
    num = sum((t[i] - tbar) * (y[i] - ybar) for i in 1:n)
    den = sum((t[i] - tbar)^2 for i in 1:n)
    slope = den == 0 ? 0.0 : num / den
    return ybar - slope * tbar, slope
end

"""
    _refit_tracklet(tracklet, timestamps, detections_per_frame, match_radius, max_speed)

Refit `tracklet`'s velocity by least squares over its own matched points
(relative to `tracklet`'s first frame's timestamp), then re-match every
frame against the refined line. Returns `nothing` if the refit velocity
exceeds `max_speed`, otherwise the refined tracklet (which may have more
or fewer points than the input), sorted by frame.
"""
function _refit_tracklet(tracklet, timestamps, detections_per_frame, match_radius::Real, max_speed::Real)
    t0 = timestamps[tracklet[1].frame]
    ts = [timestamps[p.frame] - t0 for p in tracklet]
    x0, vx = _linfit(ts, [p.x for p in tracklet])
    y0, vy = _linfit(ts, [p.y for p in tracklet])

    hypot(vx, vy) > max_speed && return nothing

    Point = eltype(tracklet)
    refined = Point[]
    for k in eachindex(detections_per_frame)
        dt = timestamps[k] - t0
        pred_x, pred_y = x0 + vx * dt, y0 + vy * dt
        best = _closest_detection(detections_per_frame[k], pred_x, pred_y, match_radius)
        best === nothing && continue
        push!(refined, _point(k, best))
    end

    return refined
end

"""
    _filter_tracklets(tracklets, detections_per_frame, timestamps;
                      min_speed=0.0, max_flux_ratio=Inf, max_residual=Inf) -> Vector

Drop tracklets that fail IASC's own "true signature" tests (its
Signature Guide: a real asteroid moves in a straight line, at constant
speed, at fairly constant brightness — "rejected because the magnitude
fluctuates by more than 1"). Straight line and constant speed are what
[`link_candidates`](@ref) already enforces; this adds the other two:

- `min_speed` (pixels per unit of `timestamps`): the least-squares speed
  over the tracklet's own points must reach it. A tracklet barely moving
  is a static source matched to itself — on a difference image, a
  residual left by a star.
- `max_flux_ratio`: brightest over faintest detection flux must not
  exceed it (2.5 is IASC's 1 magnitude). Chance alignments of noise
  peaks rarely keep a steady flux; a real object does.
- `max_residual` (pixels): for a tracklet missing some frame, the RMS
  distance of its points from their own least-squares line must not
  exceed it. `link_candidates` accepts any point within `match_radius` of
  the line — 7.8 px for IASC sets, wide enough for real motion — so with
  a frame allowed missing, bent chance alignments of three noise peaks get
  through; four on one line by chance are rare, so complete tracklets are
  left alone. On the real IASC sets, partial tracklets of real objects had
  residuals of 0.12-0.31 px against a median 2.6 px for spurious ones, and
  1.0 px cut `XY25_p10` from 744 tracklets to under 50. Complete tracklets
  of faint (V ≈ 21.3) real objects reached 1.40 px — applied to them too,
  the cut lost 2019 PK5 and the Jupiter Trojan 2019 NB9.

With all at their defaults, `tracklets` is returned unchanged.
"""
function _filter_tracklets(tracklets, detections_per_frame, timestamps;
                           min_speed::Real=0.0, max_flux_ratio::Real=Inf, max_residual::Real=Inf)
    (min_speed <= 0 && isinf(max_flux_ratio) && isinf(max_residual)) && return tracklets
    return filter(tracklets) do t
        if min_speed > 0
            ts = [timestamps[p.frame] for p in t]
            _, vx = _linfit(ts, [p.x for p in t])
            _, vy = _linfit(ts, [p.y for p in t])
            hypot(vx, vy) < min_speed && return false
        end
        if isfinite(max_flux_ratio)
            fluxes = [_detection_flux(detections_per_frame[p.frame], p) for p in t]
            all(f -> f > 0, fluxes) || return false
            maximum(fluxes) / minimum(fluxes) > max_flux_ratio && return false
        end
        # only partial tracklets: four hits on one line by chance are rare
        partial = length(t) < length(timestamps)
        partial && isfinite(max_residual) && _line_residual(t, timestamps) > max_residual && return false
        return true
    end
end

# The flux of the detection a tracklet point was taken from (`NaN` if the
# table has no flux column, or no detection at exactly that position).
function _detection_flux(detections, p)
    for d in detections
        _position(d) == (p.x, p.y) && return hasproperty(d, :flux) ? Float64(d.flux) : NaN
    end
    return NaN
end

# RMS distance (pixels) of a tracklet's points from its least-squares
# constant-velocity line.
function _line_residual(tracklet, timestamps)
    ts = [timestamps[p.frame] for p in tracklet]
    x0, vx = _linfit(ts, [p.x for p in tracklet])
    y0, vy = _linfit(ts, [p.y for p in tracklet])
    return sqrt(sum(((p.x - x0 - vx * t)^2 + (p.y - y0 - vy * t)^2) for (p, t) in zip(tracklet, ts)) / length(tracklet))
end

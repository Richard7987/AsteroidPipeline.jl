"""
    link_candidates(detections_per_frame, timestamps;
                     max_speed::Real=Inf, match_radius::Real=2.0,
                     min_frames::Integer=length(detections_per_frame))

Match source detections across a sequence of frames by consistent linear
motion, producing asteroid candidate tracklets.

`detections_per_frame` is a vector of per-frame detection tables (as
returned by [`detect_sources`](@ref), each with `x`, `y` columns);
`timestamps` gives the observation time of each frame, in the same time
unit as `max_speed` (pixels per unit time).

Every pair of detections in the first two frames is treated as a trial
linear motion vector. For frames 3 onward, the closest detection to the
linear extrapolation is kept if it falls within `match_radius` pixels.
Trial pairs whose implied speed exceeds `max_speed`, or whose resulting
tracklet spans fewer than `min_frames` frames, are discarded.

Returns a vector of tracklets; each tracklet is a vector of
`(frame, x, y)` named tuples for detections consistent with the same
linear motion.
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

    dt12 = timestamps[2] - timestamps[1]
    dt12 == 0 && throw(ArgumentError("timestamps[1] and timestamps[2] must differ"))

    for d1 in detections_per_frame[1], d2 in detections_per_frame[2]
        dx, dy = d2.x - d1.x, d2.y - d1.y
        speed = hypot(dx, dy) / abs(dt12)
        speed > max_speed && continue

        vx, vy = dx / dt12, dy / dt12
        tracklet = Point[(frame=1, x=Float64(d1.x), y=Float64(d1.y)),
                          (frame=2, x=Float64(d2.x), y=Float64(d2.y))]

        for k in 3:nframes
            dt = timestamps[k] - timestamps[1]
            pred_x, pred_y = d1.x + vx * dt, d1.y + vy * dt

            best, best_dist = nothing, match_radius
            for d in detections_per_frame[k]
                dist = hypot(d.x - pred_x, d.y - pred_y)
                if dist <= best_dist
                    best, best_dist = d, dist
                end
            end
            best === nothing && continue
            push!(tracklet, (frame=k, x=Float64(best.x), y=Float64(best.y)))
        end

        length(tracklet) >= min_frames && push!(tracklets, tracklet)
    end

    return tracklets
end

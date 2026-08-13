"""
    link_candidates(detections_per_frame, timestamps)

Match source detections across a sequence of frames by consistent linear
motion, producing asteroid candidate tracklets.

`detections_per_frame` is a vector of per-frame detection tables (as
returned by [`detect_sources`](@ref)); `timestamps` gives the observation
time of each frame.
"""
function link_candidates(detections_per_frame, timestamps)
    error("not yet implemented")
end

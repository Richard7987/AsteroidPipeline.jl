#=
Validates search_field/find_variable_sources against a real,
independently-confirmed variable star — the one significant scientific
gap left open by the rest of this session's real-data testing (which
only ever confirmed known *moving* objects via SkyBoT, never a known
*variable*).

Target: ASASSN-V J183620.31 (VSX: type EW — a contact eclipsing binary,
period 0.322 d, amplitude 0.63 mag). Data: ZTF field 487, CCD 12,
quadrant 1, zr filter, night 2019-06-10 — a real high-cadence campaign
with 144 exposures spanning 2.45 h, thinned to every 5th exposure (~29
frames) by fetch_data.sh. Not included in the repo; fetch first with
`fetch_data.sh` in this directory.

Declared before running, not adjusted after seeing results: 2.45 h
covers only ~32% of one 0.322 d period, so full period recovery from
this single night is not expected. But an EW binary varies continuously
(not just at eclipse), so this window should still show clear,
chi2-significant variability if find_variable_sources works — that's the
positive-control criterion this demo actually tests.
=#
using AsteroidPipeline

const DATA_DIR = joinpath(@__DIR__, "..", "data", "real", "variable")
const TARGET_RA = 279.085
const TARGET_DEC = 5.56101

paths = sort(readdir(DATA_DIR, join=true))
isempty(paths) && error("no frames found in $DATA_DIR — run examples/fetch_data.sh first")
println("Running search_field on $(length(paths)) frames...")

# threshold=60 (vs real_data_demo.jl's 8) is deliberate, not a default: this
# field sits near the galactic plane and has ~12,900 detections/frame at
# threshold=8 (vs field 451's ~130), which makes link_candidates's pairwise
# tracklet search (quadratic in detections/frame) impractically slow —
# confirmed directly: it did not finish in 35 min at threshold=8. The
# target itself is extremely bright (flux ~87,000 at this field's noise
# level, found within 1.7 px of its WCS-predicted position at every
# threshold tested from 10 to 100), so this loses none of the signal this
# demo actually needs — a genuine tradeoff for a genuinely dense field, not
# a hidden shortcut.
result = search_field(paths; timestamp_key="OBSMJD", threshold=60.0, match_radius=10.0, max_speed=5000.0)
n_variables = length(unique(result.variables.id))
println("$(n_variables) variable candidate(s) from $(length(result.variables)) detections.")

if n_variables == 0
    println("\nNo variable candidates recovered — find_variable_sources did not flag ",
            "ASASSN-V J183620.31 on this data. Reported as-is.")
else
    first_rows = [first(filter(r -> r.id == id, result.variables)) for id in unique(result.variables.id)]
    matches = crossmatch_catalog(first_rows, :vsx; radius=5.0)
    println("\n$(length(matches))/$(n_variables) candidate(s) match a known VSX variable:")
    for m in matches
        println("  id=$(m.id): $(m.name) ($(m.class), period=$(m.period) d), offset $(round(m.distance_arcsec, digits=1))\"")
    end

    target_distance(row) = 3600 * hypot((row.ra - TARGET_RA) * cosd(TARGET_DEC), row.dec - TARGET_DEC)
    target_row = argmin(target_distance, first_rows)
    if target_distance(target_row) <= 5.0
        println("\nASASSN-V J183620.31 recovered: id=$(target_row.id), ",
                "offset $(round(target_distance(target_row), digits=1))\" from the VSX position.")

        println("\nBuilding a forced-photometry light curve at the recovered position and ",
                "attempting period recovery (not expected to resolve the full 0.322 d period ",
                "from a 2.45 h window — reporting whatever comes out regardless)...")
        times, flux, flux_err = light_curve(paths, target_row.ra, target_row.dec; timestamp_key="OBSMJD")
        period_result = recover_rotation_period(times, flux; minimum_period=0.01, maximum_period=0.5)
        println("  best period: $(round(period_result.period, digits=4)) d, ",
                "power=$(round(period_result.power, digits=2)), ",
                "FAP=$(round(period_result.false_alarm_probability, digits=4))")
        if isapprox(period_result.period, 0.322; atol=0.02)
            println("  matches VSX's catalogued 0.322 d period.")
        else
            println("  does not match VSX's catalogued 0.322 d period (expected, given the ",
                     "partial phase coverage) — FAP above indicates whether *any* periodic ",
                     "signal was detected at all, real or aliased.")
        end
    else
        println("\nASASSN-V J183620.31 itself was NOT among the recovered candidates ",
                "(closest candidate is $(round(target_distance(target_row), digits=1))\" away, ",
                "outside the 5\" match radius) — reported as-is, not treated as a partial success.")
    end
end
